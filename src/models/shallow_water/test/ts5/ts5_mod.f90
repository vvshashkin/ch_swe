module ts5_mod


use domain_mod,               only : domain_t
use domain_factory_mod,       only : create_domain
use stvec_mod,                only : stvec_t
use stvec_swm_mod,            only : stvec_swm_t
use stvec_swm_factory_mod,    only : create_stvec_swm
use operator_mod,             only : operator_t
use operator_swm_mod,         only : operator_swm_t
use operator_swm_factory_mod, only : create_swm_operator
use timescheme_mod,           only : timescheme_t
use timescheme_factory_mod,   only : create_timescheme
use outputer_abstract_mod,    only : outputer_t, outputer_vector_t
use outputer_factory_mod,     only : create_master_paneled_outputer,&
                                     create_latlon_outputer, create_latlon_vec_outputer
use parcomm_mod,              only : parcomm_global
use generic_config_mod,       only : generic_config_t
use config_tools_mod,         only : get_config_type_by_extension
use grid_field_mod,           only : grid_field_t
use grid_field_factory_mod,   only : create_grid_field


use operator_swm_diff_mod,         only : operator_swm_diff_t
use operator_swm_diff_factory_mod, only : create_swm_diff_operator

use test_fields_mod, only : solid_rotation_t, ts2_height_generator_t, &
                            ts5_orography_generator_t,                &
                            set_scalar_test_field, set_vector_test_field

use key_value_mod, only : key_value_r8_t
use const_mod,     only : pi, Earth_omega, Earth_radii, &
                          Earth_grav, Earth_sidereal_T

implicit none

real(kind=8), parameter :: grav    = Earth_grav
real(kind=8), parameter :: h_mean  = 5960.0_8
real(kind=8), parameter :: a       = Earth_radii
real(kind=8), parameter :: u0      = 20.0_8
real(kind=8), parameter :: r_mount = pi / 9.0_8
real(kind=8), parameter :: h_mount = 2000.0_8
real(kind=8), parameter :: x_mount = 0.0_8
real(kind=8), parameter :: y_mount = 0.5_8*sqrt(3.0_8)
real(kind=8), parameter :: z_mount = 0.5_8

type(ts2_height_generator_t)    :: height_field
type(ts5_orography_generator_t) :: orography_generator
type(solid_rotation_t)          :: velocity_field

type(grid_field_t)              :: h_orog, h_total

contains

subroutine run_ts5()

    use vec_math_mod, only : l2norm
    use namelist_read_mod, only : read_namelist_as_str

    class(generic_config_t), allocatable :: config, subconfig
    type(domain_t)     :: domain
    type(grid_field_t) :: curl

    class(stvec_t),           allocatable :: state
    class(operator_t),        allocatable :: operator
    class(timescheme_t),      allocatable :: timescheme
    class(outputer_t),        allocatable :: outputer, outputer_curl
    class(outputer_vector_t), allocatable :: outputer_vec

    type(operator_swm_diff_t),   allocatable :: operator_diff
    class(timescheme_t),         allocatable :: timescheme_diff

    integer(kind=4) :: halo_width = 8

    character(:), allocatable :: config_file, config_string, &
                                 timescheme_name, name, v_components_type
    type(key_value_r8_t) :: diagnostics

    real(kind=8)     :: omega, a
    real(kind=8), allocatable :: rotation_axis(:)
    real(kind=8)     :: dt, tau_write, simulation_time, tau_diagnostics
    integer(kind=4)  :: nstep_write, nstep_diagnostics, N
    logical          :: found

    real(kind=8)    :: time, l2err, l2_ex, l2u, l2v
    integer(kind=4) :: it

    !get and parse config from file
    config_file = "swm_model.cfg"
    call read_namelist_as_str(config_string, config_file, parcomm_global%myid)
    config = get_config_type_by_extension(config_file)
    call config%parse(config_string)

    !ensure that correct default values are set:
    call config%get(omega, "domain%metric%omega",found=found,default=Earth_omega)
    if(.not. found) call config%set("domain%metric%omega",Earth_omega)
    if(omega /= Earth_omega) &
        print *, "WARNING: TS5 using non-standard omega = ", omega

    call config%get(a, "domain%metric%planet_radius",found=found,default=Earth_radii)
    if(.not. found) call config%set("domain%metric%planet_radius",Earth_radii)
    if(a /= Earth_radii) &
        print *, "WARNING: TS5 using non-standard planet radius = ", a

    call config%get(rotation_axis,"domain%metric%rotation_axis",found=found, &
                    default = [0.0_8, 0.0_8, 1.0_8])
    if(.not. found) call config%set("domain%metric%rotation_axis",rotation_axis)
    if(rotation_axis(1) /= 0.0_8 .or. rotation_axis(2) /= 0.0_8 .or. &
       rotation_axis(3) /= 1.0_8) &
        print *, "WARNING: TS5 using non-standard rotation axis = ", rotation_axis

    height_field = ts2_height_generator_t(h_mean = h_mean, &
                u0 = u0, omega = omega, a = a, &
                grav = grav)

    orography_generator = ts5_orography_generator_t(r_mount=r_mount,   &
                          h_mount=h_mount, x_mount=x_mount, &
                          y_mount=y_mount, z_mount=z_mount)

    velocity_field = solid_rotation_t(u0 = u0)

    !get time parameters
    call config%get(dt,"dt")
    call config%get(tau_write,"tau_write")
    call config%get(tau_diagnostics,"tau_diagnostics")
    nstep_write = nint(tau_write/dt)
    nstep_diagnostics = nint(tau_diagnostics/dt)
    call config%get(simulation_time,"simulation_time")

    !write CFL values for diagnostics
    call config%get(N,"domain%N")
    if (parcomm_global%myid==0) print*, "Advective CFL = ", &
            real(dt*u0/(2*pi/4/N)/a,4)

    if (parcomm_global%myid==0) print*, "Gravity Wave CFL = ", &
            real(dt*sqrt(h_mean*Earth_grav)    &
            /(2*pi/4/N)/a,4)

    !create model basic entities:
    call config%get(subconfig,"domain")
    call create_domain(domain, subconfig)

    call create_grid_field(h_orog,0,0,domain%mesh_p)
    call create_grid_field(h_total,0,0,domain%mesh_p)

    call create_stvec_swm(state,     domain, halo_width, 0)

    call config%get(name,"swm_operator_section")
    call config%get(subconfig,name)
    call subconfig%get(v_components_type,"v_components_type")
    call create_swm_operator(operator, Earth_grav, subconfig, domain, orography_generator)

    call config%get(timescheme_name,"timescheme_name")
    call create_timescheme(timescheme, state, timescheme_name)

    call config%get(name,"diffusion_section")
    call config%get(subconfig,name)
    call config%get(timescheme_name,name//"%diff_time_scheme")
    call create_swm_diff_operator(operator_diff, subconfig, domain)
    call create_timescheme(timescheme_diff, state, timescheme_name)

    print*, 4*domain%partition%Nh, 2*domain%partition%Nh+1

    if(domain%horizontal_staggering == "Ah") then
        call create_latlon_outputer(outputer, 2*domain%partition%Nh+1, 4*domain%partition%Nh, "Ah", domain)
        call create_latlon_outputer(outputer_curl,     2*domain%partition%Nh+1, 4*domain%partition%Nh, "Ah", domain)
        call create_latlon_vec_outputer(outputer_vec,  2*domain%partition%Nh+1, 4*domain%partition%Nh, "Ah", &
                                        v_components_type, domain)
    else if(domain%horizontal_staggering == "C") then
        call create_latlon_outputer(outputer, 2*domain%partition%Nh+1, 4*domain%partition%Nh, "A", domain)
        call create_latlon_outputer(outputer_curl,     2*domain%partition%Nh+1, 4*domain%partition%Nh, "Ah", domain)
        call create_latlon_vec_outputer(outputer_vec,  2*domain%partition%Nh+1, 4*domain%partition%Nh, "C", &
                                        v_components_type, domain)
    else
        call parcomm_global%abort("This staggering is not implemented in"//&
                                  " swe test 2 output:"//&
                                  domain%horizontal_staggering)
    end if

    call create_grid_field(curl, halo_width, 0, domain%mesh_q)

    select type(state)
    class is (stvec_swm_t)

        call set_scalar_test_field(state%h, height_field, domain%mesh_p, 0)
        call set_scalar_test_field(h_orog, orography_generator, domain%mesh_p, 0)
        call state%h%update(-1.0_8,h_orog,domain%mesh_p)

        call set_vector_test_field(state%u, state%v, velocity_field, &
                                   domain%mesh_u, domain%mesh_v, 0, v_components_type)

        call h_total%assign(1.0_8,state%h,1.0_8,h_orog,domain%mesh_p)
        call outputer%write(h_total, domain, 'h.dat', 1)
        call outputer_vec%write(state%u, state%v, domain, 'u.dat', 'v.dat', 1)
        select type(operator)
        class is (operator_swm_t)
            call operator%curl_op%calc_curl(curl, state%u, state%v, domain)
            call outputer_curl%write(curl, domain, 'curl.dat', 1)
        end select
    end select


    do it = 1, int(simulation_time/dt)

        call timescheme%step(state, operator, domain, dt)
        select type(state)
        class is (stvec_swm_t)
            call state%h%assign(1.0_8,state%h,1.0_8,h_orog,domain%mesh_p)
            call timescheme_diff%step(state, operator_diff, domain, dt)
            call state%h%assign(1.0_8,state%h,-1.0_8,h_orog,domain%mesh_p)
        end select

        time = it*dt

        if(mod(it, nstep_diagnostics) == 0) then
            diagnostics = operator%get_diagnostics(state, domain)
            if(parcomm_global%myid == 0) call diagnostics%print()
        end if

        if(mod(it,nstep_write) == 0) then

            select type(state)
            class is (stvec_swm_t)
                call h_total%assign(1.0_8,state%h,1.0_8,h_orog,domain%mesh_p)
                call outputer%write(h_total, domain, 'h.dat', int(it/nstep_write)+1)
                call outputer_vec%write(state%u, state%v, domain, 'u.dat', 'v.dat', int(it/nstep_write)+1)
                select type(operator)
                class is (operator_swm_t)
                    call operator%curl_op%calc_curl(curl, state%u, state%v, domain)
                    call outputer_curl%write(curl, domain, 'curl.dat', int(it/nstep_write)+1)
                end select
            end select

            if (parcomm_global%myid==0) print*, "Hours = ", real(time/3600 ,4), &
                                                "irec=",int(it/nstep_write)+1

        end if
    end do

end subroutine run_ts5

end module ts5_mod
