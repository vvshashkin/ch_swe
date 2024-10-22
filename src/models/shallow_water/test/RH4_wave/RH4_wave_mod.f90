module RH4_wave_mod

use domain_mod,               only : domain_t
use domain_factory_mod,       only : create_domain
use stvec_mod,                only : stvec_t
use stvec_swm_mod,            only : stvec_swm_t
use stvec_swm_factory_mod,    only : create_stvec_swm
use operator_mod,             only : operator_t
use operator_swm_factory_mod, only : create_swm_operator
use timescheme_mod,           only : timescheme_t
use timescheme_factory_mod,   only : create_timescheme
use outputer_abstract_mod,    only : outputer_t, outputer_vector_t
use outputer_factory_mod,     only : create_master_paneled_outputer,&
                                     create_latlon_outputer, create_latlon_vec_outputer
use parcomm_mod,              only : parcomm_global
use config_tools_mod,         only : get_config_type_by_extension
use generic_config_mod,       only : generic_config_t

use operator_swm_diff_mod,         only : operator_swm_diff_t
use operator_swm_diff_factory_mod, only : create_swm_diff_operator

use const_mod,  only : pi, Earth_grav, Earth_omega, Earth_radii, &
                       Earth_sidereal_T

use test_fields_mod, only : rh4_wave_height_generator_t, rh4_wave_wind_generator_t
use key_value_mod, only : key_value_r8_t

implicit none

real(kind=8) :: grav   = Earth_grav
real(kind=8) :: omega  = Earth_omega
real(kind=8) :: a      = Earth_radii
real(kind=8) :: h_mean = 8.0_8*10**3

type(rh4_wave_height_generator_t) :: height_field
type(rh4_wave_wind_generator_t)   :: velocity_field

contains

subroutine run_RH4_wave()

    use const_mod,    only : pi
    use vec_math_mod, only : l2norm

    use namelist_read_mod, only : read_namelist_as_str

    type(domain_t)     :: domain

    class(stvec_t),           allocatable :: state, state_err, state_ex
    class(operator_t),        allocatable :: operator
    class(timescheme_t),      allocatable :: timescheme
    class(outputer_t),        allocatable :: outputer
    class(outputer_vector_t), allocatable :: outputer_vec
    class(generic_config_t),  allocatable :: config, subconfig

    type(operator_swm_diff_t),   allocatable :: operator_diff
    class(timescheme_t),         allocatable :: timescheme_diff

    integer(kind=4) :: halo_width = 8

    character(:), allocatable :: config_file, config_string, &
                                 timescheme_name, name, v_components_type
    type(key_value_r8_t) :: diagnostics

    real(kind=8), allocatable :: rotation_axis(:)
    real(kind=8)      :: dt, tau_write, tau_diagnostics, simulation_time
    integer(kind=4)   :: nstep_write, nstep_diagnostics

    real(kind=8)    :: time, l2err, l2_ex
    integer(kind=4) :: it, N
    logical found

    !get and parse config from file
    config_file = "swm_model.cfg"
    call read_namelist_as_str(config_string, config_file, parcomm_global%myid)
    config = get_config_type_by_extension(config_file)
    call config%parse(config_string)

!ensure that correct default values are set:
    call config%get(omega, "domain%metric%omega",found=found,default=Earth_omega)
    if(.not. found) call config%set("domain%metric%omega",Earth_omega)
    if(omega /= Earth_omega) &
        print *, "WARNING: RH4 using non-standard omega = ", omega

    call config%get(a, "domain%metric%planet_radius",found=found,default=Earth_radii)
    if(.not. found) call config%set("domain%metric%planet_radius",Earth_radii)
    if(a /= Earth_radii) &
        print *, "WARNING: RH4 using non-standard planet radius = ", a

    call config%get(rotation_axis,"domain%metric%rotation_axis",found=found, &
                    default = [0.0_8, 0.0_8, 1.0_8])
    if(.not. found) call config%set("domain%metric%rotation_axis",rotation_axis)
    if(rotation_axis(1) /= 0.0_8 .or. rotation_axis(2) /= 0.0_8 .or. &
       rotation_axis(3) /= 1.0_8) &
        print *, "WARNING: RH4 using non-standard rotation axis = ", rotation_axis

    !get time parameters
    call config%get(dt,"dt")
    call config%get(tau_write,"tau_write")
    call config%get(tau_diagnostics,"tau_diagnostics")
    nstep_write = nint(tau_write/dt)
    nstep_diagnostics = nint(tau_diagnostics/dt)
    call config%get(simulation_time,"simulation_time")

    height_field = rh4_wave_height_generator_t(h_mean = h_mean, omega = omega, &
           a = a, grav = grav)

    velocity_field = rh4_wave_wind_generator_t(a = a, omega = omega)

    !write CFL values for diagnostics
    call config%get(N,"domain%N")
    if (parcomm_global%myid==0) print*, "Gravity Wave CFL = ", &
            real(dt*sqrt(h_mean*Earth_grav)    &
            /(2*pi/4/N)/a,4)

    !create model basic entities:
    call config%get(subconfig,"domain")
    call create_domain(domain, subconfig)

    call create_stvec_swm(state,     domain, halo_width, 0)
    call create_stvec_swm(state_ex,  domain, 0         , 0)
    call create_stvec_swm(state_err, domain, 0         , 0)

    call config%get(name,"swm_operator_section")
    call config%get(subconfig,name)
    call subconfig%get(v_components_type,"v_components_type")
    call create_swm_operator(operator, Earth_grav, subconfig, domain)

    call config%get(timescheme_name,"timescheme_name")
    call create_timescheme(timescheme, state, timescheme_name)

    call config%get(name,"diffusion_section")
    call config%get(subconfig,name)
    call config%get(timescheme_name,name//"%diff_time_scheme")
    call create_swm_diff_operator(operator_diff, subconfig, domain)
    call create_timescheme(timescheme_diff, state, timescheme_name)

    print*, 4*domain%partition%Nh, 2*domain%partition%Nh+1

    if(domain%horizontal_staggering == "Ah") then
        call create_latlon_outputer(outputer, 2*domain%partition%Nh+1, 4*domain%partition%Nh, "xy", domain)
        call create_latlon_vec_outputer(outputer_vec,  2*domain%partition%Nh+1, 4*domain%partition%Nh, "xy", "xy", &
                                        v_components_type, domain)
    else if(domain%horizontal_staggering == "C") then
        call create_latlon_outputer(outputer, 2*domain%partition%Nh+1, 4*domain%partition%Nh, "o", domain)
        call create_latlon_vec_outputer(outputer_vec,  2*domain%partition%Nh+1, 4*domain%partition%Nh, "x","y", &
                                       v_components_type, domain)
    else
        call parcomm_global%abort("This staggering is not implemented in"//&
                                  " barotropic instability swe test output:"//&
                                  domain%horizontal_staggering)
    end if


    call get_exact_solution(state,    domain, v_components_type)
    call get_exact_solution(state_ex, domain, v_components_type)

    select type(state)
    class is (stvec_swm_t)
        call outputer%write(state%h, domain, 'h_swm.dat',1)
        call outputer_vec%write(state%u,state%v,domain,"u.dat","v.dat",1)
    end select

    do it = 1, int(simulation_time/dt)

        call timescheme%step(state, operator, domain, dt)
        call timescheme_diff%step(state, operator_diff, domain, dt)

        time = it*dt

        ! call state_err%assign(1.0_8, state_ex, -1.0_8, state, domain)

        if(mod(it, nstep_diagnostics) == 0) then
            diagnostics = operator%get_diagnostics(state, domain)
            if(domain%parcomm%myid==0) call diagnostics%print()
        end if

        if(mod(it,nstep_write) == 0) then

            select type(state)
            class is (stvec_swm_t)
                call outputer%write(state%h, domain, 'h.dat', int(it/nstep_write)+1)
                l2err = l2norm(state%h, domain%mesh_p, domain%parcomm)
                if (parcomm_global%myid==0) print*, "Hours = ", real(time/3600 ,4), &
                                                    "L2err =", real(l2err,4)
            end select
        end if
    end do

end subroutine run_RH4_wave

subroutine get_exact_solution(state, domain, v_components_type)

    use test_fields_mod, only : set_vector_test_field, &
                                set_scalar_test_field

    class(stvec_t), intent(inout) :: state
    type(domain_t), intent(in)    :: domain
    character(*),   intent(in)    :: v_components_type

    select type(state)
    class is (stvec_swm_t)

        call set_scalar_test_field(state%h, height_field, domain%mesh_p, 0)
        call set_vector_test_field(state%u, state%v, velocity_field, &
              domain%mesh_u, domain%mesh_v, 0, v_components_type)

    end select


end subroutine get_exact_solution

end module RH4_wave_mod
