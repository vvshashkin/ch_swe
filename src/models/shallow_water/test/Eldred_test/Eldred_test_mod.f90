module Eldred_test_mod

use domain_mod,               only : domain_t
use domain_factory_mod,       only : create_domain
use grid_field_mod,           only : grid_field_t
use grid_field_factory_mod,   only : create_grid_field
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
use generic_config_mod,       only : generic_config_t
use config_tools_mod,         only : get_config_type_by_extension

use operator_swm_diff_mod,         only : operator_swm_diff_t
use operator_swm_diff_factory_mod, only : create_swm_diff_operator
use operator_swm_mod,              only : operator_swm_t

use const_mod,  only : Earth_grav, Earth_omega, Earth_radii, pi, Earth_sidereal_T

use random_friction_mod, only : random_scalar_t, initialize_random_scalar
use test_fields_mod, only : Eldred_test_height_generator, Eldred_test_wind_generator, &
                            set_scalar_test_field, set_vector_test_field
use key_value_mod,   only : key_value_r8_t

implicit none

real(kind=8) :: grav    = Earth_grav
real(kind=8) :: h_mean  = 1e3_8
real(kind=8) :: tau_forcing = 86400._8*100._8

contains

subroutine run_Eldred_test()

    use vec_math_mod, only : l2norm
    use namelist_read_mod, only : read_namelist_as_str

    class(generic_config_t), allocatable :: config, subconfig
    type(domain_t)     :: domain

    class(stvec_t),           allocatable :: state, state_eq
    class(operator_t),        allocatable :: operator
    class(timescheme_t),      allocatable :: timescheme
    class(outputer_t),        allocatable :: outputer, outputer_curl
    class(outputer_vector_t), allocatable :: outputer_vec

    type(random_scalar_t)   :: random_scalar
    type(operator_swm_diff_t),   allocatable :: operator_diff
    class(timescheme_t),         allocatable :: timescheme_diff

    integer(kind=4) :: halo_width = 8

    real(kind=8)     :: omega, a
    real(kind=8), allocatable :: rotation_axis(:)
    character(:), allocatable :: config_file, config_string, &
                                 timescheme_name, name, v_components_type
    type(key_value_r8_t) :: diagnostics

    real(kind=8)     :: dt, tau_write, simulation_time, tau_diagnostics
    integer(kind=4)  :: nstep_write, nstep_diagnostics, N

    real(kind=8)    :: time, l2err, l2_ex
    integer(kind=4) :: it
    logical :: found

    type(grid_field_t) :: curl, ut, vt, div

    !get and parse config from file
    config_file = "swm_model.cfg"
    call read_namelist_as_str(config_string, config_file, parcomm_global%myid)
    config = get_config_type_by_extension(config_file)
    call config%parse(config_string)

    !ensure that correct default values are set:
    call config%get(omega, "domain%metric%omega",found=found,default=Earth_omega)
    if(.not. found) call config%set("domain%metric%omega",Earth_omega)
    if(omega /= Earth_omega) &
        print *, "WARNING: Eldred test using non-standard omega = ", omega

    call config%get(a, "domain%metric%planet_radius",found=found,default=Earth_radii)
    if(.not. found) call config%set("domain%metric%planet_radius",Earth_radii)
    if(a /= Earth_radii) &
        print *, "WARNING: Eldred test using non-standard planet radius = ", a

    call config%get(rotation_axis,"domain%metric%rotation_axis",found=found, &
                    default = [0.0_8, 0.0_8, 1.0_8])
    if(.not. found) call config%set("domain%metric%rotation_axis",rotation_axis)
    if(rotation_axis(1) /= 0.0_8 .or. rotation_axis(2) /= 0.0_8 .or. &
       rotation_axis(3) /= 1.0_8) &
        print *, "WARNING: Eldred test using non-standard rotation axis = ", rotation_axis

    !get time parameters
    call config%get(dt,"dt")
    call config%get(tau_write,"tau_write")
    call config%get(tau_diagnostics,"tau_diagnostics")
    nstep_write = nint(tau_write/dt)
    nstep_diagnostics = nint(tau_diagnostics/dt)
    call config%get(simulation_time,"simulation_time")

    !write CFL values for diagnostics
    call config%get(N,"domain%N")
    if (parcomm_global%myid==0) print*, "Gravity Wave CFL = ", &
            real(dt*sqrt(h_mean*Earth_grav)    &
            /(2*pi/4/N)/a,4)

    !create model basic entities:
    call config%get(subconfig,"domain")
    call create_domain(domain, subconfig)

    call create_stvec_swm(state,     domain, halo_width, 0)
    call create_stvec_swm(state_eq,  domain, halo_width, 0)

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

    call initialize_random_scalar(random_scalar, l=6, tau=1e4_8, &
    amp=10._8, domain=domain)

    print*, 4*domain%partition%Nh, 2*domain%partition%Nh+1

    if(domain%horizontal_staggering == "Ah") then
        call create_latlon_outputer(outputer, 2*domain%partition%Nh+1, 4*domain%partition%Nh, "xy", domain)
        call create_latlon_outputer(outputer_curl, 2*domain%partition%Nh+1, 4*domain%partition%Nh, "xy", domain)
        call create_latlon_vec_outputer(outputer_vec,  2*domain%partition%Nh+1, 4*domain%partition%Nh, "xy", "xy", &
                                   v_components_type, domain)
    else if(domain%horizontal_staggering == "Ch") then
        call create_latlon_outputer(outputer, 2*domain%partition%Nh+1, 4*domain%partition%Nh, "xy", domain)
        call create_latlon_vec_outputer(outputer_vec,  2*domain%partition%Nh+1, 4*domain%partition%Nh, "y", "x", &
                                   v_components_type, domain)
    else if(domain%horizontal_staggering == "C") then
        call create_latlon_outputer(outputer, 2*domain%partition%Nh+1, 4*domain%partition%Nh, "o", domain)
        call create_latlon_outputer(outputer_curl, 2*domain%partition%Nh+1, 4*domain%partition%Nh, "xy", domain)
        call create_latlon_vec_outputer(outputer_vec,  2*domain%partition%Nh+1, 4*domain%partition%Nh, "x", "y", &
                                       v_components_type, domain)
    else
        call parcomm_global%abort("This staggering is not implemented in"//&
                                  " barotropic instability swe test output:"//&
                                  domain%horizontal_staggering)
    end if

    select type(state_eq)
    class is (stvec_swm_t)
        ! call init_equilibrium(state_eq, domain)
        call set_scalar_test_field(state_eq%h,Eldred_test_height_generator,domain%mesh_p,0)
        call set_vector_test_field(state_eq%u,state_eq%v,Eldred_test_wind_generator, &
                                   domain%mesh_u,domain%mesh_v,0,v_components_type)
    end select

    select type(state)
    class is (stvec_swm_t)
        call state%u%assign(0.0_8,domain%mesh_u)
        call state%v%assign(0.0_8,domain%mesh_v)
        !call state%h%assign(h_mean,domain%mesh_p)
        call random_scalar%apply_update_forcing(state%h,domain,dt)
        call state%h%update(h_mean,domain%mesh_p)
        ! call state%assign(1.0_8,state_eq,0.0_8,state_eq,domain)

        call outputer%write(state%h, domain, 'h.dat',1)
        call outputer_vec%write(state%u,state%v,domain,"u.dat","v.dat",1)
    end select

    call create_grid_field(curl, halo_width, 0, domain%mesh_q)
    call create_grid_field(ut, halo_width, 0, domain%mesh_u)
    call create_grid_field(vt, halo_width, 0, domain%mesh_v)
    call create_grid_field(div, halo_width, 0, domain%mesh_p)

    do it = 1, int(simulation_time/dt)

        select type(operator)
        type is (operator_swm_t)
            !call random_scalar%apply_update_forcing(operator%h_surf,domain,dt)
            !if(mod(it,nstep_write) == 0) then
            !	call outputer%write(operator%h_surf, domain, 'h_surf.dat', it/nstep_write)
            !end if
        end select

        call timescheme%step(state, operator, domain, dt)
        ! call timescheme_diff%step(state, operator_diff, domain, dt)

        call state%update(dt/tau_forcing, state_eq, -dt/tau_forcing, state, domain)

        time = it*dt

        if(mod(it, nstep_diagnostics) == 0) then
            diagnostics = operator%get_diagnostics(state, domain)
            if(parcomm_global%myid == 0) call diagnostics%print()
        end if

        if(mod(it,nstep_write) == 0) then


            select type(state)
            class is (stvec_swm_t)
                call outputer%write(state%h, domain, 'h.dat', int(it/nstep_write)+1)
                call outputer_vec%write(state%u,state%v,domain,"u.dat","v.dat",int(it/nstep_write)+1)

                select type(operator)
                class is (operator_swm_t)
                    call operator%curl_op%calc_curl(curl, state%u, state%v, domain)
                    call outputer_curl%write(curl, domain, 'curl.dat', int(it/nstep_write)+1)
                    call operator%co2contra_op%transform(ut,vt,state%u,state%v,domain)
                    call operator%div_op%calc_div(div,ut,vt,domain)
                    call outputer%write(div, domain, 'div.dat', int(it/nstep_write)+1)
                end select

                if (parcomm_global%myid==0) print*, "Hours = ", real(time/3600 ,4), &
                                                    "rec=",int(it/nstep_write)+1
            end select
        end if
    end do

end subroutine run_Eldred_test

end module Eldred_test_mod
