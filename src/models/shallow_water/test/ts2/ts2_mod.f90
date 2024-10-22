module ts2_mod

use domain_mod,               only : domain_t
use domain_factory_mod,       only : create_domain
use stvec_mod,                only : stvec_t
use stvec_swm_mod,            only : stvec_swm_t
use stvec_swm_factory_mod,    only : create_stvec_swm
use operator_mod,             only : operator_t
use operator_swm_factory_mod, only : create_swm_operator
use operator_swm_mod,         only : operator_swm_t
use timescheme_mod,           only : timescheme_t
use timescheme_factory_mod,   only : create_timescheme
use outputer_abstract_mod,    only : outputer_t, outputer_vector_t
use outputer_factory_mod,     only : create_master_paneled_outputer,&
                                     create_latlon_outputer, create_latlon_vec_outputer
use parcomm_mod,              only : parcomm_global
use generic_config_mod,       only : generic_config_t
use config_tools_mod,         only : get_config_type_by_extension

use timer_mod, only : init_timer_mod, print_timer, start_timer, stop_timer

use operator_swm_diff_mod,         only : operator_swm_diff_t
use operator_swm_diff_factory_mod, only : create_swm_diff_operator

use test_fields_mod, only : solid_rotation_t, ts2_height_generator_t

use key_value_mod, only : key_value_r8_t
use const_mod,     only : pi, Earth_omega, Earth_radii, &
                          Earth_grav, Earth_sidereal_T

use timer_mod, only : init_timer_mod, stop_timer, start_timer, print_timer

implicit none

real(kind=8), parameter :: h_mean = 29400/Earth_grav
real(kind=8), parameter :: u0 = pi*Earth_radii/6.0_8/Earth_sidereal_T

type(ts2_height_generator_t) :: height_field
type(solid_rotation_t)       :: velocity_field

contains

subroutine run_ts2()

    use vec_math_mod, only : l2norm

    use namelist_read_mod, only : read_namelist_as_str

    class(generic_config_t), allocatable :: config, subconfig
    type(domain_t)     :: domain

    class(stvec_t),           allocatable :: state, state_err, state_ex
    class(operator_t),        allocatable :: operator
    class(timescheme_t),      allocatable :: timescheme, timescheme_diff
    class(outputer_t),        allocatable :: outputer
    class(outputer_vector_t), allocatable :: outputer_vec

    type(operator_swm_diff_t),   allocatable :: operator_diff

    integer(kind=4) :: halo_width = 8

    character(:), allocatable :: config_file, config_string, &
                                 timescheme_name, name, v_components_type
    type(key_value_r8_t) :: diagnostics

    real(kind=8)     :: omega, a
    real(kind=8), allocatable :: rotation_axis(:)
    real(kind=8)     :: dt, tau_write, simulation_time, tau_diagnostics
    integer(kind=4)  :: nstep_write, nstep_diagnostics, N
    logical          :: found

    real(kind=8)    :: time
    real(kind=8)    :: l2_err_h, l2_ex_h, l2_err_u, l2_ex_u
    real(kind=8)    :: linf_err_h, linf_ex_h, linf_err_u, linf_ex_u
    integer(kind=4) :: it

    !get and parse config from file
    config_file = "swm_model.cfg"
    call read_namelist_as_str(config_string, config_file, parcomm_global%myid)
    config = get_config_type_by_extension(config_file)
    call config%parse(config_string)

    !ensure that correct default values are set:
    call config%get(omega, "domain%metric%omega",found=found,default=Earth_omega)
    if(.not. found) call config%set("domain%metric%omega",Earth_omega)
    if(found .and. omega /= Earth_omega) &
        print *, "WARNING: TS2 using non-standard omega = ", omega

    call config%get(a, "domain%metric%planet_radius",found=found,default=Earth_radii)
    if(.not. found) call config%set("domain%metric%planet_radius",Earth_radii)
    if(found .and. a /= Earth_radii) &
        print *, "WARNING: TS2 using non-standard planet radius = ", a

    call config%get(rotation_axis,"domain%metric%rotation_axis",found=found, &
                    default = [0.0_8, 0.0_8, 1.0_8])
    if(.not. found) call config%set("domain%metric%rotation_axis",rotation_axis)

    !tune initial field generators:
    height_field = ts2_height_generator_t(h_mean = h_mean, &
                u0 = u0, omega = omega, a = a, &
                grav = Earth_grav, axis = rotation_axis)

    velocity_field = solid_rotation_t(u0 = u0, axis = rotation_axis)

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

    call init_timer_mod(domain%parcomm%comm_w)

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
    else if(domain%horizontal_staggering == "Ch") then
        call create_latlon_outputer(outputer, 2*domain%partition%Nh+1, 4*domain%partition%Nh, "xy", domain)
        call create_latlon_vec_outputer(outputer_vec,  2*domain%partition%Nh+1, 4*domain%partition%Nh, "y", "x", &
                                   v_components_type, domain)
    else if(domain%horizontal_staggering == "C") then
        call create_latlon_outputer(outputer, 2*domain%partition%Nh+1, 4*domain%partition%Nh, "o", domain)
        call create_latlon_vec_outputer(outputer_vec,  2*domain%partition%Nh+1, 4*domain%partition%Nh, "x", "y", &
                                       v_components_type, domain)
    else
        call parcomm_global%abort("This staggering is not implemented in"//&
                                  " swe test 2 output:"//&
                                  domain%horizontal_staggering)
    end if

    call get_exact_solution(state,    domain, v_components_type)
    call get_exact_solution(state_ex, domain, v_components_type)

    select type(state_ex)
    class is (stvec_swm_t)
        l2_ex_h = l2norm(state_ex%h,  domain%mesh_p, domain%parcomm)
        l2_ex_u = sqrt(l2norm(state_ex%u,  domain%mesh_u, domain%parcomm)**2+&
                       l2norm(state_ex%v,  domain%mesh_v, domain%parcomm)**2)
        linf_ex_h = state_ex%h%maxabs(domain%mesh_p, domain%parcomm)
        linf_ex_u = max(state_ex%u%maxabs(domain%mesh_u, domain%parcomm), &
                        state_ex%v%maxabs(domain%mesh_v, domain%parcomm))
    end select

    select type(state)
    class is (stvec_swm_t)
        call outputer%write(state%h, domain, 'h.dat', 1)
        call outputer_vec%write(state%u, state%v, domain, 'u.dat', 'v.dat', 1)
    end select

    !initialize timer_mod
    call init_timer_mod(domain%parcomm%comm_w)
    call start_timer("main_loop")

    do it = 1, int(simulation_time/dt)

        call start_timer("step")
            call timescheme%step(state, operator, domain, dt)
        call stop_timer("step")

        call start_timer("diffusion")
            ! call timescheme_diff%step(state, operator_diff, domain, dt)
        call stop_timer("diffusion")

        time = it*dt

        call start_timer("IO & diagnoostics")

        call state_err%assign(1.0_8, state_ex, -1.0_8, state, domain)

        if(mod(it, nstep_diagnostics) == 0) then
            diagnostics = operator%get_diagnostics(state, domain)
            if(parcomm_global%myid == 0) call diagnostics%print()
            select type (operator)
            type is(operator_swm_t)
                call operator%apply(state_err,state,domain)
                diagnostics = operator%get_diagnostics_tend(state, state_err, domain)
                if(parcomm_global%myid == 0) call diagnostics%print()
            end select

            select type(state_err)
            class is (stvec_swm_t)
                !call outputer%write(state_err%h, domain, 'h_err.dat', int(it/nstep_write))
                !call outputer_vec%write(state_err%u, state_err%v, domain, 'u_err.dat', 'v_err.dat', int(it/nstep_write))

                l2_err_h = l2norm(state_err%h,  domain%mesh_p, domain%parcomm) / l2_ex_h
                l2_err_u = sqrt(l2norm(state_err%u,  domain%mesh_u, domain%parcomm)**2+&
                                l2norm(state_err%v,  domain%mesh_v, domain%parcomm)**2) / l2_ex_u
                linf_err_h = state_err%h%maxabs(domain%mesh_p, domain%parcomm) / linf_ex_h
                linf_err_u = max(state_err%u%maxabs(domain%mesh_u, domain%parcomm), &
                                 state_err%v%maxabs(domain%mesh_v, domain%parcomm)) / linf_ex_u
                if (parcomm_global%myid==0) print '(A,F12.4,4(A,E15.7))', &
                                                  "Hours = ", real(time/3600 ,4), &
                                                  " l2_h =", real(l2_err_h,4), " linf_h = ", real(linf_err_h,4), &
                                                  " l2_u =", real(l2_err_u,4), " linf_u = ", real(linf_err_u,4)
            end select
        end if

        if(mod(it,nstep_write) == 0) then

            select type(state)
            class is (stvec_swm_t)
                call outputer%write(state%h, domain, 'h.dat', int(it/nstep_write)+1)
                call outputer_vec%write(state%u, state%v, domain, 'u.dat', 'v.dat', int(it/nstep_write)+1)
            end select
        end if

        call stop_timer("IO & diagnoostics")
    end do

    call stop_timer("main_loop")
    call print_timer(net_section_name = "main_loop")

end subroutine run_ts2

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

end module ts2_mod
