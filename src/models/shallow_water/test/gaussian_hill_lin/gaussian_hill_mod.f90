module gaussian_hill_mod

use domain_mod,               only : domain_t
use domain_factory_mod,       only : create_domain
use stvec_mod,                only : stvec_t
use stvec_swm_mod,            only : stvec_swm_t
use stvec_swm_factory_mod,    only : create_stvec_swm
use operator_mod,             only : operator_t
use operator_swm_factory_mod, only : create_linear_swm_operator
use timescheme_mod,           only : timescheme_t
use timescheme_factory_mod,   only : create_timescheme
use outputer_abstract_mod,    only : outputer_t, outputer_vector_t
use outputer_factory_mod,     only : create_master_paneled_outputer,&
                                     create_latlon_outputer, create_latlon_vec_outputer
use parcomm_mod,              only : parcomm_global
use test_fields_mod, only : solid_rotation_t, gaussian_hill_scalar_field_generator_t

use operator_swm_factory_mod, only : create_swm_operator

use key_value_mod, only : key_value_r8_t

use const_mod,  only : Earth_grav, Earth_omega, Earth_radii, pi, Earth_sidereal_T

use generic_config_mod, only : generic_config_t
use config_tools_mod,   only : parse_config_file

implicit none

real(kind=8), parameter :: h_mean = 29400/Earth_grav
real(kind=8), parameter :: u0 = pi*Earth_radii/6.0_8/Earth_sidereal_T

type(gaussian_hill_scalar_field_generator_t) :: height_field
type(solid_rotation_t)                       :: velocity_field

contains

subroutine run_gaussian_hill()

    use const_mod,    only : pi
    use vec_math_mod, only : l2norm

    use namelist_read_mod, only : read_namelist_as_str

    type(domain_t)     :: domain

    class(stvec_t),           allocatable :: state
    class(operator_t),        allocatable :: operator
    class(timescheme_t),      allocatable :: timescheme
    class(outputer_t),        allocatable :: outputer
    class(outputer_vector_t), allocatable :: outputer_vec
    class(generic_config_t),  allocatable :: config, subconfig

    integer(kind=4) :: halo_width = 8

    character(:), allocatable :: namelist_string
    type(key_value_r8_t) :: diagnostics

    real(kind=8)     :: dt
    real(kind=8)     :: tau_write, tau_diagnostics, simulation_time
    integer(kind=4)  :: nstep_write, nstep_diagnostics

    real(kind=8)    :: time
    integer(kind=4) :: it

    real(kind=8)              :: omega, a
    real(kind=8), allocatable :: rotation_axis(:)

    character(:), allocatable :: timescheme_name, name, v_components_type

    config = parse_config_file("swm_model.cfg", parcomm_global)

    block
        logical :: found

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

    end block


    height_field = gaussian_hill_scalar_field_generator_t()

    velocity_field = solid_rotation_t(u0 = u0, axis = rotation_axis)

    !get time parameters
    call config%get(dt, "dt")
    call config%get(tau_write, "tau_write")
    call config%get(tau_diagnostics, "tau_diagnostics")
    nstep_write = nint(tau_write/dt)
    nstep_diagnostics = nint(tau_diagnostics/dt)
    call config%get(simulation_time, "simulation_time")

    cfl_diagnostics: block
    !write CFL values for diagnostics

        integer(kind=4) :: N

        call config%get(N, "domain%N")
        if (parcomm_global%myid==0) print*, "Advective CFL = ", &
                real(dt*u0/(2*pi/4/N)/a,4)

        if (parcomm_global%myid==0) print*, "Gravity Wave CFL = ", &
                real(dt*sqrt(h_mean*Earth_grav)    &
                /(2*pi/4/N)/a,4)

    end block cfl_diagnostics

    !create model basic entities:
    call config%get(subconfig,"domain")
    call create_domain(domain, subconfig)

    call create_stvec_swm(state, domain, halo_width, 0)

    call config%get(name, "swm_operator_section")
    call config%get(subconfig, name)
    call subconfig%get(v_components_type, "v_components_type")
    call create_swm_operator(operator, Earth_grav, subconfig, domain)

    call config%get(timescheme_name, "timescheme_name")
    call create_timescheme(timescheme, state, timescheme_name)

    print*, 4*domain%partition%Nh, 2*domain%partition%Nh+1

    if(domain%horizontal_staggering == "Ah") then
        call create_latlon_outputer(outputer, 2*domain%partition%Nh+1, 4*domain%partition%Nh, "xy", domain)
        call create_latlon_vec_outputer(outputer_vec,  2*domain%partition%Nh+1, 4*domain%partition%Nh, "xy", "xy", &
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

    call get_initial_state(state, domain, v_components_type)

    select type(state)
    class is (stvec_swm_t)
        call outputer%write(state%h, domain, 'h.dat', 1)
        call outputer_vec%write(state%u, state%v, domain, 'u.dat', 'v.dat', 1)
    end select


    print*, tau_write

    do it = 1, int(simulation_time/dt)

        call timescheme%step(state, operator, domain, dt)

        time = it*dt

        print*, it

        if(mod(it,nstep_write) == 0) then

            select type(state)
            class is (stvec_swm_t)
                call outputer%write(state%h, domain, 'h.dat', int(it/nstep_write)+1)
                call outputer_vec%write(state%u, state%v, domain, 'u.dat', 'v.dat', int(it/nstep_write)+1)
            end select
        end if
    end do

end subroutine run_gaussian_hill

subroutine get_initial_state(state, domain, v_components_type)

    use test_fields_mod, only : set_vector_test_field, &
                                set_scalar_test_field

    class(stvec_t), intent(inout) :: state
    type(domain_t), intent(in)    :: domain
    character(*),   intent(in)    :: v_components_type

    select type(state)
    class is (stvec_swm_t)

        call set_scalar_test_field(state%h, height_field, domain%mesh_p, 0)
        call state%u%assign(0.0_8, domain%mesh_u)
        call state%v%assign(0.0_8, domain%mesh_v)

    end select


end subroutine get_initial_state

end module gaussian_hill_mod
