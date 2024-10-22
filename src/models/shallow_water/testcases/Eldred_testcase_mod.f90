module Eldred_testcase_mod

use swm_forcing_mod,       only : swm_forcing_t
use domain_mod,            only : domain_t
use generic_config_mod,    only : generic_config_t
use stvec_flexible_mod,    only : stvec_t, stvec_flexible_t
use stvec_swm_factory_mod, only : create_stvec_swm
use const_mod,             only : Earth_omega, Earth_radii
use grid_field_mod,        only : grid_field_t
use test_fields_mod,       only : scalar_field_generator_t,    &
                                  set_vector_test_field,       &
                                  set_scalar_test_field,       &
                                  Eldred_test_height_generator,&
                                  Eldred_test_wind_generator,  &
                                  zero_scalar_field_generator_t
use random_friction_mod,   only : random_scalar_t, initialize_random_scalar

implicit none

real(kind=8), parameter, private :: forcing_tau_default = 86400._8*100._8
real(kind=8), parameter, private :: h_mean_default = 1e3_8

type, extends(swm_forcing_t) :: Eldred_forcing_t

    class(stvec_t), allocatable :: equilibrium
    real(kind=8) :: tau = forcing_tau_default

    contains
        procedure :: get_tend => forcing_get_tend
        procedure :: apply => forcing_apply

end type

contains

subroutine setup_Eldred_testcase(state, orography_generator, forcing, &
                                 config, v_components_type, domain)

    class(generic_config_t), intent(inout) :: config
    character(len=*),        intent(in)    :: v_components_type
    type(domain_t),          intent(in)    :: domain

    class(stvec_t),                               intent(inout) :: state
    class(scalar_field_generator_t), allocatable, intent(out)   :: orography_generator
    class(swm_forcing_t),            allocatable, intent(out)   :: forcing


    real(kind=8) :: h_mean, omega, a, tau
    type(grid_field_t), pointer :: h, u, v

    type(Eldred_forcing_t), allocatable :: test_forcing
    type(random_scalar_t)               :: random_scalar

    !get parameters
    call config%get(h_mean,"h_mean",default = h_mean_default)
    call config%get(tau,"forcing_tau",default = forcing_tau_default)

    omega = domain%metric%omega
    a = domain%metric%scale

    !ensure that correct default values are set:
    if(domain%parcomm%myid == 0) then

        if(h_mean /= h_mean_default) &
            print *, "WARNING: Eldred test using non-standard h_mean = ", h_mean

        if(omega /= Earth_omega) &
            print *, "WARNING: Eldred test using non-standard omega = ", omega

        if(a /= Earth_radii) &
            print *, "WARNING: Eldred test using non-standard planet radius = ", a

        if(tau /=  forcing_tau_default) &
            print *, "WARNING: Eldred test using non-standard forcing time-scale = ", tau

    end if

    select type(state)
    class is (stvec_flexible_t)

        call state%get_field(h,"h")
        call initialize_random_scalar(random_scalar, l=6, tau=1.0_8, &
                                      amp=10._8, domain=domain)
        call random_scalar%apply_update_forcing(h,domain,dt=1.0_8)
        call h%update(h_mean,domain%mesh_p)

        call state%get_field(u,"u")
        call state%get_field(v,"v")
        call u%assign(0.0_8, domain%mesh_u)
        call v%assign(0.0_8, domain%mesh_v)

    end select

    allocate(test_forcing)

    call create_stvec_swm(test_forcing%equilibrium, domain, 0, 0)

    select type(equilibrium=>test_forcing%equilibrium)
    class is (stvec_flexible_t)

        call equilibrium%get_field(h,"h")
        call equilibrium%get_field(u,"u")
        call equilibrium%get_field(v,"v")

        call set_scalar_test_field(h,  Eldred_test_height_generator,domain%mesh_p,0)
        call set_vector_test_field(u,v,Eldred_test_wind_generator, &
                                   domain%mesh_u,domain%mesh_v,0,v_components_type)


    end select

    call move_alloc(test_forcing, forcing)

    orography_generator = zero_scalar_field_generator_t()

end subroutine

subroutine forcing_get_tend(this, tend, state, domain)

    class(Eldred_forcing_t), intent(inout) :: this

    class(stvec_t),          intent(inout) :: tend

    class(stvec_t),          intent(inout) :: state
    type(domain_t),          intent(in)    :: domain

    call tend%assign(this%tau, this%equilibrium, -this%tau, state, domain)

end subroutine

subroutine forcing_apply(this, state, dt, domain)

    class(Eldred_forcing_t), intent(inout) :: this

    class(stvec_t),          intent(inout) :: state

    real(kind=8),            intent(in)    :: dt
    type(domain_t),          intent(in)    :: domain

    call state%update(dt/this%tau, this%equilibrium, -dt/this%tau, state, domain)

end subroutine


end module