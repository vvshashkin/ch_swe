module oscillating_gaussian_hill_test_mod

use test_fields_mod,     only : scalar_field_generator_t, &
                                set_vector_test_field, set_scalar_test_field, &
                                oscillating_gaussian_hill_scalar_field_generator_t
use domain_mod,          only : domain_t
use generic_config_mod,  only : generic_config_t
use stvec_flexible_mod,  only : stvec_flexible_t, stvec_t
use grid_field_mod,      only : grid_field_t
use abstract_scorer_mod, only : scorer_t

use const_mod,           only : pi, Earth_omega, Earth_radii, &
                                Earth_grav, Earth_sidereal_T

implicit none

contains

subroutine setup_oscillating_Gaussian_hill_test(state, config, v_components_type, domain)

    class(generic_config_t), intent(inout) :: config
    character(len=*),        intent(in)    :: v_components_type
    type(domain_t),          intent(in)    :: domain

    class(stvec_t),               intent(inout) :: state

    type(grid_field_t), pointer :: h, u, v
    type(oscillating_gaussian_hill_scalar_field_generator_t) :: height_field
    integer(kind=4) :: nu

    call config%get(nu,"nu",default=32)
    height_field = oscillating_gaussian_hill_scalar_field_generator_t(nu=nu)

    select type(state)
        class is (stvec_flexible_t)
    
            call state%get_field(h,"h")
            call set_scalar_test_field(h, height_field, domain%mesh_p, 0)
    
            call state%get_field(u,"u")
            call state%get_field(v,"v")
            call u%assign(0.0_8, domain%mesh_u)
            call v%assign(0.0_8, domain%mesh_v)

    end select 

end subroutine setup_oscillating_Gaussian_hill_test

end module oscillating_gaussian_hill_test_mod