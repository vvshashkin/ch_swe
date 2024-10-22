module RH4_testcase_mod

use test_fields_mod,    only : rh4_wave_height_generator_t,      &
                               rh4_wave_wind_generator_t,        &
                               set_vector_test_field, set_scalar_test_field
use domain_mod,         only : domain_t
use generic_config_mod, only : generic_config_t
use stvec_flexible_mod, only : stvec_flexible_t, stvec_t
use grid_field_mod,     only : grid_field_t

use const_mod,          only : pi, Earth_omega, Earth_radii, &
                               Earth_grav

implicit none

real(kind=8), parameter, private :: h_mean_default = 8.0e3_8

contains

subroutine setup_RH4_test(state, config, v_components_type, domain)

    class(generic_config_t), intent(inout) :: config
    character(len=*),        intent(in)    :: v_components_type
    type(domain_t),          intent(in)    :: domain

    class(stvec_t),          intent(inout) :: state

    !locals
    integer(kind=4) :: N_quad_points
    real(kind=8)    :: omega, a, rotation_axis(3), h_mean
    type(rh4_wave_height_generator_t) :: height_field
    type(rh4_wave_wind_generator_t)   :: velocity_field
    type(grid_field_t), pointer :: h, u, v

    omega = domain%metric%omega
    rotation_axis = domain%metric%rotation_axis
    a = domain%metric%scale

    call config%get(h_mean,"h_mean",default=h_mean_default)

    if(domain%parcomm%myid == 0) then

        if(h_mean /= h_mean_default) &
            print *, "WARNING: Barotropic instability test using non-standard h_mean = ", h_mean

        if(omega /= Earth_omega) &
            print *, "WARNING: Barotropic instability test using non-standard omega = ", omega

        if(a /= Earth_radii) &
            print *, "WARNING: Barotropic instability test using non-standard planet radius = ", a

        if(any(rotation_axis /= [0.0_8,0.0_8,1.0_8])) &
            print *, "WARNING: Barotropic instability test using non-standard Earth rotation axis = ", rotation_axis

    end if

    height_field = rh4_wave_height_generator_t(h_mean = h_mean, omega = omega, &
           a = a, grav = Earth_grav)

    velocity_field = rh4_wave_wind_generator_t(a = a, omega = omega)

    select type(state)
    class is (stvec_flexible_t)

        call state%get_field(h,"h")
        call state%get_field(u,"u")
        call state%get_field(v,"v")

    end select

    call set_scalar_test_field(h, height_field, domain%mesh_p, 0)
    call set_vector_test_field(u, v, velocity_field, &
                               domain%mesh_u, domain%mesh_v, 0, v_components_type)

end subroutine setup_RH4_test


end module RH4_testcase_mod