module Williamson_test5_mod

use test_fields_mod,        only : scalar_field_generator_t,                     &
                                   set_vector_test_field, set_scalar_test_field, &
                                   solid_rotation_t, ts2_height_generator_t,     &
                                   ts5_orography_generator_t
use domain_mod,             only : domain_t
use generic_config_mod,     only : generic_config_t
use stvec_flexible_mod,     only : stvec_flexible_t, stvec_t
use grid_field_mod,         only : grid_field_t
use grid_field_factory_mod, only : create_grid_field

use const_mod,              only : pi, Earth_omega, Earth_radii, &
                                   Earth_grav
implicit none

real(kind=8), private, parameter :: h_mean_default = 5960.0_8
real(kind=8), private, parameter :: u0_default = 20.0_8
real(kind=8), private, parameter :: r_mount_default = pi / 9.0_8
real(kind=8), private, parameter :: h_mount_default = 2000.0_8
real(kind=8), private, parameter :: x_mount = 0.0_8
real(kind=8), private, parameter :: y_mount = 0.5_8*sqrt(3.0_8)
real(kind=8), private, parameter :: z_mount = 0.5_8


contains

subroutine setup_Williamson_test5(state, orography_generator, config, &
                                  v_components_type, domain)

    class(generic_config_t), intent(inout) :: config
    character(len=*),        intent(in)    :: v_components_type
    type(domain_t),          intent(in)    :: domain

    class(stvec_t),                               intent(inout) :: state
    class(scalar_field_generator_t), allocatable, intent(out)   :: orography_generator


    real(kind=8) :: u0, h_mean, r_mount, h_mount, &
                    omega, a, grav, rotation_axis(3)
    type(grid_field_t), pointer :: h, u, v
    type(grid_field_t) :: h_surf
    type(ts2_height_generator_t) :: height_field
    type(solid_rotation_t)       :: velocity_field

    !get parameters
    call config%get(u0,     "u0",     default = u0_default)
    call config%get(h_mean, "h_mean", default = h_mean_default)
    call config%get(r_mount,"r_mount",default = r_mount_default)
    call config%get(h_mount,"h_mount",default = h_mount_default)

    omega = domain%metric%omega
    rotation_axis = domain%metric%rotation_axis
    a = domain%metric%scale

    !ensure that correct default values are set:
    if(domain%parcomm%myid == 0) then

        if(u0 /= u0_default) &
            print *, "WARNING: Williamson test 5 using non-standard u0 = ", u0

        if(h_mean /= h_mean_default) &
            print *, "WARNING: Williamson test 5 using non-standard h_mean = ", h_mean

        if(omega /= Earth_omega) &
            print *, "WARNING: Williamson test 5 using non-standard omega = ", omega

        if(a /= Earth_radii) &
            print *, "WARNING: Williamson test 5 using non-standard planet radius = ", a

        if(any(rotation_axis /= [0.0_8,0.0_8,1.0_8])) &
            print *, "WARNING: Williamson test 5 using non-standard Earth rotation axis = ", rotation_axis

    end if

    height_field = ts2_height_generator_t(h_mean = h_mean, &
                                          u0 = u0, omega = omega, a = a, &
                                          grav = Earth_grav, axis = rotation_axis)

    velocity_field = solid_rotation_t(u0 = u0, axis = rotation_axis)

    orography_generator = ts5_orography_generator_t(r_mount=r_mount,   &
                              h_mount=h_mount, x_mount=x_mount,        &
                              y_mount=y_mount, z_mount=z_mount)
    
    select type(state)
    class is (stvec_flexible_t)

        call state%get_field(h,"h")

        call state%get_field(u,"u")
        call state%get_field(v,"v")

    end select

    call create_grid_field(h_surf, 0 ,0, domain%mesh_p)
    call set_scalar_test_field(h_surf, orography_generator, domain%mesh_p, 0)
    call set_scalar_test_field(h, height_field, domain%mesh_p, 0)
    call h%update(-1.0_8, h_surf, domain%mesh_p)

    call set_vector_test_field(u, v, velocity_field, domain%mesh_u, &
                               domain%mesh_v, 0, v_components_type)

end subroutine

end module Williamson_test5_mod