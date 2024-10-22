module barotropic_instability_mod

use test_fields_mod,    only : barotropic_instability_height_generator_t,    &
                               barotropic_instability_wind_generator_t,      &
                               set_vector_test_field, set_scalar_test_field
use domain_mod,         only : domain_t
use generic_config_mod, only : generic_config_t
use stvec_flexible_mod, only : stvec_flexible_t, stvec_t
use grid_field_mod,     only : grid_field_t

use const_mod,          only : pi, Earth_omega, Earth_radii, &
                               Earth_grav, Earth_sidereal_T

implicit none

real(kind=8), parameter, private :: H0_default     = 11000._8
real(kind=8), parameter, private :: u0_default     = 80.0_8
real(kind=8), parameter, private :: h_pert_default = 120.0_8

contains

subroutine setup_barotropic_instability_test(state, config, &
                                             v_components_type, domain)

    class(generic_config_t), intent(inout) :: config
    character(len=*),        intent(in)    :: v_components_type
    type(domain_t),          intent(in)    :: domain

    class(stvec_t),          intent(inout) :: state

    !locals
    integer(kind=4) :: N_quad_points
    real(kind=8)    :: omega, a, rotation_axis(3), u0, H0, h_pert
    type(barotropic_instability_height_generator_t) :: height_field
    type(barotropic_instability_wind_generator_t)   :: velocity_field
    type(grid_field_t), pointer :: h, u, v

    omega = domain%metric%omega
    rotation_axis = domain%metric%rotation_axis
    a = domain%metric%scale

    call config%get(u0,    "u0",    default=u0_default)
    call config%get(H0,    "H0",    default=H0_default)
    call config%get(h_pert,"h_pert",default=h_pert_default)
    call config%get(N_quad_points, "N_quad_points",default = 10*domain%partition%Nh)

    if(domain%parcomm%myid == 0) then

        if(u0 /= u0_default) &
            print *, "WARNING: Barotropic instability test using non-standard u0 = ", u0

        if(H0 /= H0_default) &
            print *, "WARNING: Barotropic instability test using non-standard H0 = ", H0

        if(h_pert /= h_pert_default) &
            print *, "WARNING: Barotropic instability test using non-standard h_pert = ", h_pert

        if(omega /= Earth_omega) &
            print *, "WARNING: Barotropic instability test using non-standard omega = ", omega

        if(a /= Earth_radii) &
            print *, "WARNING: Barotropic instability test using non-standard planet radius = ", a

        if(any(rotation_axis /= [0.0_8,0.0_8,1.0_8])) &
            print *, "WARNING: Barotropic instability test using non-standard Earth rotation axis = ", rotation_axis

    end if

    height_field = create_barotropic_instability_height_field_generator(H0 = H0, &
                              omega = omega, grav = Earth_grav, &
                              h_pert = h_pert, u0 = u0,   &
                              a = a, Nq = N_quad_points)

    velocity_field = barotropic_instability_wind_generator_t(u0=u0)

    select type(state)
    class is (stvec_flexible_t)

        call state%get_field(h,"h")
        call state%get_field(u,"u")
        call state%get_field(v,"v")

    end select

    call set_scalar_test_field(h, height_field, domain%mesh_p, 0)
    call set_vector_test_field(u, v, velocity_field, &
                               domain%mesh_u, domain%mesh_v, 0, v_components_type)

end subroutine setup_barotropic_instability_test

function barotropic_instability_u(u0, phi) result(u)
    use const_mod, only : pi

    real(kind=8), intent(in) :: u0, phi
    real(kind=8)             :: u

    real(kind=8), parameter :: phi0 = pi/7d0, phi1 = .5d0*pi-phi0

    u = u0/exp(-4._8/(phi1-phi0)**2) * exp(1._8/min((phi-phi0)*(phi-phi1),-1e-2))

end function barotropic_instability_u

function create_barotropic_instability_height_field_generator(H0, &
                              omega, grav, h_pert, u0, a, Nq) result (height_gen)

    use const_mod, only : pi

    real(kind=8), intent(in) :: H0
    real(kind=8), intent(in) :: omega
    real(kind=8), intent(in) :: grav
    real(kind=8), intent(in) :: h_pert
    real(kind=8), intent(in) :: u0
    real(kind=8), intent(in) :: a
    integer(kind=4), intent(in) :: Nq

    type(barotropic_instability_height_generator_t) :: height_gen

    real(kind=8), parameter :: phi0 = pi/7d0, phi1 = .5d0*pi-phi0
    real(kind=8)            :: dphi
    real(kind=8)            :: phi, u, dh
    real(kind=8)            :: urhs
    integer(kind=4)         :: j, iq

    real(kind=8), parameter :: xq(3) = [-sqrt(0.6_8), 0.0_8, sqrt(0.6_8)]
    real(kind=8), parameter :: wq(3) = [5.0_8/9.0_8, 8.0_8/9.0_8, 5.0_8/9.0_8]

    height_gen%H0 = H0
    height_gen%Nq = Nq
    allocate(height_gen%H_zonal(-1:Nq+1))
    dphi = (phi1 - phi0) / real(Nq,8)
    height_gen%dphi = dphi
    height_gen%h_pert = h_pert

    height_gen%H_zonal(-1) = H0!0.0_8
    height_gen%H_zonal( 0) = H0!0.0_8
    do j=1, Nq
        dh = 0.0_8
        do iq=1, size(wq)
            phi   = phi0+(j-0.5_8)*dphi+dphi*xq(iq)*0.5_8
            u = barotropic_instability_u(u0, phi)
            urhs  =-(u**2*tan(phi)+2._8*omega*a*sin(phi)*u) / grav
            dh = dh+wq(iq)*urhs
        end do
        height_gen%H_zonal(j) = height_gen%H_zonal(j-1)+0.5_8*dh*dphi
    end do
    height_gen%H_zonal(Nq+1) = height_gen%H_zonal(Nq)

    height_gen%H_north = height_gen%H_zonal(Nq)

end function

end module barotropic_instability_mod