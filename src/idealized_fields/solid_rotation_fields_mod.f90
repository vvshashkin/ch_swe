module solid_rotation_fields_mod

use vector_field_latlon_mod,    only : vector_field_latlon_t
use vector_field_cartesian_mod, only : vector_field_cartesian_t
use scalar_field_cartesian_mod, only : scalar_field_cartesian_t
use const_mod,                  only : pref=>P_reference, kappa, Cp, rgaz, &
                                       Day24h_sec, pi, Earth_grav, Earth_radii

implicit none

type, extends(vector_field_latlon_t) :: solid_rotation_vector_field_t

    real(kind=8) :: u0 = 1.0_8
    real(kind=8) :: alpha = 0.0_8
    real(kind=8) :: w_max  = 0.0_8
    real(kind=8) :: tau   = 12.0_8*Day24h_sec
    real(kind=8) :: Lz    = 10e3

    contains
        procedure, public :: get_latlon_vec
end type

type, extends(vector_field_cartesian_t) :: solid_rotation_vector_field_cart_t

    real(kind=8) :: u0 = 1.0_8
    real(kind=8) :: alpha = 0.0_8
    real(kind=8) :: w_max  = 0.0_8
    real(kind=8) :: tau   = 12.0_8*Day24h_sec
    real(kind=8) :: Lz    = 10e3

    contains
        procedure, public :: get_cartesian_vec
end type

type, extends(scalar_field_cartesian_t) :: solid_rotation_isotermic_theta_t
    real(kind=8) :: u0, alpha, T0, omega, sphere_rad, grav, p0 = 93000.0_8
    contains
        procedure, public :: get_scalar_cartesian => get_isotermic_theta
end type

type, extends(scalar_field_cartesian_t) :: solid_rotation_isotermic_PExner_t
    real(kind=8) :: u0, alpha, T0, omega, sphere_rad, grav, p0 = 93000.0_8
    contains
        procedure, public :: get_scalar_cartesian => get_isotermic_Pexner
end type

type, extends(scalar_field_cartesian_t) :: solid_rotation_isotermic_rho_t
    real(kind=8) :: u0=1.0_8, alpha=0.0_8,T0=300.0_8,omega=0.0_8,&
                    sphere_rad=Earth_radii, grav=Earth_grav, p0 = 93000.0_8
    contains
        procedure, public :: get_scalar_cartesian => get_isotermic_rho
end type

type, extends(scalar_field_cartesian_t) :: solid_rotation_Nbw_theta_t
    real(kind=8) :: u0, alpha, T0, omega, Nbw, sphere_rad, grav, peq = 1e5_8
    contains
    procedure :: get_scalar_cartesian => get_Nbw_theta
end type

type, extends(scalar_field_cartesian_t) :: solid_rotation_Nbw_PExner_t
    real(kind=8) :: u0, alpha, T0, omega, Nbw, sphere_rad, grav, peq = 1e5_8
    contains
    procedure :: get_scalar_cartesian => get_Nbw_PExner
end type

contains

pure function get_latlon_vec(this,phi,lambda,h,klev,time) result(uvw)

    class(solid_rotation_vector_field_t), intent(in) :: this
    real(kind=8),                         intent(in) :: lambda, phi, h, time
    integer(kind=4),                      intent(in) :: klev

    real(kind=8) :: uvw(3)

    uvw(1) = this%u0*(cos(this%alpha)*cos(phi)-sin(this%alpha)*cos(lambda)*sin(phi))
    uvw(2) = this%u0*sin(this%alpha)*sin(lambda)
    uvw(3) = this%w_max*sin(pi*h/this%Lz)*cos(2._8*pi*time/this%tau)
end function

pure function get_cartesian_vec(this,x,y,z,h,klev,time) result(v)

    class(solid_rotation_vector_field_cart_t), intent(in) :: this
    real(kind=8),                              intent(in) :: x,y,z,h,time
    integer(kind=4),                           intent(in) :: klev

    real(kind=8) :: v(4)
    real(kind=8) :: rot_axis(3)

    rot_axis(1:3) = rotation_axis(this%alpha)![sin(this%alpha), 0.0_8, cos(this%alpha)]

    v(1) = this%u0*(rot_axis(2)*z-rot_axis(3)*y)
    v(2) = this%u0*(-rot_axis(1)*z+rot_axis(3)*x)
    v(3) = this%u0*(rot_axis(1)*y-rot_axis(2)*x)
    v(4) = this%w_max*sin(pi*h/this%Lz)*cos(2._8*pi*time/this%tau)
end function

pure elemental real(kind=8) function get_isotermic_theta(this,x,y,z,h,klev,&
                                                                 time) result(f)
    class(solid_rotation_isotermic_theta_t), intent(in) :: this
    real(kind=8),                            intent(in) :: x,y,z,h
    integer(kind=4),                         intent(in) :: klev
    real(kind=8),                            intent(in) :: time

    real(kind=8) :: ps, phi, axis(3)

    phi = calc_axis_phi(x,y,z,this%alpha)

    ps = get_isotermic_ps(this%u0,this%omega,this%sphere_rad,this%T0,&
                                                        this%p0,this%grav,phi)
    f = this%T0 * (pref/ps)**kappa*exp(this%grav*h / (Cp*this%T0))

end function

pure elemental real(kind=8) function get_isotermic_PExner(this,x,y,z,h,klev,&
                                                                 time) result(f)
    class(solid_rotation_isotermic_PExner_t), intent(in) :: this
    real(kind=8),                             intent(in) :: x,y,z,h
    integer(kind=4),                          intent(in) :: klev
    real(kind=8),                             intent(in) :: time

    real(kind=8) :: ps, phi

    phi = calc_axis_phi(x,y,z,this%alpha)

    ps = get_isotermic_ps(this%u0,this%omega,this%sphere_rad,this%T0,&
                                                        this%p0,this%grav,phi)
    f = (ps/pref)**kappa*exp(-this%grav*h / (Cp*this%T0))

end function

pure elemental real(kind=8) function get_isotermic_rho(this,x,y,z,h,klev,&
                                                                 time) result(f)
    class(solid_rotation_isotermic_rho_t), intent(in) :: this
    real(kind=8),                          intent(in) :: x,y,z,h
    integer(kind=4),                       intent(in) :: klev
    real(kind=8),                          intent(in) :: time

    real(kind=8) :: ps, phi

    phi = calc_axis_phi(x,y,z,this%alpha)

    ps = get_isotermic_ps(this%u0,this%omega,this%sphere_rad,this%T0,&
                                                        this%p0,this%grav,phi)
    f = ps*exp(-this%grav*h / (Rgaz*this%T0)) / (Rgaz*this%T0)

end function

pure elemental real(kind=8) function get_Nbw_theta(this,x,y,z,h,klev,&
                                                                 time) result(f)
    class(solid_rotation_Nbw_theta_t),   intent(in) :: this
    real(kind=8),                        intent(in) :: x,y,z,h
    integer(kind=4),                     intent(in) :: klev
    real(kind=8),                        intent(in) :: time

    real(kind=8) :: ps, phi, T, P

    phi = calc_axis_phi(x,y,z,this%alpha)

    T = calc_T_Nbw_const(phi,h,this%u0,this%omega,this%sphere_rad,this%Nbw,this%grav,this%T0)
    P = calc_PExner_Nbw_const(phi,h,this%u0,this%omega,this%sphere_rad,this%Nbw,&
                                                       this%grav,this%T0,this%peq)
    f = T/P

end function

pure elemental real(kind=8) function get_Nbw_PExner(this,x,y,z,h,klev,&
                                                                 time) result(f)
    class(solid_rotation_Nbw_PExner_t),  intent(in) :: this
    real(kind=8),                        intent(in) :: x,y,z,h
    integer(kind=4),                     intent(in) :: klev
    real(kind=8),                        intent(in) :: time

    real(kind=8) :: ps, phi, T, P

    phi = calc_axis_phi(x,y,z,this%alpha)

    f = calc_PExner_Nbw_const(phi,h,this%u0,this%omega,this%sphere_rad,this%Nbw,&
                                                       this%grav,this%T0,this%peq)
end function

pure elemental real(kind=8) function calc_axis_phi(x,y,z,alpha) result(phi)
    real(kind=8), intent(in) :: x,y,z,alpha

    real(kind=8) :: axis(3)

    axis(1:3) = rotation_axis(alpha)
    phi = asin((x*axis(1)+y*axis(2)+z*axis(3)) / sqrt(x**2+y**2+z**2))

end function

pure function rotation_axis(alpha) result(axis)
    real(kind=8), intent(in) :: alpha
    real(kind=8) :: axis(3)
    axis = [sin(alpha), 0.0_8, cos(alpha)]
end function

real(kind=8) pure function get_isotermic_ps(u0,omega,a,T0,p0,grav,phi) result(ps)
    real(kind=8), intent(in) ::  u0,omega,a,T0,p0,grav,phi

    real(kind=8) :: Nb2

    Nb2 = grav**2 / (Cp*T0)
    ps = p0*exp(-0.5_8*Nb2*u0/(grav**2*kappa)*(u0+2.0_8*omega*a)*(sin(phi)**2-1.0_8))
end function

pure elemental real(kind=8) function calc_PExner_Nbw_const(phi,h,u0,omega,a,Nbw,&
                                                           grav,T0,peq) result(P)
    real(kind=8), intent(in) ::  u0,omega,a,Nbw,grav,phi,T0,peq,h

    real(kind=8) :: G, T

    G = grav**2 / (Nbw**2*Cp)
    T = calc_T_Nbw_const(phi,h,u0,omega,a,Nbw,grav,T0)
    P = (peq/pref)**kappa*(T/T0)* &
        exp(0.25_8*u0/(G*Cp)*(u0+2.0_8*omega*a)*(cos(2.0_8*phi)-1.0_8) - grav*h/(G*Cp))
end function

real(kind=8) pure function calc_T_Nbw_const(phi,h,u0,omega,a,Nbw,grav,T0) result(T)
    real(kind=8), intent(in) ::  phi, h, u0,omega,a,Nbw,grav,T0

    real(kind=8) :: G, exponent, Tsurf

    exponent = exp(Nbw**2*h/grav)
    G = grav**2 / (Nbw**2*Cp)
    Tsurf = G+(T0-G)*exp(-0.25_8*u0*Nbw**2/grav**2 * (u0+2.0_8*omega*a)*(cos(2.0_8*phi)-1.0_8))
    T = G*(1._8-exponent)+Tsurf*exponent

end function

end module
