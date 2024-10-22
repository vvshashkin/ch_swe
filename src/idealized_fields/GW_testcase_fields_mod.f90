module GW_testcase_fields_mod

use scalar_field_cartesian_mod, only : scalar_field_cartesian_t
use const_mod,                  only : pi, Earth_grav, Rgaz, Cv, pref=>P_reference

implicit none

real(kind=8), private, parameter :: GW_omega = 0.0_8, GW_T0 = 300.0_8, GW_peq = 1e5_8

type, extends(scalar_field_cartesian_t) :: GW_theta_t
    real(kind=8) :: u0, d, sphere_rad, Nbw, Lz, amp
    contains
        procedure :: get_scalar_cartesian => get_theta
end type

type, extends(scalar_field_cartesian_t) :: GW_theta_linear_t
    real(kind=8) :: u0, d, sphere_rad, Nbw, Lz, amp
    contains
        procedure :: get_scalar_cartesian => get_theta_linear
end type

type, extends(scalar_field_cartesian_t) :: GW_PExner_t
    real(kind=8) :: u0, sphere_rad, Nbw
    contains
        procedure :: get_scalar_cartesian => get_PExner
end type

type, extends(scalar_field_cartesian_t) :: GW_rho_t
    real(kind=8) :: u0, d, sphere_rad, Nbw, Lz, amp
    contains
        procedure :: get_scalar_cartesian => get_rho
end type

contains

pure elemental real(kind=8) function get_theta(this,x,y,z,h,klev,time) result(f)

    use solid_rotation_fields_mod, only : calc_T_Nbw_const, calc_PExner_Nbw_const

    class(GW_theta_t),   intent(in) :: this
    real(kind=8),        intent(in) :: x,y,z,h
    integer(kind=4),     intent(in) :: klev
    real(kind=8),        intent(in) :: time

    real(kind=8) :: phi, T, P, d

    phi = asin(z/sqrt(x**2+y**2+z**2))

    T = calc_T_Nbw_const(phi,h,this%u0,GW_omega,this%sphere_rad,this%Nbw,Earth_grav,GW_T0)
    P = calc_PExner_Nbw_const(phi,h,this%u0,GW_omega,this%sphere_rad,this%Nbw,&
                                                       Earth_grav,GW_T0,GW_peq)

    d = acos(x/sqrt(x**2+y**2+z**2))*this%sphere_rad

    f = T/P + this%amp*this%d**2 / (this%d**2+d**2)*sin(pi*h/this%Lz)

end function

pure elemental real(kind=8) function get_theta_linear(this,x,y,z,h,klev,time) result(f)

    class(GW_theta_linear_t),   intent(in) :: this
    real(kind=8),               intent(in) :: x,y,z,h
    integer(kind=4),            intent(in) :: klev
    real(kind=8),               intent(in) :: time

    real(kind=8) :: d

    d = acos(x/sqrt(x**2+y**2+z**2))*this%sphere_rad

    f = this%amp*this%d**2 / (this%d**2+d**2)*sin(pi*h/this%Lz)

end function

pure elemental real(kind=8) function get_PExner(this,x,y,z,h,klev,time) result(f)

    use solid_rotation_fields_mod, only : calc_PExner_Nbw_const

    class(GW_PExner_t),  intent(in) :: this
    real(kind=8),        intent(in) :: x,y,z,h
    integer(kind=4),     intent(in) :: klev
    real(kind=8),        intent(in) :: time

    real(kind=8) :: phi, T, P, d

    phi = asin(z/sqrt(x**2+y**2+z**2))

    f = calc_PExner_Nbw_const(phi,h,this%u0,GW_omega,this%sphere_rad,this%Nbw,&
                                                       Earth_grav,GW_T0,GW_peq)

end function

pure elemental real(kind=8) function get_rho(this,x,y,z,h,klev,time) result(f)

    use solid_rotation_fields_mod, only : calc_T_Nbw_const, calc_PExner_Nbw_const

    class(GW_rho_t),     intent(in) :: this
    real(kind=8),        intent(in) :: x,y,z,h
    integer(kind=4),     intent(in) :: klev
    real(kind=8),        intent(in) :: time

    real(kind=8) :: phi, T, P, d, theta

    phi = asin(z/sqrt(x**2+y**2+z**2))

    T = calc_T_Nbw_const(phi,h,this%u0,GW_omega,this%sphere_rad,this%Nbw,Earth_grav,GW_T0)
    P = calc_PExner_Nbw_const(phi,h,this%u0,GW_omega,this%sphere_rad,this%Nbw,&
                                                       Earth_grav,GW_T0,GW_peq)

    d = acos(x/sqrt(x**2+y**2+z**2))*this%sphere_rad

    theta = T/P + this%amp*this%d**2 / (this%d**2+d**2)*sin(pi*h/this%Lz)

    f = pref*P**(Cv/Rgaz) / (Rgaz*theta)

end function

end module
