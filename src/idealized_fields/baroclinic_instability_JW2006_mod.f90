module baroclinic_instability_JW2006_mod

use vector_field_latlon_mod,    only : vector_field_latlon_t
use scalar_field_latlon_mod,    only : scalar_field_latlon_t
use const_mod,                  only : pi, grav=>Earth_grav, Rd=>Rgaz, &
                                       omega => Earth_omega, &
                                       Ra => Earth_radii, kappa

implicit none

real(kind=8),    parameter :: p0 = 1e5_8
real(kind=8),    parameter :: sigma_0 = 0.252_8 !jet core level pressure 252hPa
real(kind=8),    parameter :: sigma_t = 0.2_8   !tropopause level 200hPa
real(kind=8),    parameter :: u0 = 35.0_8
real(kind=8),    parameter :: Gamma = 0.005_8 !Temperature lapse rate K/m
real(kind=8),    parameter :: T0 = 288.0_8, delta_T = 4.8e5
!perturbation position:
real(kind=8),    parameter :: lambda_c = pi/9.0_8, phi_c = 2.0_8*pi/9.0_8
!h-> sigma solver const:
real(kind=8),    parameter :: tolerance = 1e-10
integer(kind=4), parameter :: maxit = 15

type, extends(vector_field_latlon_t) :: baroclinic_instability_wind_t

    real(kind=8) :: u_pert = 1.0_8 ! velocity perturbation amplitude

    contains
        procedure, public :: get_latlon_vec
end type

!for postprocessing test:
type, extends(vector_field_latlon_t) :: baroclinic_instability_isobaric_wind_t

    real(kind=8) :: p =1e5_8, u_pert = 1.0_8

    contains
        procedure, public :: get_latlon_vec => get_isobaric_wind
end type

type, extends(scalar_field_latlon_t) :: baroclinic_instability_PExner_t
    contains
        procedure, public :: get_scalar_latlon => get_Pexner
end type

type, extends(scalar_field_latlon_t) :: baroclinic_instability_theta_t
    contains
    procedure :: get_scalar_latlon => get_theta
end type

type, extends(scalar_field_latlon_t) :: baroclinic_instability_rho_t
    contains
    procedure :: get_scalar_latlon => get_rho
end type

contains

pure function get_latlon_vec(this,phi,lambda,h,klev,time) result(uvw)

    class(baroclinic_instability_wind_t), intent(in) :: this
    real(kind=8),                         intent(in) :: lambda, phi, h, time
    integer(kind=4),                      intent(in) :: klev

    real(kind=8) :: uvw(3)

    real(kind=8) :: r, sigma_v

    sigma_v = (solve_h_for_sigma(h,phi)-sigma_0)*0.5_8*pi

    r = acos(cos(phi_c)*cos(phi)*cos(lambda-lambda_c)+sin(phi)*sin(phi_c))

    uvw(1) = u0*cos(sigma_v)**1.5_8*sin(2.0_8*phi)**2 + &
                                    this%u_pert*exp(-100.0_8*r**2)
    uvw(2) = 0.0_8
    uvw(3) = 0.0_8

end function

pure function get_isobaric_wind(this,phi,lambda,h,klev,time) result(uvw)

    class(baroclinic_instability_isobaric_wind_t), intent(in) :: this
    real(kind=8),                                  intent(in) :: lambda, phi, h, time
    integer(kind=4),                               intent(in) :: klev

    real(kind=8) :: uvw(3)

    real(kind=8) :: r, sigma_v

    sigma_v = (this%p/p0-sigma_0)*0.5_8*pi

    r = acos(cos(phi_c)*cos(phi)*cos(lambda-lambda_c)+sin(phi)*sin(phi_c))

    uvw(1) = u0*cos(sigma_v)**1.5_8*sin(2.0_8*phi)**2 + &
                                    this%u_pert*exp(-100.0_8*r**2)
    uvw(2) = 0.0_8
    uvw(3) = 0.0_8

end function

pure elemental real(kind=8) function get_PExner(this,phi,lambda,h,klev,time) result(f)
    class(baroclinic_instability_PExner_t), intent(in) :: this
    real(kind=8),                           intent(in) :: phi, lambda, h
    integer(kind=4),                        intent(in) :: klev
    real(kind=8),                           intent(in) :: time

    f = solve_h_for_sigma(h,phi)**kappa

end function

pure elemental real(kind=8) function get_theta(this,phi,lambda,h,klev,time) result(f)
    class(baroclinic_instability_theta_t),  intent(in) :: this
    real(kind=8),                           intent(in) :: phi, lambda, h
    integer(kind=4),                        intent(in) :: klev
    real(kind=8),                           intent(in) :: time

    real(kind=8) :: sigma

    sigma = solve_h_for_sigma(h, phi)

    f = calc_temp(sigma,phi)/sigma**kappa

end function

pure elemental real(kind=8) function get_rho(this,phi,lambda,h,klev,time) result(f)
    class(baroclinic_instability_rho_t),    intent(in) :: this
    real(kind=8),                           intent(in) :: phi, lambda, h
    integer(kind=4),                        intent(in) :: klev
    real(kind=8),                           intent(in) :: time

    real(kind=8) :: sigma

    sigma = solve_h_for_sigma(h, phi)

    f = p0*sigma / (Rd*calc_temp(sigma,phi))

end function

pure function solve_h_for_sigma(h,phi) result(sigma)
    real(kind=8), intent(in) :: h, phi
    real(kind=8)             :: sigma

    real(kind=8)    :: temp, dh
    integer(kind=4) :: it

    sigma = min(1.0_8,max(max(1.0_8-Gamma*h/T0,0.0_8)**(grav/(Rd*Gamma)),1e-5_8))

    do it = 1, maxit
        temp = calc_Temp(sigma,phi)
        dh = h-calc_H(sigma,phi)
        sigma = min(1.0_8,max(1e-6,sigma*(1.0_8-dh*grav/(Rd*temp))))
        if(abs(dh) < tolerance) exit
    end do
end function

pure function calc_temp(sigma, phi) result(temp)

    real(kind=8), intent(in) :: sigma, phi
    real(kind=8) :: temp

    real(kind=8) :: A, B, sigma_v

    sigma_v = 0.5*pi*(sigma-sigma_0)
    A = (-2.0_8*sin(phi)**6*(cos(phi)**2+1.0_8/3.0_8)+10.0_8/63.0_8)
    B = (8.0_8/5.0_8*cos(phi)**3*(sin(phi)**2+2.0_8/3.0_8)-0.25_8*pi)

    temp = T0*sigma**(Rd*Gamma/grav)+ &
           delta_T*max(0.0_8,sigma_t-sigma)**5+ &
           0.75_8*sigma*pi*u0/Rd*sin(sigma_v)*sqrt(cos(sigma_v))* &
          (A*2.0_8*u0*cos(sigma_v)**1.5_8+B*Ra*omega)

end function

pure function calc_H(sigma,phi) result(H) !geopotential at level p0*sigma

    real(kind=8), intent(in)    :: sigma
    real(kind=8), intent(in)    :: phi

    real(kind=8) :: H

    real(kind=8)    :: sigma_v, A, B, u

    sigma_v = 0.5*pi*(sigma-sigma_0)
    u = u0*cos(sigma_v)**1.5_8
    A = (-2.0_8*sin(phi)**6*(cos(phi)**2+1.0_8/3.0_8)+10.0_8/63.0_8)
    B = (8.0_8/5.0_8*cos(phi)**3*(sin(phi)**2+2.0_8/3.0_8)-0.25_8*pi)

    H = T0/Gamma*(1.0_8-sigma**(Rd*Gamma/grav))+u*(A*u+B*Ra*omega)/grav

    if(sigma<sigma_t) &
        H = H-Rd*delta_T * ((log(sigma/sigma_t)+137.0_8/60.0_8)*sigma_t**5            - &
                             5.0_8*sigma_t**4*sigma+5.0_8*sigma_t**3*sigma**2         - &
                             10.0_8/3.0_8*sigma_t**2*sigma**3+1.25_8*sigma_t*sigma**4 - &
                             0.2_8*sigma**5)/grav

end function

end module
