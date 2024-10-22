module vertical_profiles_mod

    use scalar_field_cartesian_mod,        only : scalar_field_cartesian_t
    use const_mod,                         only : Earth_grav, Cp, Rgaz, P_reference, kappa

    implicit none
    
    real(kind=8),    parameter :: BW_sigma_t = 0.2_8   !tropopause level 200hPa
    real(kind=8),    parameter :: BW_Gamma = 0.005_8 !Temperature lapse rate K/m
    real(kind=8),    parameter :: BW_T0 = 288.0_8, BW_delta_T = 4.8e5
    real(kind=8),    parameter :: BW_tolerance = 1e-10_8
    integer(kind=4), parameter :: BW_maxit = 15

    !constant Brunt-Vaisala frequency profiles:

    type, extends(scalar_field_cartesian_t) :: const_Nbw_profile_theta_t
        real(kind=8) :: Nbw, t0, grav = Earth_grav, Cp = Cp
        contains
        procedure, public :: get_scalar_cartesian => get_const_Nbw_theta
    end type
    
    type, extends(const_Nbw_profile_theta_t) :: const_Nbw_profile_dtheta_dz_t
        contains
        procedure, public :: get_scalar_cartesian => get_const_Nbw_dtheta_dz
    end type

    type, extends(const_Nbw_profile_theta_t) :: const_Nbw_profile_P_t
        contains
        procedure, public :: get_scalar_cartesian => get_const_Nbw_P
    end type

    type, extends(const_Nbw_profile_theta_t) :: const_Nbw_profile_dP_dz_t
        contains
        procedure, public :: get_scalar_cartesian => get_const_Nbw_dP_dz
    end type

    !Isotermic profiles

    type, extends(scalar_field_cartesian_t) :: isotermal_profile_theta_t
        real(kind=8) :: t0, p_surf, Cp = Cp, Rgaz=Rgaz, grav = Earth_grav
        contains
        procedure, public :: get_scalar_cartesian => get_isotermal_theta
    end type
    
    type, extends(isotermal_profile_theta_t) :: isotermal_profile_dtheta_dz_t
        contains
        procedure, public :: get_scalar_cartesian => get_isotermal_dtheta_dz
    end type

    type, extends(isotermal_profile_theta_t) :: isotermal_profile_P_t
        contains
        procedure, public :: get_scalar_cartesian => get_isotermal_P
    end type

    type, extends(isotermal_profile_theta_t) :: isotermal_profile_dP_dz_t
        contains
        procedure, public :: get_scalar_cartesian => get_isotermal_dP_dz
    end type

    !Jablonowski Williamson baroclinic wave background profiles

    type, extends(scalar_field_cartesian_t) :: BW_profile_theta_t
        contains
        procedure, public :: get_scalar_cartesian => get_BW_theta
    end type

    type, extends(BW_profile_theta_t) :: BW_profile_dtheta_dz_t
       contains
      procedure, public :: get_scalar_cartesian => get_BW_dtheta_dz
    end type

    type, extends(BW_profile_theta_t) :: BW_profile_P_t
        contains
        procedure, public :: get_scalar_cartesian => get_BW_P
    end type

    type, extends(BW_profile_theta_t) :: BW_profile_dP_dz_t
        contains
        procedure, public :: get_scalar_cartesian => get_BW_dP_dz
    end type

    contains

    !constant Nbw:

    pure elemental real(kind=8) function get_const_Nbw_theta(this,x,y,z,h,klev,time) &
                                                                        result(f)
        class(const_Nbw_profile_theta_t), intent(in) :: this
        real(kind=8),                     intent(in) :: x,y,z,h
        integer(kind=4),                  intent(in) :: klev
        real(kind=8),                     intent(in) :: time

        f = this%t0*exp(this%Nbw**2*h / this%grav)

    end function
    
    pure elemental real(kind=8) function get_const_Nbw_dtheta_dz(this,x,y,z,h,klev,time) &
                                                                        result(f)
        class(const_Nbw_profile_dtheta_dz_t), intent(in) :: this
        real(kind=8),                         intent(in) :: x,y,z,h
        integer(kind=4),                      intent(in) :: klev
        real(kind=8),                         intent(in) :: time

        f = this%Nbw**2/this%grav*this%t0*exp(this%Nbw**2*h / this%grav)

    end function

    pure elemental real(kind=8) function get_const_Nbw_P(this,x,y,z,h,klev,time) &
                                                                        result(f)
        class(const_Nbw_profile_P_t), intent(in) :: this
        real(kind=8),                 intent(in) :: x,y,z,h
        integer(kind=4),              intent(in) :: klev
        real(kind=8),                 intent(in) :: time

        f = 1.0_8+this%grav**2/(this%Nbw**2*this%Cp*this%t0)*(exp(-this%Nbw**2/this%grav*h)-1.0_8)

    end function

    pure elemental real(kind=8) function get_const_Nbw_dP_dz(this,x,y,z,h,klev,time) &
                                                                        result(f)
        class(const_Nbw_profile_dP_dz_t), intent(in) :: this
        real(kind=8),                     intent(in) :: x,y,z,h
        integer(kind=4),                  intent(in) :: klev
        real(kind=8),                     intent(in) :: time

        f = -this%grav/(this%Cp*this%t0)*exp(-this%Nbw**2/this%grav*h)

    end function

    !Isotermal

    pure elemental real(kind=8) function get_isotermal_theta(this,x,y,z,h,klev,time) &
                                                                        result(f)
        class(isotermal_profile_theta_t), intent(in) :: this
        real(kind=8),                     intent(in) :: x,y,z,h
        integer(kind=4),                  intent(in) :: klev
        real(kind=8),                     intent(in) :: time

        f = this%t0*(P_reference/this%p_surf*exp(h*this%grav / (this%t0*this%rgaz)))**(this%rgaz/this%Cp)

    end function

    pure elemental real(kind=8) function get_isotermal_dtheta_dz(this,x,y,z,h,klev,time) &
                                                                        result(f)
        class(isotermal_profile_dtheta_dz_t), intent(in) :: this
        real(kind=8),                         intent(in) :: x,y,z,h
        integer(kind=4),                      intent(in) :: klev
        real(kind=8),                         intent(in) :: time

        real(kind=8) :: kappa

        kappa = this%Rgaz / this%Cp

        f = this%t0*(P_reference/this%p_surf)**(kappa)* &
                exp(kappa*h*this%grav / (this%t0*this%rgaz))*&
                    kappa*this%grav / (this%t0*this%rgaz)

    end function

    pure elemental real(kind=8) function get_isotermal_P(this,x,y,z,h,klev,time) &
                                                                        result(f)
        class(isotermal_profile_P_t), intent(in) :: this
        real(kind=8),                 intent(in) :: x,y,z,h
        integer(kind=4),              intent(in) :: klev
        real(kind=8),                 intent(in) :: time

        real(kind=8) :: kappa

        kappa = this%Rgaz / this%Cp

        f = (this%p_surf/P_reference*exp(-h*this%grav / (this%rgaz*this%t0)))**kappa

    end function

    pure elemental real(kind=8) function get_isotermal_dP_dz(this,x,y,z,h,klev,time) &
                                                                        result(f)
        class(isotermal_profile_dP_dz_t), intent(in) :: this
        real(kind=8),                     intent(in) :: x,y,z,h
        integer(kind=4),                  intent(in) :: klev
        real(kind=8),                     intent(in) :: time

        real(kind=8) :: kappa

        kappa = this%Rgaz / this%Cp

        f = -(this%p_surf/P_reference)**kappa*exp(-kappa*h*this%grav / (this%rgaz*this%t0))*&
                                               kappa*this%grav / (this%rgaz*this%t0)

    end function

    !Baroclinic wave background

    pure elemental real(kind=8) function get_BW_theta(this,x,y,z,h,klev,time) &
                                                                        result(f)
        class(BW_profile_theta_t), intent(in) :: this
        real(kind=8),              intent(in) :: x,y,z,h
        integer(kind=4),           intent(in) :: klev
        real(kind=8),              intent(in) :: time

        real(kind=8) :: sigma

        sigma = solve_BW_h_for_sigma(h)

        f = calc_BW_temp(sigma) / sigma**kappa

    end function

    pure elemental real(kind=8) function get_BW_dtheta_dz(this,x,y,z,h,klev,time) &
                                                                        result(f)
        class(BW_profile_dtheta_dz_t), intent(in) :: this
        real(kind=8),                  intent(in) :: x,y,z,h
        integer(kind=4),               intent(in) :: klev
        real(kind=8),                  intent(in) :: time

        real(kind=8) :: sigma, temp, temp1, dsigma_dz

        sigma = solve_BW_h_for_sigma(h)

        temp = calc_BW_temp(sigma)
        temp1 = calc_BW_dtemp_dsigma(sigma)

        dsigma_dz = -Earth_grav*sigma / (Rgaz*temp)

        f = dsigma_dz*(temp1 - kappa*temp/sigma) / sigma**kappa

    end function

    pure elemental real(kind=8) function get_BW_P(this,x,y,z,h,klev,time) &
                                                                    result(f)
        class(BW_profile_P_t), intent(in) :: this
        real(kind=8),          intent(in) :: x,y,z,h
        integer(kind=4),       intent(in) :: klev
        real(kind=8),          intent(in) :: time

        f = solve_BW_h_for_sigma(h)**kappa

    end function

    pure elemental real(kind=8) function get_BW_dP_dz(this,x,y,z,h,klev,time) &
                                                                      result(f)
        class(BW_profile_dP_dz_t), intent(in) :: this
        real(kind=8),              intent(in) :: x,y,z,h
        integer(kind=4),           intent(in) :: klev
        real(kind=8),              intent(in) :: time

        real(kind=8) :: sigma, temp, dsigma_dz

        sigma = solve_BW_h_for_sigma(h)

        temp = calc_BW_temp(sigma)

        dsigma_dz = -Earth_grav*sigma / (Rgaz*temp)

        f = kappa*sigma**(kappa-1.0_8)*dsigma_dz

    end function

    pure elemental function solve_BW_h_for_sigma(h) result(sigma)
        real(kind=8), intent(in) :: h
        real(kind=8)             :: sigma

        real(kind=8)    :: temp, dh
        integer(kind=4) :: it

        sigma = min(1.0_8,max(max(1.0_8-BW_Gamma*h/BW_T0,0.0_8)**(Earth_grav/(Rgaz*BW_Gamma)),1e-5_8))

        do it = 1, BW_maxit
            temp = calc_BW_Temp(sigma)
            dh = h-calc_BW_H(sigma)
            sigma = min(1.0_8,max(1e-6,sigma*(1.0_8-dh*Earth_grav/(Rgaz*temp))))
            if(abs(dh) < BW_tolerance) exit
        end do
    end function

    pure elemental function calc_BW_temp(sigma) result(temp)

        real(kind=8), intent(in) :: sigma
        real(kind=8) :: temp

        temp = BW_T0*sigma**(Rgaz*BW_Gamma/Earth_grav)+ BW_delta_T*max(0.0_8,BW_sigma_t-sigma)**5

    end function

    pure elemental function calc_BW_dtemp_dsigma(sigma) result(temp)

        real(kind=8), intent(in) :: sigma
        real(kind=8) :: temp

        temp = BW_T0*(Rgaz*BW_Gamma/Earth_grav)*sigma**(Rgaz*BW_Gamma/Earth_grav-1.0_8) - 5.0_8*BW_delta_T*max(0.0_8,BW_sigma_t-sigma)**4

    end function

    pure elemental function calc_BW_H(sigma) result(H) !geopotential at level p0*sigma

        real(kind=8), intent(in)    :: sigma

        real(kind=8) :: H

        H = BW_T0/BW_Gamma*(1.0_8-sigma**(Rgaz*BW_Gamma/Earth_grav))

        if(sigma<BW_sigma_t) &
            H = H-Rgaz*BW_delta_T * ((log(sigma/BW_sigma_t)+137.0_8/60.0_8)*BW_sigma_t**5            - &
                                 5.0_8*BW_sigma_t**4*sigma+5.0_8*BW_sigma_t**3*sigma**2         - &
                                 10.0_8/3.0_8*BW_sigma_t**2*sigma**3+1.25_8*BW_sigma_t*sigma**4 - &
                                 0.2_8*sigma**5)/Earth_grav

    end function

end module