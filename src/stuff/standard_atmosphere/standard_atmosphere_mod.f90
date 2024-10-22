module standard_atmosphere_mod

use const_mod, only : Rgaz, grav=>Earth_grav, Cp

implicit none

private
public :: get_standard_atm_temp_z, get_standard_atm_pres_z, get_standard_atm_height_p, &
          get_standard_atm_temp_p, get_standard_atm_theta_z, get_standard_atm_theta_p

real(kind=8), parameter :: z0(7) = [-610.0_8,  11.0e3_8, 20.0e3_8, 32.0e3_8, &
                                     47.0e3_8, 51.0e3_8, 71.0e3_8]
real(kind=8), parameter :: dTdz(7) = [-6.5e-3_8, 0.0_8, 1.0e-3_8, 2.8e-3_8, &
                                       0.0_8, -2.8e-3_8, -2.0e-3_8]
real(kind=8), parameter :: T1 = 292.150_8, T2 = T1+dTdz(1)*(z0(2)-z0(1)), &
                           T3 = T2+dTdz(2)*(z0(3)-z0(2)), T4 = T3+dTdz(3)*(z0(4)-z0(3)), &
                           T5 = T4+dTdz(4)*(z0(5)-z0(4)), T6 = T5+dTdz(5)*(z0(6)-z0(5)), &
                           T7 = T6+dTdz(6)*(z0(7)-z0(6))
real(kind=8), parameter :: T0(7) = [T1, T2, T3, T4, T5, T6, T7]

real(kind=8), parameter :: pref = 1e5_8 !for potential temp calculation
real(kind=8), parameter :: p1 = 108900.0_8
real(kind=8), parameter :: p2 = p1*(1.0_8+dTdz(1)/T0(1)*(z0(2)-z0(1)))**(-grav/(Rgaz*dTdz(1)))
real(kind=8), parameter :: p3 = p2*exp(-grav*(z0(3)-z0(2))/(Rgaz*T0(2)))
real(kind=8), parameter :: p4 = p3*(1.0_8+dTdz(3)/T0(3)*(z0(4)-z0(3)))**(-grav/(Rgaz*dTdz(3)))
real(kind=8), parameter :: p5 = p4*(1.0_8+dTdz(4)/T0(4)*(z0(5)-z0(4)))**(-grav/(Rgaz*dTdz(4)))
real(kind=8), parameter :: p6 = p5*exp(-grav*(z0(6)-z0(5))/(Rgaz*T0(5)))
real(kind=8), parameter :: p7 = p6*(1.0_8+dTdz(6)/T0(6)*(z0(7)-z0(6)))**(-grav/(Rgaz*dTdz(6)))

real(kind=8), parameter :: p0(7) = [p1,p2,p3,p4,p5,p6,p7]

contains

pure elemental real(kind=8) function get_standard_atm_temp_z(z) result(temp)
    real(kind=8), intent(in) :: z

    integer(kind=4) :: k, ikk

    ikk = size(z0,1)
    do k=1, size(z0,1)-1
        if(z < z0(k+1)) then
            ikk = k
            exit
        end if
    end do

    temp = T0(ikk)+dTdz(ikk)*(z-z0(ikk))

end function

pure elemental real(kind=8) function get_standard_atm_pres_z(z) result(pres)
    real(kind=8), intent(in) :: z

    integer(kind=4) :: k, ikk

    ikk = size(z0,1)
    do k=1, size(z0,1)-1
        if(z < z0(k+1)) then
            ikk = k
            exit
        end if
    end do

    if(dTdz(ikk) == 0.0_8) then
        pres = p0(ikk)*exp(-grav*(z-z0(ikk))/(rgaz*T0(ikk)))
    else
        pres = p0(ikk)*(1.0_8+dTdz(ikk)/T0(ikk)*(z-z0(ikk)))**(-grav/(rgaz*dTdz(ikk)))
    end if

end function

pure elemental real(kind=8) function get_standard_atm_height_p(p) result(z)
    real(kind=8), intent(in) :: p

    integer(kind=4) :: k, ikk

    ikk = size(z0,1)
    do k=1, size(z0,1)-1
        if(p >= p0(k+1)) then
            ikk = k
            exit
        end if
    end do

    if(dTdz(ikk) == 0.0_8) then
        z = z0(ikk)+rgaz*T0(ikk)/grav*log(p0(ikk)/p)
    else
        z = z0(ikk)+T0(ikk)/dTdz(ikk)*((p0(ikk)/p)**(rgaz*dTdz(ikk)/grav)-1._8)
    end if

end function

pure elemental real(kind=8) function get_standard_atm_temp_p(p) result(temp)
    real(kind=8), intent(in) :: p

    temp = get_standard_atm_temp_z(get_standard_atm_height_p(p))

end function

pure elemental real(kind=8) function get_standard_atm_theta_z(z) result(theta)
    real(kind=8), intent(in) :: z

    theta = get_standard_atm_temp_z(z)*(pref/get_standard_atm_pres_z(z))**(rgaz/Cp)

end function

pure elemental real(kind=8) function get_standard_atm_theta_p(p) result(theta)
    real(kind=8), intent(in) :: p

    theta = get_standard_atm_temp_p(p)*(pref/p)**(rgaz/Cp)

end function

end module
