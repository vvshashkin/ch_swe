module baroclinic_instability_test_therm_mod

use baroclinic_instability_test_parameters_mod, only : u0, sigma_0, sigma_t, &
                                                       Gamma, T0, delta_T
use const_mod, only : grav=>Earth_grav, omega=>Earth_omega, Ra=>Earth_radii, &
                      Rd=>rgaz, pi, kappa

implicit none

contains

subroutine calc_Temp(T,sigma,phi)
    real(kind=8), intent(inout) :: T(:,:,:)
    real(kind=8), intent(in)    :: sigma(:,:,:)
    real(kind=8), intent(in)    :: phi(:,:)

    integer(kind=4) :: i, j, k
    real(kind=8)    :: sigma_v, A, B

    do k = 1, size(T,3)
        do j = 1, size(T,2)
            do i = 1, size(T,1)
                sigma_v = 0.5*pi*(sigma(i,j,k)-sigma_0)
                A = (-2.0_8*sin(phi(i,j))**6*(cos(phi(i,j))**2+1.0_8/3.0_8)+10.0_8/63.0_8)
                B = (8.0_8/5.0_8*cos(phi(i,j))**3*(sin(phi(i,j))**2+2.0_8/3.0_8)-0.25_8*pi)

                T(i,j,k) = T0*sigma(i,j,k)**(Rd*Gamma/grav)+ &
                           delta_T*max(0.0_8,sigma_t-sigma(i,j,k))**5+ &
                           0.75_8*sigma(i,j,k)*pi*u0/Rd*sin(sigma_v)*sqrt(cos(sigma_v))* &
                           (A*2.0_8*u0*cos(sigma_v)**1.5_8+B*Ra*omega)
            end do
        end do
    end do
end subroutine

subroutine calc_H(H,sigma,phi) !geopotential at level p0*sigma
    real(kind=8), intent(inout) :: H(:,:,:)
    real(kind=8), intent(in)    :: sigma(:,:,:)
    real(kind=8), intent(in)    :: phi(:,:)

    integer(kind=4) :: i, j, k
    real(kind=8)    :: sigma_v, A, B, u

    do k = 1, size(H,3)
        do j = 1, size(H,2)
            do i = 1, size(H,1)
                sigma_v = 0.5*pi*(sigma(i,j,k)-sigma_0)
                u = u0*cos(sigma_v)**1.5_8
                A = (-2.0_8*sin(phi(i,j))**6*(cos(phi(i,j))**2+1.0_8/3.0_8)+10.0_8/63.0_8)
                B = (8.0_8/5.0_8*cos(phi(i,j))**3*(sin(phi(i,j))**2+2.0_8/3.0_8)-0.25_8*pi)

                H(i,j,k) = T0/Gamma*(1.0_8-sigma(i,j,k)**(Rd*Gamma/grav))+&
                                                       u*(A*u+B*Ra*omega)/grav
                if(sigma(i,j,k)<sigma_t) &
                    H(i,j,k) = H(i,j,k)-Rd*delta_T * &
                    ((log(sigma(i,j,k)/sigma_t)+137.0_8/60.0_8)*sigma_t**5            - &
                      5.0_8*sigma_t**4*sigma(i,j,k)+5.0_8*sigma_t**3*sigma(i,j,k)**2         - &
                      10.0_8/3.0_8*sigma_t**2*sigma(i,j,k)**3+1.25_8*sigma_t*sigma(i,j,k)**4 - &
                      0.2_8*sigma(i,j,k)**5)/grav
            end do
        end do
    end do
end subroutine

pure elemental real(kind=8) function calc_Pexner(sigma) result(P)
    real(kind=8), intent(in)    :: sigma

    P = sigma**kappa
end function

subroutine calc_theta(theta,sigma,phi)
    real(kind=8), intent(inout) :: theta(:,:,:)
    real(kind=8), intent(in)    :: sigma(:,:,:)
    real(kind=8), intent(in)    :: phi(:,:)

    integer(kind=4) :: ie, je, ke

    ie = size(theta,1); je = size(theta,2); ke = size(theta,3)

    call calc_Temp(theta,sigma,phi)
    theta(1:ie,1:je,1:ke) = theta(1:ie,1:je,1:ke) / calc_Pexner(sigma(1:ie,1:je,1:ke))
end subroutine

end module
