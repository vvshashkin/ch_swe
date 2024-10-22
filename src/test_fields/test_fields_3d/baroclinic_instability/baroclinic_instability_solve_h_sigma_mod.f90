module baroclinic_instability_solve_h_sigma_mod

use baroclinic_instability_test_parameters_mod, only : Gamma, T0
use baroclinic_instability_test_therm_mod,      only : calc_H, calc_temp
use const_mod, only : grav=>Earth_grav, Rd=>rgaz
implicit none

real(kind=8),    parameter :: tolerance = 1e-10
integer(kind=4), parameter :: maxit = 15

contains

subroutine solve_h_for_sigma(sigma,h,phi)
    real(kind=8), intent(inout) :: sigma(:,:,:)
    real(kind=8), intent(in)    :: h(:,:,:)
    real(kind=8), intent(in)    :: phi(:,:)

    integer(kind=4) :: i, j, k, it, ie, je, ke
    real(kind=8) :: temp(1:size(sigma,1),1:size(sigma,2),1:size(sigma,3))
    real(kind=8) :: dh(1:size(sigma,1),1:size(sigma,2),1:size(sigma,3))

    ie = size(sigma,1); je = size(sigma,2); ke = size(sigma,3)
    !first guess
    do k = 1, ke
        do j = 1, je
            do i = 1, ie
                sigma(i,j,k) = min(1.0_8,max(max(1.0_8-Gamma*h(i,j,k)/T0,0.0_8)**(grav/(Rd*Gamma)),1e-5_8))
            end do
        end do
    enddo
    do it = 1, maxit
        call calc_Temp(temp,sigma,phi)
        call calc_H(dh,sigma,phi)
        do k=1,ke; do j=1,je; do i=1,ie
            dh(i,j,k) = h(i,j,k)-dh(i,j,k)
            sigma(i,j,k) = min(1.0_8,max(1e-6,sigma(i,j,k)*(1.0_8-dh(i,j,k)*grav/(Rd*temp(i,j,k)))))
        end do; end do; end do
        if(maxval(abs(dh(:,:,:))) < tolerance) exit
    end do
end subroutine

end module
