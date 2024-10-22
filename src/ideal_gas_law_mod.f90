module ideal_gas_law_mod

use grid_field_mod, only : grid_field_t, tile_field_t
use mesh_mod,       only : mesh_t, tile_mesh_t
use const_mod,      only : Rgaz, Cp, Cv, kappa, P_reference

implicit none

contains

subroutine calc_PExner_from_p(P, pres, pref, mesh)

    type(grid_field_t), intent(inout) :: P
    type(grid_field_t), intent(in)    :: pres
    real(kind=8),       intent(in)    :: pref
    type(mesh_t),       intent(in)    :: mesh

    integer(kind=4) :: i, j, k, t

    do t = mesh%ts, mesh%te
        do k = mesh%tile(t)%ks, mesh%tile(t)%ke
            do j = mesh%tile(t)%js, mesh%tile(t)%je
                do i = mesh%tile(t)%is, mesh%tile(t)%ie
                    P%tile(t)%p(i,j,k) = (pres%tile(t)%p(i,j,k)/P_reference)**kappa
                end do
            end do
        end do
    end do

end subroutine

subroutine calc_PExner(P, rho, theta, mesh)
    type(grid_field_t), intent(inout) :: P
    type(grid_field_t), intent(in)    :: rho, theta
    type(mesh_t),       intent(in)    :: mesh

    integer(kind=4) :: i, j, k, t

    do t = mesh%ts, mesh%te
        do k = mesh%tile(t)%ks, mesh%tile(t)%ke
            do j = mesh%tile(t)%js, mesh%tile(t)%je
                do i = mesh%tile(t)%is, mesh%tile(t)%ie
                    P%tile(t)%p(i,j,k) = (Rgaz*rho%tile(t)%p(i,j,k)*theta%tile(t)%p(i,j,k)/P_reference)**(Rgaz/Cv)
                end do
            end do
        end do
    end do

end subroutine

subroutine calc_PExner_prime(P, rho0, rho1, theta0, theta1, P0, mesh)
    type(grid_field_t), intent(inout) :: P
    type(grid_field_t), intent(in)    :: rho0, rho1, theta0, theta1, P0
    type(mesh_t),       intent(in)    :: mesh

    integer(kind=4) :: i, j, k, t

    do t = mesh%ts, mesh%te
        do k = mesh%tile(t)%ks, mesh%tile(t)%ke
            do j = mesh%tile(t)%js, mesh%tile(t)%je
                do i = mesh%tile(t)%is, mesh%tile(t)%ie
                    P%tile(t)%p(i,j,k) = Rgaz/Cv *P0%tile(t)%p(i,j,k) * &
                                         (rho1%tile(t)%p(i,j,k)   / rho0%tile(t)%p(i,j,k)+  &
                                          theta1%tile(t)%p(i,j,k) / theta0%tile(t)%p(i,j,k))
                end do
            end do
        end do
    end do

end subroutine

end module
