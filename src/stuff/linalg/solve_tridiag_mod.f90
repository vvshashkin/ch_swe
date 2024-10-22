module solve_tridiag_mod

use grid_field_mod, only : grid_field_t, tile_field_t
use mesh_mod,       only : mesh_t, tile_mesh_t

implicit none

contains

!solves a*sol(k-1)+b*sol(k)+c*sol(k+1) = rhs(k) in last index for 3rd-rank vectors
!a(1) assumed 0, so a is (:,:,2:) -shaped array

subroutine solve_tridiag(sol,rhs,a,b,c)
    real(kind=8), intent(out)   :: sol(:,:,:)
    real(kind=8), intent(in)    :: rhs(:,:,:)
    real(kind=8), intent(inout) :: a(:,:,2:), b(:,:,:), c(:,:,:)

    real(kind=8)    :: b1(size(rhs,1),size(rhs,2),size(rhs,3))
    integer(kind=4) :: i, j, k, nx, ny, nz

    nx = size(sol,1)
    ny = size(sol,2)
    nz = size(sol,3)

    sol(1:nx,1:ny,1)  = rhs(1:nx,1:ny,1)
    b1(1:nx,1:ny,1) = b(1:nx,1:ny,1)

    do k = 2, nz
        do j = 1, ny
            do i = 1, nx
                b1(i,j,k) = b(i,j,k)-a(i,j,k) / b1(i,j,k-1) * c(i,j,k-1)
                sol(i,j,k)  = rhs(i,j,k)-a(i,j,k) / b1(i,j,k-1) * sol(i,j,k-1)
            end do
        end do
    end do

    do j = 1, ny
        do i = 1, nx
            sol(i,j,nz) = sol(i,j,nz) / b1(i,j,nz)
        end do
    end do

    do k = nz-1,1,-1
        do j = 1, ny
            do i = 1, nx
                sol(i,j,k) = (sol(i,j,k)-c(i,j,k)*sol(i,j,k+1)) / b1(i,j,k)
            end do
        end do
    end do

end subroutine

!Factored tridiag to accelerate multiple solves with the same matrix
!solves a_k*sol(k-1) + b_k*sol(k) + c_k*sol(k+1) = rhs(k)
!in vertical for grid_fields
!
!Fisrt step (factor_tridiag) is some kind of LU factorization with a-b-c matrix:
!Second step (solve_factored) is solution of A*sol=rhs with previously prepared a-b-c coefficients

subroutine factor_tridiag(a,b,c,mesh) !prepare coefficients for solution
    type(grid_field_t), intent(inout) :: a, b, c
    type(mesh_t),       intent(in)    :: mesh

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        call factor_tridiag_tile(a%tile(t),b%tile(t), c%tile(t), mesh%tile(t))
    end do
end subroutine

subroutine factor_tridiag_tile(a,b,c,mesh)
    type(tile_field_t), intent(inout) :: a, b, c
    type(tile_mesh_t),  intent(in)    :: mesh

    integer(kind=4) :: i, j, k

    do k = 2, mesh%nz
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                a%p(i,j,k) = a%p(i,j,k) / b%p(i,j,k-1)
                b%p(i,j,k) = b%p(i,j,k) - c%p(i,j,k-1)*a%p(i,j,k)
            end do
        end do
    end do

    do k = 1, mesh%nz
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                b%p(i,j,k) = 1.0_8 / b%p(i,j,k)
            end do
        end do
    end do

end subroutine

subroutine solve_factored_tridiag(sol,rhs,a,b,c,mesh) !use prepared coefficients for solution
    type(grid_field_t), intent(inout) :: sol
    type(grid_field_t), intent(in)    :: rhs, a, b, c
    type(mesh_t),       intent(in)    :: mesh

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        call solve_factored_tridiag_tile(sol%tile(t), rhs%tile(t), a%tile(t), &
                                         b%tile(t), c%tile(t), mesh%tile(t))
    end do
end subroutine

subroutine solve_factored_tridiag_tile(sol,rhs,a,b,c,mesh)
    type(tile_field_t), intent(inout) :: sol
    type(tile_field_t), intent(in)    :: rhs, a, b, c
    type(tile_mesh_t),  intent(in)    :: mesh

    integer(kind=4) :: i, j, k

    sol%p(mesh%is:mesh%ie,mesh%js:mesh%je,1) = rhs%p(mesh%is:mesh%ie,mesh%js:mesh%je,1)

    do k = 2, mesh%nz
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                sol%p(i,j,k) = rhs%p(i,j,k) - a%p(i,j,k) * sol%p(i,j,k-1)
            end do
        end do
    end do

    do j = mesh%js, mesh%je
        do i = mesh%is, mesh%ie
            sol%p(i,j,mesh%nz) = sol%p(i,j,mesh%nz) * b%p(i,j,mesh%nz)
        end do
    end do

    do k = mesh%nz-1,1,-1
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                sol%p(i,j,k) = b%p(i,j,k)*(sol%p(i,j,k)-c%p(i,j,k)*sol%p(i,j,k+1))
            end do
        end do
    end do

end subroutine

end module