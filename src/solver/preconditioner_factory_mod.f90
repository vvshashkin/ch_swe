module preconditioner_factory_mod

use preconditioner_mod,          only : preconditioner_t
use grid_field_based_vector_mod, only : grid_field_based_vector_t
use linear_operator_mod,         only : linear_operator_t
use domain_mod,                  only : domain_t
use abstract_quadrature_mod,     only : quadrature_t

implicit none

contains

subroutine create_preconditioner(preconditioner, linear_operator, precond_name, domain)

    class(preconditioner_t), allocatable, intent(out)   :: preconditioner
    class(linear_operator_t),             intent(inout) :: linear_operator
    character(len=*),                     intent(in)    :: precond_name
    type(domain_t),                       intent(in)    :: domain

    select case(precond_name)
    case("jacobi")
        call create_grid_field_based_jacobi_preconditioner(preconditioner, linear_operator, domain)
    case("none")

    case default
        call domain%parcomm%abort("Unknown preconditioner name: "//precond_name)
    end select

end subroutine create_preconditioner

subroutine create_grid_field_based_jacobi_preconditioner(preconditioner, linear_operator, domain)

    use diagonal_preconditioner_mod, only : diagonal_preconditioner_t

    class(preconditioner_t), allocatable, intent(out)   :: preconditioner
    class(linear_operator_t),             intent(inout) :: linear_operator
    type(domain_t),                       intent(in)    :: domain

    type(diagonal_preconditioner_t), allocatable :: diag_prec

    type(grid_field_based_vector_t) :: tmp, unit_vector, inv_diag

    real(kind=8) :: sum

    integer(kind=4) :: t, k, j, i, ii, jj, kk, tt, ind

    call         tmp%init(8, 0, domain%mesh_p, quadrature_t())
    call    inv_diag%init(0, 0, domain%mesh_p, quadrature_t())
    call unit_vector%init(8, 0, domain%mesh_p, quadrature_t())

    call unit_vector%set_scalar(0.0_8, domain)

    ! open(123123 , file = 'mtrx.dat',access="direct",recl=1)

    ! ind = 1

    do t = 1, domain%partition%Nt
        do k = 1, domain%partition%tiles_p%Nz
            do j = 1, domain%partition%tiles_p%Ny
                do i = 1, domain%partition%tiles_p%Nx
                    if (domain%mesh_p%ts <= t .and. t <= domain%mesh_p%te) then
                    if (domain%mesh_p%tile(t)%ks <= k .and. k <= domain%mesh_p%tile(t)%ke) then
                    if (domain%mesh_p%tile(t)%js <= j .and. j <= domain%mesh_p%tile(t)%je) then
                    if (domain%mesh_p%tile(t)%is <= i .and. i <= domain%mesh_p%tile(t)%ie) then
                        unit_vector%grid_field%tile(t)%p(i, j, k) = 1.0_8
                    end if
                    end if
                    end if
                    end if

                    call linear_operator%apply(tmp, unit_vector, domain)

                    ! do tt = domain%mesh_p%ts, domain%mesh_p%te
                    !     do kk = domain%mesh_p%tile(t)%ks, domain%mesh_p%tile(t)%ke
                    !         do jj = domain%mesh_p%tile(t)%js, domain%mesh_p%tile(t)%je
                    !             do ii = domain%mesh_p%tile(t)%is, domain%mesh_p%tile(t)%ie
                    !                 print*, ind
                    !                 write(123123,rec=ind) real(tmp%grid_field%tile(tt)%p(ii, jj, kk),4)
                    !                 ind = ind+1
                    !             end do
                    !         end do
                    !     end do
                    ! end do

                    if (domain%mesh_p%ts <= t .and. t <= domain%mesh_p%te) then
                    if (domain%mesh_p%tile(t)%ks <= k .and. k <= domain%mesh_p%tile(t)%ke) then
                    if (domain%mesh_p%tile(t)%js <= j .and. j <= domain%mesh_p%tile(t)%je) then
                    if (domain%mesh_p%tile(t)%is <= i .and. i <= domain%mesh_p%tile(t)%ie) then
                        inv_diag%grid_field%tile(t)%p(i, j, k) = 1.0_8/tmp%grid_field%tile(t)%p(i, j, k)
                        unit_vector%grid_field%tile(t)%p(i, j, k) = 0.0_8
                    end if
                    end if
                    end if
                    end if
                end do
            end do
        end do
    end do

    allocate(diag_prec)

    call inv_diag%create_similar(diag_prec%diag, domain)

    call diag_prec%diag%copy(inv_diag, domain)

    call move_alloc(diag_prec, preconditioner)
! close(123123)
end subroutine create_grid_field_based_jacobi_preconditioner

end module preconditioner_factory_mod
