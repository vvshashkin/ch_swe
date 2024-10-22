module sbp_diff_c2i_42_mod

use sbp_diff_mod,   only : sbp_diff_t
use grid_field_mod, only : grid_field_t, tile_field_t
use mesh_mod,       only : tile_mesh_t, mesh_t
use tile_mod,       only : tile_t
use parcomm_mod,    only : parcomm_global


implicit none
! dot product matrix coeffs
real(kind=8), parameter :: h11 = 7._8/18._8, h22 = 9._8/8._8, h33 = 1._8, &
                           h44 = 71._8/72._8
! boundary projection operator coeffs
real(kind=8), parameter :: R11 = 15.0_8/8.0_8, R12 = -5.0_8/4.0_8, R13 = 3.0_8/8.0_8

real(kind=8), parameter :: q11 = -2.0_8, q12 = 3.0_8, q13 = -1.0_8, &
                           q21 = -1.0_8, q22 = 1.0_8,               &
                           q31 = 1._8/24._8, q32 = -9._8/8._8,  &
                           q33 = 9._8/8._8,  q34 = -1._8/24._8, &
                           q41 = -1._8/71._8, q42 = 6._8/71._8, q43 = -83._8/71._8, &
                           q44 = 81._8/71._8, q45 = -3._8/71._8

type, public, extends(sbp_diff_t) :: sbp_diff_c2i_42_t
contains
    procedure, public :: apply_tile
    procedure, public :: add_SAT_correction
end type sbp_diff_c2i_42_t

contains

subroutine apply_tile(this, f_out, f_in, mesh, direction)
    class(sbp_diff_c2i_42_t), intent(in)    :: this
    type(tile_field_t),       intent(inout) :: f_out
    type(tile_field_t),       intent(in)    :: f_in
    type(tile_mesh_t),        intent(in)    :: mesh
    character(len=*),         intent(in)    :: direction

    integer(kind=4) :: k, i, j

    select case(direction)
    case("x")
        do k = mesh%ks, mesh%ke
            do j = mesh%js, mesh%je
                if (mesh%is<=1 .and. mesh%ie>=1) then
                    f_out%p(1,j,k) = q11*f_in%p(1,j,k)+q12*f_in%p(2,j,k)+q13*f_in%p(3,j,k)
                end if
                if (mesh%is<=2 .and. mesh%ie>=2) then
                    f_out%p(2,j,k) = q21*f_in%p(1,j,k)+q22*f_in%p(2,j,k)
                end if
                if (mesh%is<=3 .and. mesh%ie>=3) then
                    f_out%p(3,j,k) = q31*f_in%p(1,j,k)+q32*f_in%p(2,j,k)+q33*f_in%p(3,j,k)+q34*f_in%p(4,j,k)
                end if
                if (mesh%is<=4 .and. mesh%ie>=4) then
                    f_out%p(4,j,k) = q41*f_in%p(1,j,k)+q42*f_in%p(2,j,k)+q43*f_in%p(3,j,k)+q44*f_in%p(4,j,k)+q45*f_in%p(5,j,k)
                end if
                do i = max(mesh%is, 5), min(mesh%ie, mesh%nx-4)
                    f_out%p(i,j,k) = (27.0_8*(f_in%p(i,j,k) - f_in%p(i-1,j,k)) &
                                      -(f_in%p(i+1,j,k) - f_in%p(i-2,j,k)))/24.0_8
                end do
                if (mesh%is<=mesh%nx-3 .and. mesh%ie>=mesh%nx-3) then
                    f_out%p(mesh%nx-3,j,k) = -q41*f_in%p(mesh%nx-1,j,k)-q42*f_in%p(mesh%nx-2,j,k)-q43*f_in%p(mesh%nx-3,j,k)-q44*f_in%p(mesh%nx-4,j,k)-q45*f_in%p(mesh%nx-5,j,k)
                end if
                if (mesh%is<=mesh%nx-2 .and. mesh%ie>=mesh%nx-2) then
                    f_out%p(mesh%nx-2,j,k) = -q31*f_in%p(mesh%nx-1,j,k)-q32*f_in%p(mesh%nx-2,j,k)-q33*f_in%p(mesh%nx-3,j,k)-q34*f_in%p(mesh%nx-4,j,k)
                end if
                if (mesh%is<=mesh%nx-1 .and. mesh%ie>=mesh%nx-1) then
                    f_out%p(mesh%nx-1,j,k) = -q21*f_in%p(mesh%nx-1,j,k)-q22*f_in%p(mesh%nx-2,j,k)
                end if
                if (mesh%is<=mesh%nx .and. mesh%ie>=mesh%nx) then
                    f_out%p(mesh%nx,j,k) = -q11*f_in%p(mesh%nx-1,j,k)-q12*f_in%p(mesh%nx-2,j,k)-q13*f_in%p(mesh%nx-3,j,k)
                end if
            end do
        end do
    case("y")
        do k = mesh%ks, mesh%ke
            if (mesh%js<=1 .and. mesh%je>=1) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,1,k) = q11*f_in%p(i,1,k)+q12*f_in%p(i,2,k)+q13*f_in%p(i,3,k)
                end do
            end if
            if (mesh%js<=2 .and. mesh%je>=2) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,2,k) = q21*f_in%p(i,1,k)+q22*f_in%p(i,2,k)
                end do
            end if
            if (mesh%js<=3 .and. mesh%je>=3) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,3,k) = q31*f_in%p(i,1,k)+q32*f_in%p(i,2,k)+q33*f_in%p(i,3,k)+q34*f_in%p(i,4,k)
                end do
            end if
            if (mesh%js<=4 .and. mesh%je>=4) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,4,k) = q41*f_in%p(i,1,k)+q42*f_in%p(i,2,k)+q43*f_in%p(i,3,k)+q44*f_in%p(i,4,k)+q45*f_in%p(i,5,k)
                end do
            end if
            do j = max(mesh%js, 5), min(mesh%je, mesh%ny-4)
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,k) = (27.0_8*(f_in%p(i,j,k) - f_in%p(i,j-1,k)) &
                                      -(f_in%p(i,j+1,k) - f_in%p(i,j-2,k)))/24.0_8
                end do
            end do
            if (mesh%js<=mesh%ny-3 .and. mesh%je>=mesh%ny-3) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-3,k) = -q41*f_in%p(i,mesh%ny-1,k)-q42*f_in%p(i,mesh%ny-2,k)-q43*f_in%p(i,mesh%ny-3,k)-q44*f_in%p(i,mesh%ny-4,k)-q45*f_in%p(i,mesh%ny-5,k)
                end do
            end if
            if (mesh%js<=mesh%ny-2 .and. mesh%je>=mesh%ny-2) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-2,k) = -q31*f_in%p(i,mesh%ny-1,k)-q32*f_in%p(i,mesh%ny-2,k)-q33*f_in%p(i,mesh%ny-3,k)-q34*f_in%p(i,mesh%ny-4,k)
                end do
            end if
            if (mesh%js<=mesh%ny-1 .and. mesh%je>=mesh%ny-1) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-1,k) = -q21*f_in%p(i,mesh%ny-1,k)-q22*f_in%p(i,mesh%ny-2,k)
                end do
            end if
            if (mesh%js<=mesh%ny .and. mesh%je>=mesh%ny) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny,k) = -q11*f_in%p(i,mesh%ny-1,k)-q12*f_in%p(i,mesh%ny-2,k)-q13*f_in%p(i,mesh%ny-3,k)
                end do
            end if
        end do
    case("z")
        if (mesh%ks<=1 .and. mesh%ke>=1) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,1) = q11*f_in%p(i,j,1)+q12*f_in%p(i,j,2)+q13*f_in%p(i,j,3)
                end do
            end do
        end if
        if (mesh%ks<=2 .and. mesh%ke>=2) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,2) = q21*f_in%p(i,j,1)+q22*f_in%p(i,j,2)
                end do
            end do
        end if
        if (mesh%ks<=3 .and. mesh%ke>=3) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,3) = q31*f_in%p(i,j,1)+q32*f_in%p(i,j,2)+q33*f_in%p(i,j,3)+q34*f_in%p(i,j,4)
                end do
            end do
        end if
        if (mesh%ks<=4 .and. mesh%ke>=4) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,4) = q41*f_in%p(i,j,1)+q42*f_in%p(i,j,2)+q43*f_in%p(i,j,3)+q44*f_in%p(i,j,4)+q45*f_in%p(i,j,5)
                end do
            end do
        end if
        do k = max(mesh%ks, 5), min(mesh%ke, mesh%nz-4)
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,k) = (27.0_8*(f_in%p(i,j,k) - f_in%p(i,j,k-1)) &
                                      -(f_in%p(i,j,k+1) - f_in%p(i,j,k-2)))/24.0_8
                end do
            end do
        end do
        if (mesh%ks<=mesh%nz-3 .and. mesh%ke>=mesh%nz-3) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-3) = -q41*f_in%p(i,j,mesh%nz-1)-q42*f_in%p(i,j,mesh%nz-2)-q43*f_in%p(i,j,mesh%nz-3)-q44*f_in%p(i,j,mesh%nz-4)-q45*f_in%p(i,j,mesh%nz-5)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz-2 .and. mesh%ke>=mesh%nz-2) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-2) = -q31*f_in%p(i,j,mesh%nz-1)-q32*f_in%p(i,j,mesh%nz-2)-q33*f_in%p(i,j,mesh%nz-3)-q34*f_in%p(i,j,mesh%nz-4)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz-1 .and. mesh%ke>=mesh%nz-1) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-1) = -q21*f_in%p(i,j,mesh%nz-1)-q22*f_in%p(i,j,mesh%nz-2)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz .and. mesh%ke>=mesh%nz) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz) = -q11*f_in%p(i,j,mesh%nz-1)-q12*f_in%p(i,j,mesh%nz-2)-q13*f_in%p(i,j,mesh%nz-3)
                end do
            end do
        end if
    case default
        call parcomm_global%abort("Wrong direction in sbp_diff_c2i_42_mod")
    end select

end subroutine apply_tile

subroutine add_SAT_correction(this, f_out, f_in, mesh, direction)
    class(sbp_diff_c2i_42_t), intent(in)    :: this
    type(tile_field_t),       intent(inout) :: f_out
    type(tile_field_t),       intent(in)    :: f_in
    type(tile_mesh_t),        intent(in)    :: mesh
    character(len=*),         intent(in)    :: direction

    integer(kind=4) :: k, i, j

    select case(direction)
    case("x")
        if (mesh%is<=1 .and. mesh%ie>=1) then
            do k = mesh%ks, mesh%ke
                do j = mesh%js, mesh%je
                    f_out%p(1,j,k) = f_out%p(1,j,k) &
                                   + 0.5_8/h11*(R11*(f_in%p(1,j,k)-f_in%p( 0,j,k)) + &
                                                R12*(f_in%p(2,j,k)-f_in%p(-1,j,k)) + &
                                                R13*(f_in%p(3,j,k)-f_in%p(-2,j,k)))
                end do
            end do
        end if
        if (mesh%is<=mesh%nx .and. mesh%ie>=mesh%nx) then
            do k = mesh%ks, mesh%ke
                do j = mesh%js, mesh%je
                    f_out%p(mesh%nx,j,k) = f_out%p(mesh%nx,j,k) &
                                   - 0.5_8/h11*(R11*(f_in%p(mesh%nx-1,j,k)-f_in%p(mesh%nx  ,j,k)) + &
                                                R12*(f_in%p(mesh%nx-2,j,k)-f_in%p(mesh%nx+1,j,k)) + &
                                                R13*(f_in%p(mesh%nx-3,j,k)-f_in%p(mesh%nx+2,j,k)))
                end do
            end do
        end if
    case("y")
        if (mesh%js<=1 .and. mesh%je>=1) then
            do k = mesh%ks, mesh%ke
                do i = mesh%is, mesh%ie
                    f_out%p(i,1,k) = f_out%p(i,1,k) &
                                   + 0.5_8/h11*(R11*(f_in%p(i,1,k)-f_in%p(i, 0,k)) + &
                                                R12*(f_in%p(i,2,k)-f_in%p(i,-1,k)) + &
                                                R13*(f_in%p(i,3,k)-f_in%p(i,-2,k)))
                end do
            end do
        end if
        if (mesh%js<=mesh%ny .and. mesh%je>=mesh%ny) then
            do k = mesh%ks, mesh%ke
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny,k) = f_out%p(i,mesh%ny,k) &
                                   - 0.5_8/h11*(R11*(f_in%p(i,mesh%ny-1,k)-f_in%p(i,mesh%ny  ,k)) + &
                                                R12*(f_in%p(i,mesh%ny-2,k)-f_in%p(i,mesh%ny+1,k)) + &
                                                R13*(f_in%p(i,mesh%ny-3,k)-f_in%p(i,mesh%ny+2,k)))
                end do
            end do
        end if
    case default
        call parcomm_global%abort("Wrong direction in add_SAT_correction in sbp_diff_c2i_42_mod.f90")
    end select

end subroutine add_SAT_correction

end module sbp_diff_c2i_42_mod
