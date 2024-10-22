module sbp_diff_i2c_21_mod

use sbp_diff_mod,   only : sbp_diff_t
use grid_field_mod, only : grid_field_t, tile_field_t
use mesh_mod,       only : tile_mesh_t, mesh_t
use tile_mod,       only : tile_t
use parcomm_mod,    only : parcomm_global


implicit none
! dot product matrix coeffs
real(kind=8), parameter :: h11 = 1.0_8, h22 = 1.0_8
! boundary projection operator coeffs
real(kind=8), parameter :: R11 = 3.0_8 / 2.0_8, R12 = -1.0_8 / 2.0_8

type, public, extends(sbp_diff_t) :: sbp_diff_i2c_21_t
contains
    procedure, public :: apply_tile
    procedure, public :: add_SAT_correction
end type sbp_diff_i2c_21_t

contains

subroutine apply_tile(this, f_out, f_in, mesh, direction)
    class(sbp_diff_i2c_21_t), intent(in)    :: this
    type(tile_field_t),       intent(inout) :: f_out
    type(tile_field_t),       intent(in)    :: f_in
    type(tile_mesh_t),        intent(in)    :: mesh
    character(len=*),         intent(in)    :: direction

    integer(kind=4) :: k, i, j

    select case(direction)
    case("x")
        do k = mesh%ks, mesh%ke
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,k) = f_in%p(i+1,j,k) - f_in%p(i,j,k)
                end do
            end do
        end do
    case("y")
        do k = mesh%ks, mesh%ke
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,k) = f_in%p(i,j+1,k) - f_in%p(i,j,k)
                end do
            end do
        end do
    case("z")
        do k = mesh%ks, mesh%ke
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,k) = f_in%p(i,j,k+1) - f_in%p(i,j,k)
                end do
            end do
        end do
    case default
        call parcomm_global%abort("Wrong direction in sbp_diff_i2c_21_mod")
    end select

end subroutine apply_tile

subroutine add_SAT_correction(this, f_out, f_in, mesh, direction)
    class(sbp_diff_i2c_21_t), intent(in)    :: this
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
                                   + 0.5_8*R11/h11*(f_in%p(1,j,k)-f_in%p(0,j,k))
                end do
            end do
        end if
        if (mesh%is<=2 .and. mesh%ie>=2) then
            do k = mesh%ks, mesh%ke
                do j = mesh%js, mesh%je
                    f_out%p(2,j,k) = f_out%p(2,j,k) &
                                   + 0.5_8*R12/h22*(f_in%p(1,j,k)-f_in%p(0,j,k))
                end do
            end do
        end if
        if (mesh%is<=mesh%nx .and. mesh%ie>=mesh%nx) then
            do k = mesh%ks, mesh%ke
                do j = mesh%js, mesh%je
                    f_out%p(mesh%nx,j,k) = f_out%p(mesh%nx,j,k) &
                                   + 0.5_8*R11/h11*(f_in%p(mesh%nx+2,j,k)-f_in%p(mesh%nx+1,j,k))
                end do
            end do
        end if
        if (mesh%is<=mesh%nx-1 .and. mesh%ie>=mesh%nx-1) then
            do k = mesh%ks, mesh%ke
                do j = mesh%js, mesh%je
                    f_out%p(mesh%nx-1,j,k) = f_out%p(mesh%nx-1,j,k) &
                                   + 0.5_8*R12/h22*(f_in%p(mesh%nx+2,j,k)-f_in%p(mesh%nx+1,j,k))
                end do
            end do
        end if
    case("y")
        if (mesh%js<=1 .and. mesh%je>=1) then
            do k = mesh%ks, mesh%ke
                do i = mesh%is, mesh%ie
                    f_out%p(i,1,k) = f_out%p(i,1,k) &
                                   + 0.5_8*R11/h11*(f_in%p(i,1,k)-f_in%p(i,0,k))
                end do
            end do
        end if
        if (mesh%js<=2 .and. mesh%je>=2) then
            do k = mesh%ks, mesh%ke
                do i = mesh%is, mesh%ie
                    f_out%p(i,2,k) = f_out%p(i,2,k) &
                                   + 0.5_8*R12/h22*(f_in%p(i,1,k)-f_in%p(i,0,k))
                end do
            end do
        end if
        if (mesh%js<=mesh%ny .and. mesh%je>=mesh%ny) then
            do k = mesh%ks, mesh%ke
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny,k) = f_out%p(i,mesh%ny,k) &
                                   + 0.5_8*R11/h11*(f_in%p(i,mesh%ny+2,k)-f_in%p(i,mesh%ny+1,k))
                end do
            end do
        end if
        if (mesh%js<=mesh%ny-1 .and. mesh%je>=mesh%ny-1) then
            do k = mesh%ks, mesh%ke
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-1,k) = f_out%p(i,mesh%ny-1,k) &
                                   + 0.5_8*R12/h22*(f_in%p(i,mesh%ny+2,k)-f_in%p(i,mesh%ny+1,k))
                end do
            end do
        end if
    case default
        call parcomm_global%abort("Wrong direction in add_interface_SAT_correction in sbp_diff_i2c_21_mod.f90")
    end select

end subroutine add_SAT_correction

end module sbp_diff_i2c_21_mod
