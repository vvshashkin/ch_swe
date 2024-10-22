module sbp_interp_c2i_21_mod

use sbp_interp_mod, only : sbp_interp_t
use grid_field_mod, only : grid_field_t, tile_field_t
use mesh_mod,       only : tile_mesh_t, mesh_t
use tile_mod,       only : tile_t
use parcomm_mod,    only : parcomm_global


implicit none

type, public, extends(sbp_interp_t) :: sbp_interp_c2i_21_t
contains
    procedure, public :: apply_tile
end type sbp_interp_c2i_21_t

contains

subroutine apply_tile(this, f_out, f_in, mesh, direction)
    class(sbp_interp_c2i_21_t), intent(in)    :: this
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
                    f_out%p(1,j,k) = f_in%p(1,j,k)
                end if
                do i = max(mesh%is, 2), min(mesh%ie, mesh%nx-1)
                    f_out%p(i,j,k) = 0.5_8*(f_in%p(i,j,k) + f_in%p(i-1,j,k))
                end do
                if (mesh%is<=mesh%nx .and. mesh%ie>=mesh%nx) then
                    f_out%p(mesh%nx,j,k) = f_in%p(mesh%nx-1,j,k)
                end if
            end do
        end do
    case("y")
        do k = mesh%ks, mesh%ke
            if (mesh%js<=1 .and. mesh%je>=1) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,1,k) = f_in%p(i,1,k)
                end do
            end if
            do j = max(mesh%js, 2), min(mesh%je, mesh%ny-1)
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,k) =  0.5_8*(f_in%p(i,j,k) + f_in%p(i,j-1,k))
                end do
            end do
            if (mesh%js<=mesh%ny .and. mesh%je>=mesh%ny) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny,k) = f_in%p(i,mesh%ny-1,k)
                end do
            end if
        end do
    case("z")
        if (mesh%ks<=1 .and. mesh%ke>=1) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,1) = f_in%p(i,j,1)
                end do
            end do
        end if
        do k = max(mesh%ks, 2), min(mesh%ke, mesh%nz-1)
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,k) =  0.5_8*(f_in%p(i,j,k) + f_in%p(i,j,k-1))
                end do
            end do
        end do
        if (mesh%ks<=mesh%nz .and. mesh%ke>=mesh%nz) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz) = f_in%p(i,j,mesh%nz-1)
                end do
            end do
        end if
    case default
        call parcomm_global%abort("Wrong direction in sbp_interp_c2i_21_mod")
    end select

end subroutine apply_tile

end module sbp_interp_c2i_21_mod
