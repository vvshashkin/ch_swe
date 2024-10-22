module sbp_interp_i2c_21_mod

use sbp_interp_mod, only : sbp_interp_t
use grid_field_mod, only : grid_field_t, tile_field_t
use mesh_mod,       only : tile_mesh_t, mesh_t
use tile_mod,       only : tile_t
use parcomm_mod,    only : parcomm_global


implicit none

type, public, extends(sbp_interp_t) :: sbp_interp_i2c_21_t
contains
    procedure, public :: apply_tile
end type sbp_interp_i2c_21_t

contains

subroutine apply_tile(this, f_out, f_in, mesh, direction)
    class(sbp_interp_i2c_21_t), intent(in)    :: this
    type(tile_field_t),         intent(inout) :: f_out
    type(tile_field_t),         intent(in)    :: f_in
    type(tile_mesh_t),          intent(in)    :: mesh
    character(len=*),           intent(in)    :: direction

    integer(kind=4) :: k, i, j

    select case(direction)
    case("x")
        do k = mesh%ks, mesh%ke
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,k) = 0.5_8*(f_in%p(i+1,j,k) + f_in%p(i,j,k))
                end do
            end do
        end do
    case("y")
        do k = mesh%ks, mesh%ke
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,k) = 0.5_8*(f_in%p(i,j+1,k) + f_in%p(i,j,k))
                end do
            end do
        end do
    case("z")
        do k = mesh%ks, mesh%ke
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,k) = 0.5_8*(f_in%p(i,j,k+1) + f_in%p(i,j,k))
                end do
            end do
        end do
    case default
        call parcomm_global%abort("Wrong direction in sbp_interp_i2c_21_mod")
    end select

end subroutine apply_tile

end module sbp_interp_i2c_21_mod
