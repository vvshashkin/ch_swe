module grad_perp_sbp_sat_mod

use domain_mod,             only : domain_t
use abstract_grad_perp_mod, only : grad_perp_operator_t
use grid_field_mod,         only : grid_field_t, tile_field_t
use exchange_abstract_mod,  only : exchange_t
use sbp_diff_mod,           only : sbp_diff_t
use mesh_mod,               only : mesh_t

implicit none

private

type, public, extends(grad_perp_operator_t) :: grad_perp_sbp_sat_t
    class(exchange_t), allocatable :: exch_halo
    class(sbp_diff_t), allocatable :: sbp_diff
    type(mesh_t),      pointer     :: mesh_u, mesh_v
contains
    procedure, public :: calc_grad_perp => calc_grad_perp_sbp_sat
end type grad_perp_sbp_sat_t

contains

subroutine calc_grad_perp_sbp_sat(this, gu, gv, w, domain)
    class(grad_perp_sbp_sat_t), intent(inout) :: this
    type(domain_t),             intent(in)    :: domain
    type(grid_field_t),         intent(inout) :: w
    type(grid_field_t),         intent(inout) :: gu, gv

    integer(kind=4) :: t

    call this%exch_halo%do(w, domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call calc_grad_perp_on_tile(gu%tile(t), gv%tile(t), w%tile(t),    &
                               this%mesh_u%tile(t), this%mesh_v%tile(t), this%sbp_diff,&
                               domain%mesh_u%scale)
    end do

end subroutine calc_grad_perp_sbp_sat

subroutine calc_grad_perp_on_tile(gx, gy, f, mesh_u, mesh_v, sbp_diff, scale)

    use mesh_mod, only : tile_mesh_t

    type(tile_field_t), intent(inout) :: gx, gy
    type(tile_field_t), intent(in)    :: f
    type(tile_mesh_t),  intent(in)    :: mesh_u, mesh_v
    class(sbp_diff_t),  intent(in)    :: sbp_diff
    real(kind=8),       intent(in)    :: scale

    integer(kind=4) :: i, j, k

    call sbp_diff%apply_tile(gx, f, mesh_u, "y")
    call sbp_diff%add_SAT_correction(gx, f, mesh_u, "y")

    do k = mesh_u%ks, mesh_u%ke
        do j = mesh_u%js, mesh_u%je
            do i = mesh_u%is, mesh_u%ie
                gx%p(i,j,k) = - gx%p(i,j,k) / (mesh_u%hy * scale * mesh_u%J(i,j,k))
            end do
        end do
    end do

    call sbp_diff%apply_tile(gy, f, mesh_v, "x")
    call sbp_diff%add_SAT_correction(gy, f, mesh_v, "x")

    do k = mesh_v%ks, mesh_v%ke
        do j = mesh_v%js, mesh_v%je
            do i = mesh_v%is, mesh_v%ie
                gy%p(i,j,k) = gy%p(i,j,k) / (mesh_v%hx*scale * mesh_v%J(i,j,k))
            end do
        end do
    end do

end subroutine calc_grad_perp_on_tile

end module grad_perp_sbp_sat_mod
