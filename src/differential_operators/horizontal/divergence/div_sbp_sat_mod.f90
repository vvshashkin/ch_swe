module div_sbp_SAT_mod

use abstract_div_mod,   only : div_operator_t
use domain_mod,         only : domain_t
use mesh_mod,           only : tile_mesh_t
use grid_field_mod,     only : grid_field_t, tile_field_t
use exchange_halo_mod,  only : exchange_t
use mesh_mod,           only : mesh_t
use sbp_diff_mod,       only : sbp_diff_t
use vec_math_mod,       only : multiply_by_J
use halo_mod,           only : halo_t

implicit none

type, public, extends(div_operator_t) :: div_sbp_SAT_t
    type(grid_field_t)              :: Ju, Jv
    type(grid_field_t)              :: Dx, Dy
    class(exchange_t), allocatable  :: exch_halo
    class(sbp_diff_t), allocatable  :: sbp_diff
    type(mesh_t),      pointer      :: mesh_u, mesh_v, mesh_p
    class(halo_t),     allocatable  :: proj_op
contains
    procedure, public :: calc_div
    procedure, public :: calc_div_tile
end type div_sbp_SAT_t

contains

subroutine calc_div(this, div, u, v, domain)

    class(div_sbp_SAT_t), intent(inout) :: this
    type(domain_t),       intent(in)    :: domain
    type(grid_field_t),   intent(inout) :: u, v
    type(grid_field_t),   intent(inout) :: div

    integer(kind=4) :: i, j, k, t
    real(kind=8) hx

    call multiply_by_J(this%Ju, u, this%mesh_u)
    call multiply_by_J(this%Jv, v, this%mesh_v)

    call this%exch_halo%do_vec(this%Ju, this%Jv, domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call this%calc_div_tile(div%tile(t),                      &
                                this%Ju%tile(t), this%Jv%tile(t), &
                                this%Dx%tile(t), this%Dy%tile(t), &
                                this%mesh_p%tile(t), this%mesh_p%scale)
    end do

    if(allocated(this%proj_op)) call this%proj_op%get_halo_scalar(div, domain, 1)
end subroutine calc_div

subroutine calc_div_tile(this, div, Ju, Jv, Dx, Dy, mesh, scale)

    class(div_sbp_SAT_t), intent(inout) :: this
    type(tile_field_t),   intent(inout) :: div
    type(tile_field_t),   intent(in)    :: Ju, Jv
    type(tile_field_t),   intent(inout) :: Dx, Dy
    type(tile_mesh_t),    intent(in)    :: mesh
    real(kind=8),         intent(in)    :: scale

    integer(kind=4) :: i, j, k

    call this%sbp_diff%apply_tile(Dx, Ju, mesh, "x")
    call this%sbp_diff%add_SAT_correction(Dx, Ju, mesh, "x")

    call this%sbp_diff%apply_tile(Dy, Jv, mesh, "y")
    call this%sbp_diff%add_SAT_correction(Dy, Jv, mesh, "y")

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                div%p(i,j,k) = (Dx%p(i,j,k) / mesh%hx + Dy%p(i,j,k)/ mesh%hy) &
                             / (mesh%J(i,j,k)*scale)
            end do
        end do
    end do

end subroutine calc_div_tile

end module div_sbp_sat_mod
