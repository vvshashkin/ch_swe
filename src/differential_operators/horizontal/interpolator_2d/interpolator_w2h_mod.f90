module interpolator_w2h_mod

use grid_field_mod,               only : grid_field_t, tile_field_t
use domain_mod,                   only : domain_t
use exchange_abstract_mod,        only : exchange_t
use mesh_mod,                     only : tile_mesh_t
use sbp_interp_mod,               only : sbp_interp_t
use abstract_interpolators2d_mod, only : interpolator2d_vec2vec_t

implicit none

type, public :: interpolator_w2h_t
    class(exchange_t),               allocatable :: exchange
    class(sbp_interp_t),             allocatable :: sbp_interp_w2v
    class(interpolator2d_vec2vec_t), allocatable :: interp_v2h_op
    type(grid_field_t)                           :: fu, fv, fuh, fvh
contains
    procedure, public :: interp_w2h
end type interpolator_w2h_t

contains

subroutine interp_w2h(this, fh, fw, domain)

    class(interpolator_w2h_t), intent(inout) :: this
    type(grid_field_t),        intent(inout) :: fh, fw
    type(domain_t),            intent(in)    :: domain

    integer(kind=4) :: t

    call this%exchange%do(fw, domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call interp_w2v_tile(this%fu%tile(t), this%fv%tile(t), fw%tile(t), &
                    this%sbp_interp_w2v, domain%mesh_u%tile(t), domain%mesh_v%tile(t))
    end do

    call this%interp_v2h_op%interp2d_vec2vec(this%fuh, this%fvh, this%fu, this%fv, domain)

    call fh%assign(0.5_8, this%fuh, 0.5_8, this%fvh, domain%mesh_p)

end subroutine interp_w2h

subroutine interp_w2v_tile(fu, fv, fw, sbp_interp_w2v, &
                           mesh_u, mesh_v)

    type(tile_field_t),  intent(inout) :: fu, fv
    type(tile_field_t),  intent(inout) :: fw
    class(sbp_interp_t), intent(in)    :: sbp_interp_w2v
    type(tile_mesh_t),   intent(in)    :: mesh_u, mesh_v

    call sbp_interp_w2v%apply_tile(fv, fw, mesh_v, "x")
    call sbp_interp_w2v%apply_tile(fu, fw, mesh_u, "y")

end subroutine interp_w2v_tile
end module interpolator_w2h_mod
