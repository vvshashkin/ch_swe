module interpolator_q2uv_sbp_mod

use grid_field_mod,               only : grid_field_t, tile_field_t
use domain_mod,                   only : domain_t
use exchange_abstract_mod,        only : exchange_t
use mesh_mod,                     only : tile_mesh_t
use sbp_interp_mod,               only : sbp_interp_t
use abstract_interpolators2d_mod, only : interpolator2d_vec2vec_t

implicit none

type, public, extends(interpolator2d_vec2vec_t) :: interpolator_q2uv_sbp_Ch_t
    class(exchange_t),   allocatable :: exchange
    class(sbp_interp_t), allocatable :: sbp_interp_q2uv
contains
    procedure, public :: interp2d_vec2vec => interp_q2uv_Ch_sbp
end type interpolator_q2uv_sbp_Ch_t

contains

subroutine interp_q2uv_Ch_sbp(this, u, v, u_source, v_source, domain)

    class(interpolator_q2uv_sbp_Ch_t), intent(inout) :: this
    type(grid_field_t),                intent(inout) :: u_source, v_source, u, v
    type(domain_t),                    intent(in)    :: domain

    integer(kind=4) :: t

    call this%exchange%do_vec(u_source, v_source, domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call interp_q2uv_tile(u%tile(t), v%tile(t), u_source%tile(t), v_source%tile(t), &
                    this%sbp_interp_q2uv, domain%mesh_y%tile(t), &
                    domain%mesh_x%tile(t))
    end do

end subroutine interp_q2uv_Ch_sbp

subroutine interp_q2uv_tile(u, v, u_w, v_w, sbp_interp_q2uv, &
                           mesh_u, mesh_v)

    type(tile_field_t),  intent(inout) :: u_w, v_w
    type(tile_field_t),  intent(inout) :: u, v
    class(sbp_interp_t), intent(in)    :: sbp_interp_q2uv
    type(tile_mesh_t),   intent(in)    :: mesh_u, mesh_v

    call sbp_interp_q2uv%apply_tile(u, u_w, mesh_u, "y")
    call sbp_interp_q2uv%apply_tile(v, v_w, mesh_v, "x")

end subroutine interp_q2uv_tile

end module interpolator_q2uv_sbp_mod
