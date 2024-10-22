module interpolator_uv2q_sbp_mod

use abstract_interpolators2d_mod, only : interpolator2d_vec2vec_t
use grid_field_mod,               only : grid_field_t, tile_field_t
use domain_mod,                   only : domain_t
use exchange_abstract_mod,        only : exchange_t
use mesh_mod,                     only : tile_mesh_t
use sbp_interp_mod,               only : sbp_interp_t

implicit none

type, public, extends(interpolator2d_vec2vec_t) :: interpolator_uv2q_sbp_t
    logical :: is_Ch = .false.
    class(exchange_t),   allocatable :: exchange
    class(sbp_interp_t), allocatable :: sbp_interp_uv2q
contains
    procedure, public :: interp2d_vec2vec => interp_uv2q_sbp
end type interpolator_uv2q_sbp_t

contains

subroutine interp_uv2q_sbp(this, u, v, u_source, v_source, domain)

    class(interpolator_uv2q_sbp_t), intent(inout) :: this
    type(grid_field_t),                intent(inout) :: u_source, v_source, u, v
    type(domain_t),                    intent(in)    :: domain

    integer(kind=4) :: t

    if(this%is_Ch) then
        call this%exchange%do_vec(v_source, u_source, domain%parcomm)
    else
        call this%exchange%do_vec(u_source, v_source, domain%parcomm)
    end if

    do t = domain%partition%ts, domain%partition%te
        call interp_uv2q_tile(u%tile(t), v%tile(t), u_source%tile(t), v_source%tile(t), &
                              this%sbp_interp_uv2q, domain%mesh_q%tile(t))
    end do

end subroutine interp_uv2q_sbp

subroutine interp_uv2q_tile(u_q, v_q, u, v, sbp_interp_uv2q, mesh_q)

    type(tile_field_t),  intent(inout) :: u_q, v_q
    type(tile_field_t),  intent(inout) :: u, v
    class(sbp_interp_t), intent(in)    :: sbp_interp_uv2q
    type(tile_mesh_t),   intent(in)    :: mesh_q

    call sbp_interp_uv2q%apply_tile(u_q, u, mesh_q, "y")
    call sbp_interp_uv2q%apply_tile(v_q, v, mesh_q, "x")

end subroutine interp_uv2q_tile

end module interpolator_uv2q_sbp_mod
