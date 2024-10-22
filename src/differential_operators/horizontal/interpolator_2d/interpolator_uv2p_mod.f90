module interpolator_uv2p_mod

use abstract_interpolators2d_mod, only : interpolator2d_vec2vec_t
use grid_field_mod,               only : grid_field_t, tile_field_t
use domain_mod,                   only : domain_t
use exchange_abstract_mod,        only : exchange_t
use mesh_mod,                     only : mesh_t, tile_mesh_t
use sbp_interp_mod,               only : sbp_interp_t

implicit none

type, public, extends(interpolator2d_vec2vec_t) :: interpolator2d_uv2p_sbp_C_t
    class(exchange_t),   allocatable :: exchange
    class(sbp_interp_t), allocatable :: sbp_interp_uv2p
    type(mesh_t),        pointer     :: mesh_p
contains
    procedure, public :: interp2d_vec2vec => interp_uv2p_C_sbp
end type interpolator2d_uv2p_sbp_C_t

contains

subroutine interp_uv2p_C_sbp(this, u, v, u_source, v_source, domain)

    class(interpolator2d_uv2p_sbp_C_t), intent(inout) :: this
    type(grid_field_t),                 intent(inout) :: u, v, u_source, v_source
    type(domain_t),                     intent(in)    :: domain

    integer(kind=4) :: t

    call this%exchange%do_vec(u_source, v_source, domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call interp_uv2p_tile(u%tile(t), v%tile(t), u_source%tile(t), v_source%tile(t), &
                               this%sbp_interp_uv2p, this%mesh_p%tile(t))
    end do

end subroutine interp_uv2p_C_sbp

subroutine interp_uv2p_tile(u_p, v_p, u, v, sbp_interp_uv2p, mesh_p)

    type(tile_field_t),  intent(inout) :: u_p, v_p
    type(tile_field_t),  intent(inout) :: u, v
    class(sbp_interp_t), intent(in)    :: sbp_interp_uv2p
    type(tile_mesh_t),   intent(in)    :: mesh_p

    call sbp_interp_uv2p%apply_tile(u_p, u, mesh_p, "x")
    call sbp_interp_uv2p%apply_tile(v_p, v, mesh_p, "y")

end subroutine interp_uv2p_tile

end module interpolator_uv2p_mod
