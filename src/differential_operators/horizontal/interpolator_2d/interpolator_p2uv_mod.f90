module interpolator_p2uv_sbp_C_mod

use abstract_interpolators2d_mod,  only : interpolator2d_scalar2vec_t, &
                                          interpolator2d_vec2vec_t
use grid_field_mod,                only : grid_field_t, tile_field_t
use domain_mod,                    only : domain_t
use exchange_abstract_mod,         only : exchange_t
use mesh_mod,                      only : mesh_t, tile_mesh_t
use sbp_interp_mod,                only : sbp_interp_t

implicit none

type, public, extends(interpolator2d_scalar2vec_t) :: interpolator_p2uv_sbp_C_t
    class(exchange_t),   allocatable :: exchange
    class(sbp_interp_t), allocatable :: sbp_interp_p2uv
contains
    procedure, public :: interp2d_scalar2vec => interp_p2uv_C
end type interpolator_p2uv_sbp_C_t

type, public, extends(interpolator2d_vec2vec_t) :: interpolator_pvec2uv_sbp_C_t
    class(exchange_t),   allocatable :: exchange
    class(sbp_interp_t), allocatable :: sbp_interp_p2uv
    type(mesh_t),        pointer     :: mesh_u, mesh_v
contains
    procedure, public :: interp2d_vec2vec => interp_pvec2uv_C
end type interpolator_pvec2uv_sbp_C_t

type, public, extends(interpolator2d_vec2vec_t) :: interpolator_pvec2uv_sbp_Ch_t
    class(exchange_t),   allocatable :: exchange
    class(sbp_interp_t), allocatable :: sbp_interp_p2uv
contains
    procedure, public :: interp2d_vec2vec => interp_pvec2uv_Ch
end type interpolator_pvec2uv_sbp_Ch_t

contains

subroutine interp_p2uv_C(this, u, v, p, domain)

    class(interpolator_p2uv_sbp_C_t), intent(inout) :: this
    type(grid_field_t),               intent(inout) :: p, u, v
    type(domain_t),                   intent(in)    :: domain

    integer(kind=4) :: t

    call this%exchange%do(p, domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call interp_p2uv_tile(u%tile(t), v%tile(t), p%tile(t), p%tile(t), &
                    this%sbp_interp_p2uv,  domain%mesh_u%tile(t), domain%mesh_v%tile(t))
    end do

end subroutine interp_p2uv_C

subroutine interp_pvec2uv_C(this, u, v, u_source, v_source, domain)

    class(interpolator_pvec2uv_sbp_C_t), intent(inout) :: this
    type(grid_field_t),                  intent(inout) :: u_source, v_source, u, v
    type(domain_t),                      intent(in)    :: domain

    integer(kind=4) :: t

    call this%exchange%do_vec(u_source, v_source, domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call interp_p2uv_tile(u%tile(t), v%tile(t), u_source%tile(t), v_source%tile(t), &
                              this%sbp_interp_p2uv,  this%mesh_u%tile(t), this%mesh_v%tile(t))
    end do

end subroutine interp_pvec2uv_C

subroutine interp_pvec2uv_Ch(this, u, v, u_source, v_source, domain)

    class(interpolator_pvec2uv_sbp_Ch_t), intent(inout) :: this
    type(grid_field_t),                   intent(inout) :: u_source, v_source, u, v
    type(domain_t),                       intent(in)    :: domain

    integer(kind=4) :: t

    call this%exchange%do_vec(u_source, v_source, domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call interp_p2uv_tile(u%tile(t), v%tile(t), u_source%tile(t), v_source%tile(t), &
                              this%sbp_interp_p2uv, domain%mesh_y%tile(t), domain%mesh_x%tile(t))
    end do

end subroutine interp_pvec2uv_Ch

subroutine interp_p2uv_tile(u, v, u_p, v_p, sbp_interp_p2uv, &
                           mesh_u, mesh_v)

    type(tile_field_t),  intent(inout) :: u_p, v_p
    type(tile_field_t),  intent(inout) :: u, v
    class(sbp_interp_t), intent(in)    :: sbp_interp_p2uv
    type(tile_mesh_t),   intent(in)    :: mesh_u, mesh_v

    call sbp_interp_p2uv%apply_tile(u, u_p, mesh_u, "x")
    call sbp_interp_p2uv%apply_tile(v, v_p, mesh_v, "y")

end subroutine interp_p2uv_tile

end module interpolator_p2uv_sbp_C_mod
