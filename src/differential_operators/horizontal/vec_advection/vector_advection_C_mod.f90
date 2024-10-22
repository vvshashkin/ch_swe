module vector_advection_C_mod

use abstract_v_nabla_mod,          only : v_nabla_operator_t
use abstract_vector_advection_mod, only : vector_advection_operator_t
use grid_field_mod,                only : grid_field_t, tile_field_t
use domain_mod,                    only : domain_t
use halo_mod,                      only : halo_vec_t
use abstract_interpolators2d_mod,  only : interpolator2d_vec2vec_t
use abstract_hor_Christofel_mod,   only : hor_Christofel_t
use parcomm_mod,                   only : parcomm_global

implicit none

type, public, extends(vector_advection_operator_t) :: vector_advection_C_t
    class(v_nabla_operator_t), allocatable       :: v_nabla_op
    type(grid_field_t)                           :: u_at_v, v_at_u, uh, vh
    type(grid_field_t)                           :: u_halo, v_halo
    class(interpolator2d_vec2vec_t), allocatable :: interp_v2h_op
    class(interpolator2d_vec2vec_t), allocatable :: interp_h2v_op
    class(hor_Christofel_t),         allocatable :: hor_Christofel
    class(halo_vec_t), allocatable               :: halo_uv, tendency_edge_sync
contains
    procedure :: calc_vec_advection
    procedure :: calc_vec_advection_contra
end type vector_advection_C_t

contains

subroutine calc_vec_advection(this, u_tend, v_tend, u, v, ut, vt, domain)
    class(vector_advection_C_t),  intent(inout) :: this
    type(domain_t),               intent(in)    :: domain
    type(grid_field_t),           intent(inout) :: u,  v!covariant components
    type(grid_field_t),           intent(inout) :: ut, vt!contravariant components
    type(grid_field_t),           intent(inout) :: u_tend, v_tend!advective tendencies

    call parcomm_global%abort("vector_advection_c_t, calc_vec_advection not implemented")

end subroutine calc_vec_advection

subroutine calc_vec_advection_contra(this, u_tend, v_tend, ut, vt, domain)
    class(vector_advection_C_t),  intent(inout) :: this
    type(domain_t),               intent(in)    :: domain
    type(grid_field_t),           intent(inout) :: ut, vt!contravariant components
    type(grid_field_t),           intent(inout) :: u_tend, v_tend!advective tendencies

    integer(kind=4) :: t

    call this%interp_v2h_op%interp2d_vec2vec(this%uh, this%vh, ut, vt, domain)
    call this%interp_h2v_op%interp2d_vec2vec(this%v_at_u, this%u_at_v, this%vh, this%uh, domain)

    call this%u_halo%assign(ut,domain%mesh_u)
    call this%v_halo%assign(vt,domain%mesh_v)
    call this%halo_uv%get_halo_vector(this%u_halo, this%v_halo, domain, 3)

    call this%v_nabla_op%calc_v_nabla(u_tend, this%u_halo, this%u_halo, this%v_at_u, domain%mesh_u)
    call this%v_nabla_op%calc_v_nabla(v_tend, this%v_halo, this%u_at_v, this%v_halo, domain%mesh_v)


    call this%hor_Christofel%add(u_tend, this%u_halo, this%v_at_u, domain%mesh_u, 'u')
    call this%hor_Christofel%add(v_tend, this%u_at_v, this%v_halo, domain%mesh_v, 'v')

    call this%tendency_edge_sync%get_halo_vector(u_tend, v_tend, domain, 1)
end subroutine calc_vec_advection_contra

end module vector_advection_C_mod
