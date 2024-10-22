module vector_advection_Ah_mod

use grid_field_mod,                only : grid_field_t, tile_field_t
use domain_mod,                    only : domain_t
use halo_mod,                      only : halo_vec_t
use exchange_abstract_mod,         only : exchange_t
use abstract_vector_advection_mod, only : vector_advection_operator_t
use abstract_v_nabla_mod,          only : v_nabla_operator_t
use abstract_hor_Christofel_mod,   only : hor_Christofel_t

implicit none

type, public, extends(vector_advection_operator_t) :: vector_advection_Ah_t
    class(halo_vec_t),         allocatable :: sync_edges_cov
    class(halo_vec_t),         allocatable :: sync_edges_contra
    class(exchange_t),         allocatable :: exch_uv_interior
    class(v_nabla_operator_t), allocatable :: v_nabla_op
    class(hor_Christofel_t),   allocatable :: hor_Christofel
contains
    procedure :: calc_vec_advection
    procedure :: calc_vec_advection_contra
end type vector_advection_Ah_t

contains

subroutine calc_vec_advection(this, u_tend, v_tend, u, v, ut, vt, domain)
    class(vector_advection_Ah_t), intent(inout) :: this
    type(domain_t),               intent(in)    :: domain
    type(grid_field_t),           intent(inout) :: u,  v!covariant components
    type(grid_field_t),           intent(inout) :: ut, vt!contravariant components
    type(grid_field_t),           intent(inout) :: u_tend, v_tend!advective tendencies

    integer(kind=4) :: t

    call this%exch_uv_interior%do_vec(u,v,domain%parcomm)

    call this%v_nabla_op%calc_v_nabla(u_tend, u, ut, vt, domain%mesh_u)
    call this%v_nabla_op%calc_v_nabla(v_tend, v, ut, vt, domain%mesh_v)

    ! call this%hor_Christofel%add(u_tend, u, v, ut, vt, domain%mesh_u, 'u')
    ! call this%hor_Christofel%add(v_tend, u, v, ut, vt, domain%mesh_v, 'v')
    call this%hor_Christofel%add_cov_wind_advection(u_tend, ut, vt, domain%mesh_u, 'u')
    call this%hor_Christofel%add_cov_wind_advection(v_tend, ut, vt, domain%mesh_v, 'v')

    call this%sync_edges_cov%get_halo_vector(u_tend,v_tend,domain,0)
end subroutine calc_vec_advection

subroutine calc_vec_advection_contra(this, u_tend, v_tend, ut, vt, domain)
    class(vector_advection_Ah_t), intent(inout) :: this
    type(domain_t),               intent(in)    :: domain
    type(grid_field_t),           intent(inout) :: ut, vt!contravariant components
    type(grid_field_t),           intent(inout) :: u_tend, v_tend!advective tendencies

    integer(kind=4) :: t

    call this%exch_uv_interior%do_vec(ut,vt,domain%parcomm)

    call this%v_nabla_op%calc_v_nabla(u_tend, ut, ut, vt, domain%mesh_u)
    call this%v_nabla_op%calc_v_nabla(v_tend, vt, ut, vt, domain%mesh_v)

    call this%hor_Christofel%add(u_tend, ut, vt, domain%mesh_u, 'u')
    call this%hor_Christofel%add(v_tend, ut, vt, domain%mesh_v, 'v')

    call this%sync_edges_contra%get_halo_vector(u_tend,v_tend,domain,0)
end subroutine calc_vec_advection_contra

end module vector_advection_Ah_mod
