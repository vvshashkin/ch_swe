module vector_advection_Ch_mod

use abstract_v_nabla_mod,          only : v_nabla_operator_t
use abstract_vector_advection_mod, only : vector_advection_operator_t
use grid_field_mod,                only : grid_field_t, tile_field_t
use domain_mod,                    only : domain_t
use halo_mod,                      only : halo_vec_t
use abstract_interpolators2d_mod,  only : interpolator2d_vec2vec_t
use abstract_hor_Christofel_mod,   only : hor_Christofel_t
use parcomm_mod,                   only : parcomm_global
use exchange_abstract_mod,         only : exchange_t
use mesh_mod,                      only : tile_mesh_t

implicit none

type, public, extends(vector_advection_operator_t) :: vector_advection_Ch_t
    class(v_nabla_operator_t), allocatable       :: v_nabla_op
    type(grid_field_t)                           :: u_at_v, v_at_u, uh, vh
    type(grid_field_t)                           :: u_halo, v_halo
    class(interpolator2d_vec2vec_t), allocatable :: interp_v2h_op
    class(interpolator2d_vec2vec_t), allocatable :: interp_h2v_op
    class(hor_Christofel_t),         allocatable :: hor_Christofel
    class(halo_vec_t), allocatable               :: halo_uv!, tendency_edge_sync
    class(exchange_t), allocatable               :: CH_uv_exchange
    
contains
    procedure :: calc_vec_advection
    procedure :: calc_vec_advection_contra
end type vector_advection_Ch_t

contains

subroutine calc_vec_advection(this, u_tend, v_tend, u, v, ut, vt, domain)

    class(vector_advection_Ch_t), intent(inout) :: this
    type(domain_t),               intent(in)    :: domain
    type(grid_field_t),           intent(inout) :: u,  v !covariant components
    type(grid_field_t),           intent(inout) :: ut, vt!contravariant components
    type(grid_field_t),           intent(inout) :: u_tend, v_tend!advective tendencies

    integer(kind=4) :: t

    call this%interp_v2h_op%interp2d_vec2vec(this%uh, this%vh, ut, vt, domain)
    call this%interp_h2v_op%interp2d_vec2vec(this%v_at_u, this%u_at_v, this%vh, this%uh, domain)

    call this%u_halo%assign(u,domain%mesh_u)
    call this%v_halo%assign(v,domain%mesh_v)
    call this%halo_uv%get_halo_vector(this%u_halo, this%v_halo, domain, 3)

    call this%v_nabla_op%calc_v_nabla(u_tend, this%u_halo, ut, this%v_at_u, domain%mesh_u)
    call this%v_nabla_op%calc_v_nabla(v_tend, this%v_halo, this%u_at_v, vt, domain%mesh_v)


    call this%hor_Christofel%add_cov_wind_advection(u_tend, ut, this%v_at_u, domain%mesh_u, 'u')
    call this%hor_Christofel%add_cov_wind_advection(v_tend, this%u_at_v, vt, domain%mesh_v, 'v')

    ! call this%tendency_edge_sync%get_halo_vector(u_tend, v_tend, domain, 1)

    call this%CH_uv_exchange%do_vec(u_tend, v_tend, domain%parcomm)
    do t = domain%partition%ts, domain%partition%te
        call sync_ch_cov(u_tend%tile(t), v_tend%tile(t), domain%mesh_y%tile(t), domain%mesh_x%tile(t))
    end do
    
end subroutine calc_vec_advection

subroutine calc_vec_advection_contra(this, u_tend, v_tend, ut, vt, domain)
    class(vector_advection_Ch_t), intent(inout) :: this
    type(domain_t),               intent(in)    :: domain
    type(grid_field_t),           intent(inout) :: ut, vt!contravariant components
    type(grid_field_t),           intent(inout) :: u_tend, v_tend!advective tendencies

    call parcomm_global%abort("vector_advection_ch_t, calc_vec_advection_contra not implemented")
    
end subroutine calc_vec_advection_contra

subroutine sync_ch_cov(u,v,mesh_u,mesh_v)

    type(tile_field_t), intent(inout) :: u, v
    type(tile_mesh_t),  intent(in)    :: mesh_u, mesh_v

    integer(kind=4) :: i, j, k

    if(mesh_v%is == 1) then
        do k = mesh_v%ks, mesh_v%ke
            do j = mesh_v%js, mesh_v%je
                v%p(1,j,k) = 0.5_8*(v%p(1,j,k)+v%p(0,j,k))
            end do
        end do
    end if
    if(mesh_v%ie == mesh_v%nx) then
        i = mesh_v%nx
        do k = mesh_v%ks, mesh_v%ke
            do j = mesh_v%js, mesh_v%je
                v%p(i,j,k) = 0.5_8*(v%p(i+1,j,k)+v%p(i,j,k))
            end do
        end do
    end if

    if(mesh_u%js == 1) then
        do k = mesh_u%ks, mesh_u%ke
            do i = mesh_u%is, mesh_u%ie
                u%p(i,1,k) = 0.5_8*(u%p(i,1,k)+u%p(i,0,k))
            end do
        end do
    end if
    if(mesh_u%je == mesh_u%ny) then
        j = mesh_u%ny
        do k = mesh_u%ks, mesh_u%ke
            do i = mesh_u%is, mesh_u%ie
                u%p(i,j,k) = 0.5_8*(u%p(i,j+1,k)+u%p(i,j,k))
            end do
        end do
    end if

end subroutine

end module vector_advection_Ch_mod
