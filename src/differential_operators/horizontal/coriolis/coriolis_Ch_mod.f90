module coriolis_Ch_mod

use grid_field_mod,                only : grid_field_t, tile_field_t
use domain_mod,                    only : domain_t
use abstract_coriolis_mod,         only : coriolis_operator_t
use mesh_mod,                      only : tile_mesh_t
use interpolator_w2h_mod,          only : interpolator_w2h_t
use abstract_interpolators2d_mod,  only : interpolator2d_vec2vec_t
use vec_math_mod,                  only : divide_by_J_self, multiply_by_J_self, multiply_by_J
use halo_mod,                      only : halo_vec_t
use exchange_abstract_mod,         only : exchange_t

implicit none

type, public, extends(coriolis_operator_t) :: coriolis_Ch_t

    integer(kind=4) :: option = 1

    type(grid_field_t) :: fJ2  !coriolis parameter multiplied by J**2
    type(grid_field_t) :: up, vp, curl_u, curl_v, curl_pu, curl_pv

    class(interpolator2d_vec2vec_t), allocatable :: interp_v2h_op
    class(interpolator2d_vec2vec_t), allocatable :: interp_h2v_op
    class(interpolator2d_vec2vec_t), allocatable :: interp_q2v_op

    class(halo_vec_t), allocatable :: uv_sync

    class(exchange_t), allocatable :: CH_uv_exchange

contains
    procedure, public :: calc_coriolis
    procedure, public :: calc_coriolis_vec_inv
end type coriolis_Ch_t

type, extends(coriolis_Ch_t) :: coriolis_CH_opt2_t
contains
    procedure, public :: calc_coriolis => calc_coriolis_opt2
end type coriolis_CH_opt2_t

contains

subroutine calc_coriolis(this, cor_u, cor_v, ut, vt, domain)
    class(coriolis_Ch_t),    intent(inout) :: this
    type(domain_t),          intent(in)    :: domain
    type(grid_field_t),      intent(inout) :: ut, vt!contravariant components
    type(grid_field_t),      intent(inout) :: cor_u, cor_v

    integer(kind=4) :: t

    call this%interp_v2h_op%interp2d_vec2vec(this%up, this%vp, ut, vt, domain)

    do t = domain%partition%ts, domain%partition%te
        call calc_cori_at_p_tile(this%up%tile(t), this%vp%tile(t), &
                                 this%fJ2%tile(t), domain%mesh_xy%tile(t))
    end do

    call this%uv_sync%get_halo_vector(this%up, this%vp, domain, 0)

    call this%interp_h2v_op%interp2d_vec2vec(cor_u, cor_v, this%up, this%vp, domain)

    call divide_by_J_self(cor_u, domain%mesh_y)
    call divide_by_J_self(cor_v, domain%mesh_x)

    ! call this%CH_uv_exchange%do_vec(cor_u, cor_v, domain%parcomm)
    ! do t = domain%partition%ts, domain%partition%te
    !     call sync_ch_cov(cor_u%tile(t), cor_v%tile(t), domain%mesh_y%tile(t), domain%mesh_x%tile(t))
    ! end do

end subroutine calc_coriolis

subroutine calc_coriolis_opt2(this, cor_u, cor_v, ut, vt, domain)
    class(coriolis_Ch_opt2_t), intent(inout) :: this
    type(domain_t),            intent(in)    :: domain
    type(grid_field_t),        intent(inout) :: ut, vt!contravariant components
    type(grid_field_t),        intent(inout) :: cor_u, cor_v

    integer(kind=4) :: t


    call multiply_by_J(cor_u, ut, domain%mesh_y)
    call multiply_by_J(cor_v, vt, domain%mesh_x)
    call this%interp_v2h_op%interp2d_vec2vec(this%up, this%vp, cor_u, cor_v, domain)

    call this%uv_sync%get_halo_vector(this%up, this%vp, domain, 0)

    do t = domain%partition%ts, domain%partition%te
        call calc_cori_at_p_tile(this%up%tile(t), this%vp%tile(t), &
                                 this%fJ2%tile(t), domain%mesh_xy%tile(t))
    end do

    ! call this%uv_sync%get_halo_vector(this%up, this%vp, domain, 0)

    call this%interp_h2v_op%interp2d_vec2vec(cor_u, cor_v, this%up, this%vp, domain)

end subroutine calc_coriolis_opt2

subroutine calc_cori_at_p_tile(u, v, f, mesh)

    type(tile_field_t), intent(inout) :: u, v
    type(tile_field_t), intent(in)    :: f
    type(tile_mesh_t),  intent(in)    :: mesh

    integer(kind=4) :: i, j, k
    real(kind=8)    :: tmp_u, tmp_v

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie

                tmp_u = f%p(i,j,1)*v%p(i,j,k)
                tmp_v =-f%p(i,j,1)*u%p(i,j,k)
                
                u%p(i,j,k) = tmp_u
                v%p(i,j,k) = tmp_v
            
            end do
        end do
    end do

end subroutine calc_cori_at_p_tile

subroutine calc_coriolis_vec_inv(this, cor_u, cor_v, hu, hv, h, curl, domain)
    class(coriolis_Ch_t),    intent(inout) :: this
    type(domain_t),          intent(in)    :: domain
    type(grid_field_t),      intent(inout) :: hu, hv! massflux contravariant components
    type(grid_field_t),      intent(inout) :: h, curl
    type(grid_field_t),      intent(inout) :: cor_u, cor_v

    integer(kind=4) :: t

    call this%interp_q2v_op%interp2d_vec2vec(this%curl_u, this%curl_v, curl, curl, domain)
    call this%interp_v2h_op%interp2d_vec2vec(this%curl_pu, this%curl_pv, this%curl_u, this%curl_v, domain)
    call this%curl_pu%assign(0.5_8, this%curl_pu, 0.5_8, this%curl_pv, domain%mesh_xy)
    call multiply_by_J_self(this%curl_pu,domain%mesh_xy)
    call multiply_by_J_self(this%curl_pu,domain%mesh_xy)
    call this%curl_pu%update(1.0_8,this%fJ2,domain%mesh_xy)

    call this%interp_v2h_op%interp2d_vec2vec(this%up, this%vp, hu, hv, domain)

    do t = domain%partition%ts, domain%partition%te
        call calc_cori_at_p_tile(this%up%tile(t), this%vp%tile(t), &
                                 this%curl_pu%tile(t), domain%mesh_xy%tile(t))
    end do

    call this%uv_sync%get_halo_vector(this%up, this%vp, domain, 0)

    call this%interp_h2v_op%interp2d_vec2vec(cor_u, cor_v, this%up, this%vp, domain)
    call this%interp_h2v_op%interp2d_vec2vec(this%up, this%vp, h, h, domain)

    call divide_by_J_self(cor_u, domain%mesh_y)
    call divide_by_J_self(cor_v, domain%mesh_x)

    call cor_u%assign_ratio(1.0_8, cor_u, this%up, domain%mesh_y)
    call cor_v%assign_ratio(1.0_8, cor_v, this%vp, domain%mesh_x)

    ! call this%CH_uv_exchange%do_vec(cor_u, cor_v, domain%parcomm)
    ! do t = domain%partition%ts, domain%partition%te
    !     call sync_ch_cov(cor_u%tile(t), cor_v%tile(t), domain%mesh_y%tile(t), domain%mesh_x%tile(t))
    ! end do

end subroutine

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

end module