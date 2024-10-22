module mixvec_transform_staggered_mod

use abstract_mixvec_transform_mod,  only : mixvec_transform_t
use domain_mod,                     only : domain_t
use mesh_mod,                       only : tile_mesh_t
use grid_field_mod,                 only : grid_field_t, tile_field_t
use abstract_vertical_operator_mod, only : vertical_operator_t
use abstract_interpolators2d_mod,   only : interpolator2d_vec2vec_t


implicit none

type, extends(mixvec_transform_t) :: mixvec_transform_staggered_t
    class(vertical_operator_t),      allocatable :: w2p_interp, p2w_interp
    class(interpolator2d_vec2vec_t), allocatable :: uv2p_interp, p2uv_interp
    type(grid_field_t)                           :: up, vp, wp
contains
    procedure :: transform_co2mix
    procedure :: transform_mix2contra
end type mixvec_transform_staggered_t

contains

subroutine transform_co2mix(this, u_contra, v_contra, w, u_cov, v_cov, w_cov, domain)
    class(mixvec_transform_staggered_t), intent(inout) :: this
    type(domain_t),                      intent(in)    :: domain
    type(grid_field_t),                  intent(inout) :: u_cov, v_cov, w_cov
    !output:
    type(grid_field_t),                  intent(inout) :: u_contra, v_contra, w

    integer(kind=4) :: t

    call this%w2p_interp%apply(this%wp,w_cov,domain)
    call this%uv2p_interp%interp2d_vec2vec(this%up,this%vp,u_cov,v_cov,domain)
    do t= domain%mesh_p%ts, domain%mesh_p%te
        call calculate_uv_cross_metric_terms(this%up%tile(t), this%vp%tile(t), &
                                          this%wp%tile(t), domain%mesh_p%tile(t))
    end do
    call this%p2uv_interp%interp2d_vec2vec(u_contra, v_contra, this%up, this%vp,domain)
    do t = domain%mesh_p%ts, domain%mesh_p%te
        call get_uv_contra_tile(u_contra%tile(t),v_contra%tile(t), &
                                u_cov%tile(t),v_cov%tile(t),       &
                                domain%mesh_u%tile(t), domain%mesh_v%tile(t))
        call get_real_w_tile(w%tile(t), w_cov%tile(t), domain%mesh_w%tile(t))
    end do

end subroutine transform_co2mix

subroutine calculate_uv_cross_metric_terms(u,v,w,mesh)
    type(tile_field_t), intent(inout) :: u, v
    type(tile_field_t), intent(in)    :: w
    type(tile_mesh_t),  intent(in)    :: mesh

    integer(kind=4) :: i, j, k
    real(kind=8)    :: u0, v0, w0

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                u0 = u%p(i,j,k); v0 = v%p(i,j,k); w0 = w%p(i,j,k)
                u%p(i,j,k) = (mesh%Qi(2,i,j,k)*v0+mesh%Qi(4,i,j,k)*w0)*mesh%J(i,j,k)
                v%p(i,j,k) = (mesh%Qi(2,i,j,k)*u0+mesh%Qi(5,i,j,k)*w0)*mesh%J(i,j,k)
            end do
        end do
    end do
end subroutine calculate_uv_cross_metric_terms

subroutine get_real_w_tile(w, w_cov, mesh)
    type(tile_field_t), intent(inout) :: w
    type(tile_field_t), intent(in)    :: w_cov
    type(tile_mesh_t),  intent(in)    :: mesh

    integer(kind=4) :: i, j, k

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                w%p(i,j,k) = w_cov%p(i,j,k) / mesh%a3(4,i,j,k)
            end do
        end do
    end do
end subroutine get_real_w_tile

subroutine get_uv_contra_tile(u_contra, v_contra, u_cov, v_cov, mesh_u, mesh_v)
    type(tile_field_t), intent(inout) :: u_contra, v_contra
    type(tile_field_t), intent(in)    :: u_cov, v_cov
    type(tile_mesh_t),  intent(in)    :: mesh_u, mesh_v

    integer(kind=4) :: i, j, k

    do k = mesh_u%ks, mesh_u%ke
        do j = mesh_u%js, mesh_u%je
            do i = mesh_u%is, mesh_u%ie
                u_contra%p(i,j,k) = mesh_u%Qi(1,i,j,k)*u_cov%p(i,j,k)+&
                                    u_contra%p(i,j,k) / mesh_u%J(i,j,k)
            end do
        end do
    end do
    do k = mesh_v%ks, mesh_v%ke
        do j = mesh_v%js, mesh_v%je
            do i = mesh_v%is, mesh_v%ie
                v_contra%p(i,j,k) = mesh_v%Qi(3,i,j,k)*v_cov%p(i,j,k)+&
                                    v_contra%p(i,j,k) / mesh_v%J(i,j,k)
            end do
        end do
    end do

end subroutine get_uv_contra_tile

subroutine transform_mix2contra(this, w_contra, u_contra, v_contra, w, domain)
    class(mixvec_transform_staggered_t), intent(inout) :: this
    type(grid_field_t),                  intent(inout) :: u_contra, v_contra, w
    type(domain_t),                      intent(in)    :: domain
    !output:
    type(grid_field_t),                  intent(inout) :: w_contra

    integer(kind=4) :: t

    call this%uv2p_interp%interp2d_vec2vec(this%up,this%vp,u_contra,v_contra,domain)
    do t = domain%mesh_p%ts, domain%mesh_p%te
        call calculate_Quv_tile(this%wp%tile(t), this%up%tile(t), &
                                this%vp%tile(t), domain%mesh_p%tile(t))
    end do

    call this%p2w_interp%apply(w_contra,this%wp,domain)
    do t = domain%mesh_w%ts, domain%mesh_w%te
        call calculate_w_contra_tile(w_contra%tile(t), w%tile(t), domain%mesh_w%tile(t))
    end do

end subroutine transform_mix2contra

subroutine calculate_w_contra_tile(w_contra, w, mesh)
    type(tile_field_t), intent(inout) :: w_contra
    type(tile_field_t), intent(in)    :: w
    type(tile_mesh_t),  intent(in)    :: mesh

    integer(kind=4) :: i, j, k

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                w_contra%p(i,j,k) = w%p(i,j,k) / mesh%a3(4,i,j,k) + &
                                    w_contra%p(i,j,k) / mesh%J(i,j,k)
            end do
        end do
    end do
end subroutine calculate_w_contra_tile

subroutine calculate_Quv_tile(Quv, u_contra, v_contra, mesh)
    type(tile_field_t), intent(inout) :: Quv
    type(tile_field_t), intent(in)    :: u_contra, v_contra
    type(tile_mesh_t),   intent(in)    :: mesh

    integer(kind=4) :: i, j, k

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                Quv%p(i,j,k) = -(u_contra%p(i,j,k)*mesh%a1(4,i,j,k)+   &
                                 v_contra%p(i,j,k)*mesh%a2(4,i,j,k)) * &
                                 mesh%J(i,j,k) / mesh%a3(4,i,j,k)
            end do
        end do
    end do
end subroutine calculate_Quv_tile

end module mixvec_transform_staggered_mod
