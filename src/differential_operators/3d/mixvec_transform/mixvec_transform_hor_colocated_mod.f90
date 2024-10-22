module mixvec_transform_hor_colocated_mod

use abstract_mixvec_transform_mod,  only : mixvec_transform_t
use domain_mod,                     only : domain_t
use mesh_mod,                       only : tile_mesh_t
use grid_field_mod,                 only : grid_field_t, tile_field_t
use abstract_vertical_operator_mod, only : vertical_operator_t

implicit none

type, extends(mixvec_transform_t) :: mixvec_transform_hor_colocated_t
    class(vertical_operator_t), allocatable :: w2p_interp, p2w_interp
    type(grid_field_t)                      :: wp
contains
    procedure :: transform_co2mix
    procedure :: transform_mix2contra
end type

contains

subroutine transform_co2mix(this, u_contra, v_contra, w, u_cov, v_cov, w_cov, domain)
    class(mixvec_transform_hor_colocated_t), intent(inout) :: this
    type(domain_t),                          intent(in)    :: domain
    type(grid_field_t),                      intent(inout) :: u_cov, v_cov, w_cov
    !output:
    type(grid_field_t),                      intent(inout) :: u_contra, v_contra, w

    integer(kind=4) :: t

    call this%w2p_interp%apply(this%wp,w_cov,domain)

    do t = domain%mesh_p%ts, domain%mesh_p%te
        call get_uv_contra_tile(u_contra%tile(t), v_contra%tile(t),             &
                                u_cov%tile(t), v_cov%tile(t), this%wp%tile(t),  &
                                domain%mesh_p%tile(t))
        call get_real_w_tile(w%tile(t), w_cov%tile(t), domain%mesh_w%tile(t))
    end do
end subroutine transform_co2mix

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

subroutine get_uv_contra_tile(u_contra, v_contra, u_cov, v_cov, wp_cov, mesh)
    type(tile_field_t), intent(inout) :: u_contra, v_contra
    type(tile_field_t), intent(in)    :: u_cov, v_cov, wp_cov
    type(tile_mesh_t),  intent(in)    :: mesh

    integer(kind=4) :: i, j, k
    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                u_contra%p(i,j,k) = mesh%Qi(1,i,j,k)*u_cov%p(i,j,k)+&
                                    mesh%Qi(2,i,j,k)*v_cov%p(i,j,k)+&
                                    mesh%Qi(4,i,j,k)*wp_cov%p(i,j,k)
                v_contra%p(i,j,k) = mesh%Qi(2,i,j,k)*u_cov%p(i,j,k)+&
                                    mesh%Qi(3,i,j,k)*v_cov%p(i,j,k)+&
                                    mesh%Qi(5,i,j,k)*wp_cov%p(i,j,k)
            end do
        end do
    end do

end subroutine get_uv_contra_tile

subroutine transform_mix2contra(this, w_contra, u_contra, v_contra, w, domain)
    class(mixvec_transform_hor_colocated_t), intent(inout) :: this
    type(grid_field_t),                      intent(inout) :: u_contra, v_contra, w
    type(domain_t),                          intent(in)    :: domain
    !output:
    type(grid_field_t),                      intent(inout) :: w_contra

    integer(kind=4) :: t

    ! call w_contra%assign(0.0_8,domain%mesh_w)
    do t = domain%mesh_p%ts, domain%mesh_p%te
        call calculate_Quv_tile(this%wp%tile(t), u_contra%tile(t), &
                                       v_contra%tile(t), domain%mesh_p%tile(t))
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

end module mixvec_transform_hor_colocated_mod
