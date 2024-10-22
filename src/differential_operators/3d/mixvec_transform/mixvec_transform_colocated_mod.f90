module mixvec_transform_colocated_mod

use abstract_mixvec_transform_mod, only : mixvec_transform_t
use domain_mod,                    only : domain_t
use mesh_mod,                      only : tile_mesh_t
use grid_field_mod,                only : grid_field_t, tile_field_t

implicit none

type, extends(mixvec_transform_t) :: mixvec_transform_colocated_t
contains
    procedure :: transform_co2mix
    procedure :: transform_mix2contra
end type

contains

subroutine transform_co2mix(this, u_contra, v_contra, w, u_cov, v_cov, w_cov, domain)
    class(mixvec_transform_colocated_t), intent(inout) :: this
    type(domain_t),                      intent(in)    :: domain
    type(grid_field_t),                  intent(inout) :: u_cov, v_cov, w_cov
    !output:
    type(grid_field_t),                  intent(inout) :: u_contra, v_contra, w

    integer(kind=4) :: t

    do t = domain%mesh_w%ts, domain%mesh_w%te
        call transform_co2mix_tile(u_contra%tile(t), v_contra%tile(t), w%tile(t), &
                                   u_cov%tile(t), v_cov%tile(t), w_cov%tile(t),   &
                                   domain%mesh_w%tile(t))
    end do
end subroutine transform_co2mix

subroutine transform_co2mix_tile(u_contra, v_contra, w, u_cov, v_cov, w_cov, mesh)
    type(tile_field_t), intent(inout) :: u_contra, v_contra, w
    type(tile_field_t), intent(in)    :: u_cov, v_cov, w_cov
    type(tile_mesh_t),  intent(in)    :: mesh

    integer(kind=4) :: i, j, k
    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                u_contra%p(i,j,k) = mesh%Qi(1,i,j,k)*u_cov%p(i,j,k)+&
                                    mesh%Qi(2,i,j,k)*v_cov%p(i,j,k)+&
                                    mesh%Qi(4,i,j,k)*w_cov%p(i,j,k)
                v_contra%p(i,j,k) = mesh%Qi(2,i,j,k)*u_cov%p(i,j,k)+&
                                    mesh%Qi(3,i,j,k)*v_cov%p(i,j,k)+&
                                    mesh%Qi(5,i,j,k)*w_cov%p(i,j,k)
                w%p(i,j,k) = w_cov%p(i,j,k) / mesh%a3(4,i,j,k)
            end do
        end do
    end do

end subroutine transform_co2mix_tile

subroutine transform_mix2contra(this, w_contra, u_contra, v_contra, w, domain)
    class(mixvec_transform_colocated_t), intent(inout) :: this
    type(grid_field_t),                  intent(inout) :: u_contra, v_contra, w
    type(domain_t),                      intent(in)    :: domain
    !output:
    type(grid_field_t),                  intent(inout) :: w_contra

    integer(kind=4) :: t

    do t = domain%mesh_w%ts, domain%mesh_w%te
        call transform_mix2contra_tile(w_contra%tile(t), u_contra%tile(t), &
                                       v_contra%tile(t), w%tile(t), domain%mesh_w%tile(t))
    end do
end subroutine transform_mix2contra

subroutine transform_mix2contra_tile(w_contra, u_contra, v_contra, w, mesh)
    type(tile_field_t), intent(inout) :: w_contra
    type(tile_field_t), intent(in)    :: u_contra, v_contra, w
    type(tile_mesh_t),  intent(in)    :: mesh

    integer(kind=4) :: i, j, k

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                w_contra%p(i,j,k) = (w%p(i,j,k)-u_contra%p(i,j,k)*mesh%a1(4,i,j,k)-&
                                                v_contra%p(i,j,k)*mesh%a2(4,i,j,k))&
                                                / mesh%a3(4,i,j,k)
            end do
        end do
    end do
end subroutine transform_mix2contra_tile

end module mixvec_transform_colocated_mod
