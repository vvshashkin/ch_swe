module shallow_atm_coriolis_C_mod

use abstract_coriolis_3d_mod,     only : coriolis3d_operator_t
use grid_field_mod,               only : grid_field_t, tile_field_t
use mesh_mod,                     only : tile_mesh_t
use domain_mod,                   only : domain_t
use abstract_interpolators2d_mod, only : interpolator2d_vec2vec_t

implicit none

type, public, extends(coriolis3d_operator_t) :: shallow_atm_coriolis_C_t
    class(interpolator2d_vec2vec_t), allocatable :: interp2d_p2uv, interp2d_uv2p
    type(grid_field_t)                           :: Q2d_uu, Q2d_vv, Q2d_uv
    type(grid_field_t)                           :: Jp2d, Ju2d, Jv2d
    type(grid_field_t)                           :: fcori_p, up, vp, u_cov, v_cov
contains
    procedure :: calc_coriolis
end type shallow_atm_coriolis_C_t

contains

subroutine calc_coriolis(this, u_tend, v_tend, w_tend, u, v, w, domain)
    class(shallow_atm_coriolis_C_t), intent(inout) :: this
    type(domain_t),                  intent(in)    :: domain
    type(grid_field_t),              intent(inout) :: u, v, w
    type(grid_field_t),              intent(inout) :: u_tend, v_tend, w_tend

    integer(kind=4) :: t

    call w_tend%assign(0.0_8,domain%mesh_w)

    call transform_contra2co_2d(this%u_cov, this%v_cov, u, v, this%up,this%vp, &
                                this%Q2d_uu,this%Q2d_uv,this%Q2d_vv,           &
                                this%Ju2d, this%Jv2d, this%Jp2d,               &
                                this%interp2d_p2uv, this%interp2d_uv2p, domain)

    call this%interp2d_uv2p%interp2d_vec2vec(this%up,this%vp,this%u_cov, &
                                                             this%v_cov,domain)
    do t = domain%mesh_p%ts, domain%mesh_p%te
        call calc_coriolis_tend_p_tile(this%up%tile(t),this%vp%tile(t),  &
                                       this%fcori_p%tile(t), domain%mesh_p%tile(t))
    end do
    call this%interp2d_p2uv%interp2d_vec2vec(u_tend,v_tend,this%up,this%vp,domain)

    call u_tend%divide_by_2dfield(this%Ju2d,domain%mesh_u)
    call v_tend%divide_by_2dfield(this%Jv2d,domain%mesh_v)
end subroutine calc_coriolis

subroutine transform_contra2co_2d(u_cov,v_cov,u,v,up,vp,Quu,Quv,Qvv,Ju,Jv,Jp,&
                                  interp2d_p2uv, interp2d_uv2p, domain)
    type(grid_field_t),              intent(inout) :: u_cov, v_cov
    type(grid_field_t),              intent(inout) :: u, v
    type(grid_field_t),              intent(inout) :: up, vp
    type(grid_field_t),              intent(in)    :: Quu, Quv, Qvv, Ju, Jv, Jp
    class(interpolator2d_vec2vec_t), intent(inout) :: interp2d_p2uv, interp2d_uv2p
    type(domain_t),                  intent(in)    :: domain

    integer(kind=4) :: t

    call interp2d_uv2p%interp2d_vec2vec(up,vp,u,v,domain)

    do t = domain%mesh_p%ts, domain%mesh_p%te
        call calc_contra2co_p_part_tile(up%tile(t),vp%tile(t),Quv%tile(t), &
                                        Jp%tile(t), domain%mesh_p%tile(t))
    end do

    call interp2d_p2uv%interp2d_vec2vec(u_cov, v_cov, up, vp, domain)

    do t = domain%mesh_u%ts, domain%mesh_u%te
        call finalize_uv_cov(u_cov%tile(t),u%tile(t),Quu%tile(t), &
                             Ju%tile(t),domain%mesh_u%tile(t))
        call finalize_uv_cov(v_cov%tile(t),v%tile(t),Qvv%tile(t), &
                             Jv%tile(t),domain%mesh_v%tile(t))
    end do

end subroutine

subroutine calc_contra2co_p_part_tile(u,v,Q_uv,Jp,mesh)
    type(tile_field_t), intent(inout) :: u, v
    type(tile_field_t), intent(in)    :: Q_uv, Jp
    type(tile_mesh_t),  intent(in)    :: mesh

    integer(kind=4) :: i, j, k
    real(kind=8)    :: up, vp

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                up = u%p(i,j,k); vp = v%p(i,j,k)
                u%p(i,j,k) = Q_uv%p(i,j,1)*Jp%p(i,j,1)*vp
                v%p(i,j,k) = Q_uv%p(i,j,1)*Jp%p(i,j,1)*up
            end do
        end do
    end do
end subroutine

subroutine finalize_uv_cov(u_cov,u,Q,Ju,mesh)
    type(tile_field_t), intent(inout) :: u_cov
    type(tile_field_t), intent(in)    :: u, Q, Ju
    type(tile_mesh_t),  intent(in)    :: mesh

    integer(kind=4) :: i, j, k

    do k=mesh%ks,mesh%ke; do j=mesh%js,mesh%je; do i=mesh%is,mesh%ie
        u_cov%p(i,j,k) = Q%p(i,j,1)*u%p(i,j,k)+u_cov%p(i,j,k)/Ju%p(i,j,1)
    end do; end do; end do
end subroutine

subroutine calc_coriolis_tend_p_tile(u,v,fcori,mesh)
    type(tile_field_t), intent(inout) :: u, v
    type(tile_field_t), intent(in)    :: fcori
    type(tile_mesh_t),  intent(in)    :: mesh

    integer(kind=4) :: i, j, k
    real(kind=8)    :: u_tend, v_tend

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                u_tend = fcori%p(i,j,1)*v%p(i,j,k)
                v_tend =-fcori%p(i,j,1)*u%p(i,j,k)
                u%p(i,j,k) = u_tend; v%p(i,j,k) = v_tend
            end do
        end do
    end do

end subroutine

end module shallow_atm_coriolis_C_mod
