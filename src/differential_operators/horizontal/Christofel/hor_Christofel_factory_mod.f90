module hor_Christofel_factory_mod

use abstract_hor_Christofel_mod, only : hor_Christofel_t
use hor_Christofel_shallow_mod,  only : hor_Christofel_shallow_t
use mesh_mod,                    only : mesh_t
use grid_field_mod,              only : grid_field_t
use grid_field_factory_mod,      only : create_grid_field, create_grid_field_2d
use metric_mod,                  only : metric_t
use parcomm_mod,                 only : parcomm_global

implicit none

contains

subroutine create_hor_Christofel(hor_Christofel,mesh_u,mesh_v,metric,vector_type)
    class(hor_Christofel_t), allocatable, intent(out) :: hor_Christofel
    type(mesh_t),     intent(in) :: mesh_u, mesh_v
    class(metric_t),  intent(in) :: metric
    character(len=*), intent(in) :: vector_type

    call create_hor_Christofel_shallow(hor_Christofel,mesh_u,mesh_v,metric,vector_type)
end subroutine

subroutine create_hor_Christofel_shallow(hor_Christofel,mesh_u,mesh_v,metric,vector_type)
    class(hor_Christofel_t), allocatable, intent(out) :: hor_Christofel
    type(mesh_t),     intent(in) :: mesh_u, mesh_v
    class(metric_t),  intent(in) :: metric
    character(len=*), intent(in) :: vector_type

    type(hor_Christofel_shallow_t), allocatable :: Christofel

    allocate(Christofel)

    Christofel%vector_type = vector_type
    call init_Christofel_terms(Christofel%Gu11,Christofel%Gu12,Christofel%Gu21, &
                               Christofel%Gu22,mesh_u,1,vector_type, metric)

    call init_Christofel_terms(Christofel%Gv11,Christofel%Gv12,Christofel%Gv21, &
                               Christofel%Gv22,mesh_v,2,vector_type, metric)

    if(vector_type == "covariant") then
        call init_Christofel_S(Christofel%Su11, Christofel%Su12, Christofel%Su22, mesh_u, 1, metric)
        call init_Christofel_S(Christofel%Sv11, Christofel%Sv12, Christofel%Sv22, mesh_v, 2, metric)
    end if

    call move_alloc(Christofel,hor_Christofel)
end subroutine

subroutine init_Christofel_terms(g11,g12,g21,g22,mesh,comp_num,vector_type, metric)
    type(grid_field_t), intent(inout) :: g11, g12, g21, g22
    type(mesh_t),       intent(in)    :: mesh
    integer(kind=4),    intent(in)    :: comp_num
    character(len=*),   intent(in)    :: vector_type
    class(metric_t),    intent(in)    :: metric

    real(kind=8) :: G(3,3,3), alpha, beta
    integer(kind=4) :: t, i, j, k, ks, ke

    ! select type(metric)
    ! class is (deep_atm)
    !     call parcomm_global%abort("hor_Christofel_factory_mod, unknown metric type: "//
    !                               metric_type)
    ! class default ! treat as shallow atm
        call create_grid_field_2d(g11,0,mesh)
        call create_grid_field_2d(g12,0,mesh)
        call create_grid_field_2d(g21,0,mesh)
        call create_grid_field_2d(g22,0,mesh)
        ks = 1; ke = 1
    ! end select

    select case(vector_type)
    case("contravariant")
        do t=mesh%ts, mesh%te; do k = ks, ke
            do j = mesh%tile(t)%js, mesh%tile(t)%je
                do i = mesh%tile(t)%is, mesh%tile(t)%ie
                    alpha = mesh%tile(t)%get_alpha(i)
                    beta = mesh%tile(t)%get_beta(j)
                    G(1:3,1:3,1:3) = metric%calculate_G(mesh%tile(t)%panel_ind, &
                                                                    alpha, beta)
                    g11%tile(t)%p(i,j,k) = G(1,1,comp_num)
                    g12%tile(t)%p(i,j,k) = G(1,2,comp_num)
                    g21%tile(t)%p(i,j,k) = G(2,1,comp_num)
                    g22%tile(t)%p(i,j,k) = G(2,2,comp_num)
                end do
            end do
        end do; end do
    case("covariant")
        do t=mesh%ts, mesh%te; do k = ks, ke
            do j = mesh%tile(t)%js, mesh%tile(t)%je
                do i = mesh%tile(t)%is, mesh%tile(t)%ie
                    alpha = mesh%tile(t)%get_alpha(i)
                    beta = mesh%tile(t)%get_beta(j)
                    G(1:3,1:3,1:3) = metric%calculate_G(mesh%tile(t)%panel_ind, &
                                                                    alpha, beta)
                    g11%tile(t)%p(i,j,k) = G(comp_num,1,1)
                    g12%tile(t)%p(i,j,k) = G(comp_num,1,2)
                    g21%tile(t)%p(i,j,k) = G(comp_num,2,1)
                    g22%tile(t)%p(i,j,k) = G(comp_num,2,2)
                end do
            end do
        end do; end do
    case default
        call parcomm_global%abort("hor_Christofel_factory_mod, unknown vector_type: "// vector_type)
    end select

end subroutine

subroutine init_Christofel_S(S11, S12, S22, mesh, comp_num, metric)

    type(grid_field_t), intent(out) :: S11, S12, S22
    type(mesh_t),       intent(in)  :: mesh
    integer(kind=4),    intent(in)  :: comp_num
    class(metric_t),    intent(in)  :: metric

    integer(kind=4) :: i, j, t
    real(kind=8)    :: alpha, beta, S(3,2)

    call create_grid_field_2d(S11,0,mesh)
    call create_grid_field_2d(S12,0,mesh)
    call create_grid_field_2d(S22,0,mesh)

    do t=mesh%ts, mesh%te; 
        do j = mesh%tile(t)%js, mesh%tile(t)%je
            do i = mesh%tile(t)%is, mesh%tile(t)%ie

                alpha = mesh%tile(t)%get_alpha(i)
                beta = mesh%tile(t)%get_beta(j)
                S(1:3,1:2) = metric%calculate_S(mesh%tile(t)%panel_ind, &
                                                              alpha, beta)
                S11%tile(t)%p(i,j,1) = S(1,comp_num)
                S12%tile(t)%p(i,j,1) = S(2,comp_num)
                S22%tile(t)%p(i,j,1) = S(3,comp_num)
            end do
        end do
    end do


end subroutine init_Christofel_S

end module hor_Christofel_factory_mod
