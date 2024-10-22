module hor_Christofel_shallow_mod

use abstract_hor_Christofel_mod, only : hor_Christofel_t
use grid_field_mod,              only : grid_field_t, tile_field_t
use mesh_mod,                    only : mesh_t, tile_mesh_t
use parcomm_mod,                 only : parcomm_global

implicit none

type, extends(hor_Christofel_t) :: hor_Christofel_shallow_t
    character(len=:), allocatable :: vector_type
    type(grid_field_t) :: Gu11, Gu12, Gu21, Gu22
    type(grid_field_t) :: Gv11, Gv12, Gv21, Gv22
    type(grid_field_t) :: Su11, Su12, Su22
    type(grid_field_t) :: Sv11, Sv12, Sv22
contains
    procedure :: add_vector_advection
    procedure :: add_cov_wind_advection
    procedure :: add_contra_wind_advection
end type hor_Christofel_shallow_t

contains

subroutine add_vector_advection(this,tend,s1,s2,ut,vt,mesh,component)
    class(hor_Christofel_shallow_t), intent(inout) :: this
    type(grid_field_t),              intent(inout) :: tend
    type(grid_field_t),              intent(in)    :: s1, s2, ut, vt
    type(mesh_t),                    intent(in)    :: mesh
    character(len=*),                intent(in)    :: component

    integer(kind=4) :: t

    select case(component)
    case("u","1")
        do t=mesh%ts, mesh%te
            call add_vector_advection_tile(tend%tile(t),this%Gu11%tile(t), this%Gu12%tile(t),  &
                                           this%Gu21%tile(t), this%Gu22%tile(t),               &
                                           s1%tile(t), s2%tile(t),ut%tile(t), vt%tile(t),      &
                                           mesh%tile(t), mesh%scale)
            ! call add_vector_advection_tile(tend%tile(t),this%Gu11%tile(t), this%Gv11%tile(t),  &
            !                                this%Gu12%tile(t), this%Gv12%tile(t),               &
            !                                s1%tile(t), s2%tile(t),ut%tile(t), vt%tile(t),      &
            !                                mesh%tile(t), mesh%scale)
        end do
    case("v","2")
        do t=mesh%ts, mesh%te
            call add_vector_advection_tile(tend%tile(t),this%Gv11%tile(t), this%Gv12%tile(t),   &
                                           this%Gv21%tile(t), this%Gv22%tile(t),                &
                                           s1%tile(t), s2%tile(t), ut%tile(t), vt%tile(t),      &
                                           mesh%tile(t), mesh%scale)
            ! call add_vector_advection_tile(tend%tile(t),this%Gu12%tile(t), this%Gv12%tile(t),   &
            !                                this%Gu22%tile(t), this%Gv22%tile(t),                &
            !                                s1%tile(t), s2%tile(t), ut%tile(t), vt%tile(t),      &
            !                                mesh%tile(t), mesh%scale)
        end do
    case default
        call parcomm_global%abort("shallow hor Christofel, incorrect component: "// component)
    end select

end subroutine add_vector_advection

subroutine add_contra_wind_advection(this,tend,ut,vt,mesh,component)
    class(hor_Christofel_shallow_t), intent(inout) :: this
    type(grid_field_t),              intent(inout) :: tend
    type(grid_field_t),              intent(in)    :: ut, vt
    type(mesh_t),                    intent(in)    :: mesh
    character(len=*),                intent(in)    :: component

    integer(kind=4) :: t

    if(this%vector_type /= "contravariant") &
        call parcomm_global%abort("Christofel%add_contra_wind_advection is called, but vector type is "// this%vector_type)

    select case(component)
    case("u","1")
        do t=mesh%ts, mesh%te
            call add_contra_wind_advection_tile(tend%tile(t),this%Gu11%tile(t),       &
                                                this%Gu12%tile(t), this%Gu22%tile(t), &
                                                ut%tile(t), vt%tile(t), mesh%tile(t), mesh%scale)
        end do
    case("v","2")
        do t=mesh%ts, mesh%te
            call add_contra_wind_advection_tile(tend%tile(t),this%Gv11%tile(t),       &
                                                this%Gv12%tile(t), this%Gv22%tile(t), &
                                                ut%tile(t), vt%tile(t), mesh%tile(t), mesh%scale)
        end do
    case default
        call parcomm_global%abort("shallow hor Christofel, incorrect component: "// component)
    end select
end subroutine add_contra_wind_advection

subroutine add_cov_wind_advection(this,tend,ut,vt,mesh,component)
    class(hor_Christofel_shallow_t), intent(inout) :: this
    type(grid_field_t),              intent(inout) :: tend
    type(grid_field_t),              intent(in)    :: ut, vt
    type(mesh_t),                    intent(in)    :: mesh
    character(len=*),                intent(in)    :: component

    integer(kind=4) :: t

    if(this%vector_type /= "covariant") &
        call parcomm_global%abort("Christofel%add_cov_wind_advection is called, but vector type is "// this%vector_type)

    select case(component)
    case("u","1")
        do t=mesh%ts, mesh%te
            call add_contra_wind_advection_tile(tend%tile(t),this%Su11%tile(t),       &
                                                this%Su12%tile(t), this%Su22%tile(t), &
                                                ut%tile(t), vt%tile(t), mesh%tile(t), -mesh%scale)
        end do
    case("v","2")
        do t=mesh%ts, mesh%te
            call add_contra_wind_advection_tile(tend%tile(t),this%Sv11%tile(t),       &
                                                this%Sv12%tile(t), this%Sv22%tile(t), &
                                                ut%tile(t), vt%tile(t), mesh%tile(t), -mesh%scale)
        end do
    case default
        call parcomm_global%abort("shallow hor Christofel, incorrect component: "// component)
    end select
end subroutine add_cov_wind_advection

subroutine add_vector_advection_tile(tend,g11,g12,g21,g22,s1,s2,ut,vt,mesh,scale)
    type(tile_field_t), intent(inout) :: tend
    type(tile_field_t), intent(in)    :: g11, g12, g21, g22
    type(tile_field_t), intent(in)    :: s1, s2, ut, vt
    type(tile_mesh_t),  intent(in)    :: mesh
    real(kind=8),       intent(in)    :: scale

    integer(kind=4) :: i, j, k

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                tend%p(i,j,k) =  tend%p(i,j,k)+&
                                 (s1%p(i,j,k)*ut%p(i,j,k)*g11%p(i,j,1)+ &
                                  s2%p(i,j,k)*ut%p(i,j,k)*g12%p(i,j,1)+ &
                                  s1%p(i,j,k)*vt%p(i,j,k)*g21%p(i,j,1)+ &
                                  s2%p(i,j,k)*vt%p(i,j,k)*g22%p(i,j,1)) / scale
            end do
        end do
    end do
end subroutine add_vector_advection_tile

subroutine add_contra_wind_advection_tile(tend,g11,g21,g22,ut,vt,mesh,scale)
    type(tile_field_t), intent(inout) :: tend
    type(tile_field_t), intent(in)    :: g11, g21, g22
    type(tile_field_t), intent(in)    :: ut, vt
    type(tile_mesh_t),  intent(in)    :: mesh
    real(kind=8),       intent(in)    :: scale

    integer(kind=4) :: i, j, k

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                tend%p(i,j,k) = tend%p(i,j,k)-(ut%p(i,j,k)*ut%p(i,j,k)*g11%p(i,j,1)+ &
                                        2.0_8*ut%p(i,j,k)*vt%p(i,j,k)*g21%p(i,j,1)+ &
                                              vt%p(i,j,k)*vt%p(i,j,k)*g22%p(i,j,1)) / scale
            end do
        end do
    end do
end subroutine add_contra_wind_advection_tile

end module hor_Christofel_shallow_mod
