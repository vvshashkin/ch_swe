module vector_idealized_field_mod

!enchancing idealized_field_t interface for setting one vector-field component

use idealized_field_mod, only : idealized_field_t
use mesh_mod,            only : mesh_t, tile_mesh_t
use grid_field_mod,      only : grid_field_t, tile_field_t
use parcomm_mod,         only : parcomm_global

implicit none

type, abstract, extends(idealized_field_t) :: vector_idealized_field_t

    character(len=:), allocatable :: component_direction, component_type

    contains
        procedure, public                              :: get_field
        procedure, public                              :: get_field_tile
        procedure(get_cartesian_vec), deferred, public :: get_cartesian_vec_tile
end type

abstract interface
    subroutine get_cartesian_vec(this,v,mesh,halo_width,time)

        import vector_idealized_field_t, tile_mesh_t

        class(vector_idealized_field_t), intent(in)  :: this
        type(tile_mesh_t),               intent(in)  :: mesh
        integer(kind=4),                 intent(in)  :: halo_width
        real(kind=8),                    intent(in)  :: time

        !Shallow atmosphere case: v(1:3) are vx,vy,vz of 'horizontal sphere', v(4) is vertical component
        !Deep atmosphere case: v(1:3) are x,y,z components of true 3d vector, v(4) should be 0
        real(kind=8),                    intent(out) :: v(1:4,mesh%is-halo_width:mesh%ie+halo_width, &
                                                              mesh%js-halo_width:mesh%je+halo_width, &
                                                              mesh%ks-halo_width:mesh%ke+halo_width)

    end subroutine
end interface

contains

subroutine get_field(this, field, mesh, halo_width, fill_value, time)

    class(vector_idealized_field_t), intent(in)    :: this
    type(grid_field_t),              intent(inout) :: field
    type(mesh_t),                    intent(in)    :: mesh
    integer(kind=4),                 intent(in), optional :: halo_width
    real(kind=8),                    intent(in), optional :: fill_value, time

    integer(kind=4) :: t, hw

    hw = 0
    if(present(halo_width)) hw = halo_width

    do t = mesh%ts, mesh%te
        call this%get_field_tile(field%tile(t),mesh%tile(t),hw,fill_value,time)
    end do

end subroutine

subroutine get_field_tile(this, field, mesh, halo_width, fill_value,time)

    class(vector_idealized_field_t), intent(in)           :: this
    type(tile_field_t),              intent(inout)        :: field
    type(tile_mesh_t),               intent(in), target   :: mesh
    integer(kind=4),                 intent(in)           :: halo_width
    real(kind=8),                    intent(in), optional :: fill_value, time

    integer(kind=4) :: i, j, k, isv, iev, jsv, jev, ksv, kev, ncomp
    real(kind=8), pointer :: basis_vec(:,:,:,:)
    real(kind=8) :: v(1:4, mesh%is-halo_width:mesh%ie+halo_width, &
                           mesh%js-halo_width:mesh%je+halo_width, &
                           mesh%ks-halo_width:mesh%ke+halo_width)
   real(kind=8)  :: time_loc

   time_loc = this%time_default
   if(present(time)) time_loc = time

    call this%get_cartesian_vec_tile(v,mesh,halo_width,time_loc)

    if(present(fill_value)) then
        isv = field%is; iev = field%ie
        jsv = field%js; jev = field%je
        ksv = field%ks; kev = field%ke
        field%p(isv:iev,jsv:jev,ksv:kev) = fill_value
    end if

    if(this%component_direction == "x" .and. this%component_type == "covariant") then
        basis_vec => mesh%a1
    else if(this%component_direction == "y" .and. this%component_type == "covariant") then
        basis_vec => mesh%a2
    else if(this%component_direction == "z" .and. &
            (this%component_type == "covariant" .or. this%component_type == "real")) then
        basis_vec => mesh%a3
    else if(this%component_direction == "x" .and. this%component_type == "contravariant") then
        basis_vec => mesh%b1
    else if(this%component_direction == "y" .and. this%component_type == "contravariant") then
        basis_vec => mesh%b2
    else if(this%component_direction == "z" .and. this%component_type == "contravariant") then
        basis_vec => mesh%b3
    else
        call parcomm_global%abort(__FILE__//": get_field_tile, unknown combination of component_direction and component type: "// this%component_direction//" "//this%component_type)
    end if

    ncomp = size(basis_vec,1)
    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                field%p(i,j,k) = sum(v(1:ncomp,i,j,k)*basis_vec(1:ncomp,i,j,k))
            end do
        end do
    end do

    if(this%component_type == "real") then
        do k = mesh%ks, mesh%ke
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    field%p(i,j,k) = field%p(i,j,k) / sqrt(sum(basis_vec(1:ncomp,i,j,k)**2))
                end do
            end do
        end do
    end if

end subroutine

end module
