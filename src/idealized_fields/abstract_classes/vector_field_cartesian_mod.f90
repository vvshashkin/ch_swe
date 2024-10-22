module vector_field_cartesian_mod

!shallow atmosphere vector field driver using cartesian coordinates

use vector_idealized_field_mod, only : vector_idealized_field_t
use mesh_mod,                   only : mesh_t, tile_mesh_t
use grid_field_mod,             only : grid_field_t, tile_field_t
use parcomm_mod,                only : parcomm_global
use sph_coords_mod,             only : cart2sph, sph2cart_vec

implicit none

type, abstract, extends(vector_idealized_field_t) :: vector_field_cartesian_t

    contains
        procedure, public :: get_cartesian_vec_tile
        procedure(get_cartesian_vec), deferred, public :: get_cartesian_vec
end type

abstract interface
    pure function get_cartesian_vec(this,x,y,z,h,klev,time) result(v)
        import vector_field_cartesian_t
        class(vector_field_cartesian_t), intent(in) :: this
        real(kind=8),                    intent(in) :: x, y, z, h, time
        integer(kind=4),                 intent(in) :: klev

        !Shallow atmosphere case: v(1:3) are vx,vy,vz of 'horizontal sphere', v(4) is vertical component
        !Deep atmosphere case: v(1:3) are x,y,z components of true 3d vector, v(4) should be 0
        real(kind=8) :: v(4)
    end function
end interface

contains

subroutine get_cartesian_vec_tile(this,v,mesh,halo_width,time)

    class(vector_field_cartesian_t), intent(in)  :: this
    type(tile_mesh_t),               intent(in)  :: mesh
    integer(kind=4),                 intent(in)  :: halo_width
    real(kind=8),                    intent(in)  :: time

    real(kind=8),                    intent(out) :: v(1:4,mesh%is-halo_width:mesh%ie+halo_width, &
                                                          mesh%js-halo_width:mesh%je+halo_width, &
                                                          mesh%ks-halo_width:mesh%ke+halo_width)

    integer(kind=4) :: i, j, k

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                v(1:4,i,j,k) = this%get_cartesian_vec(mesh%rx(i,j,k),mesh%ry(i,j,k),&
                                                      mesh%rz(i,j,k),mesh%h(i,j,k),k,time)
            end do
        end do
    end do

end subroutine

end module
