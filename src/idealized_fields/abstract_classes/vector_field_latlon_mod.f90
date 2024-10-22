module vector_field_latlon_mod

!shallow atmosphere latlon vector field driver

use vector_idealized_field_mod, only : vector_idealized_field_t
use mesh_mod,                   only : mesh_t, tile_mesh_t
use grid_field_mod,             only : grid_field_t, tile_field_t
use parcomm_mod,                only : parcomm_global
use sph_coords_mod,             only : cart2sph, sph2cart_vec

implicit none

type, abstract, extends(vector_idealized_field_t) :: vector_field_latlon_t

    contains
        procedure, public :: get_cartesian_vec_tile
        procedure(get_latlon_vec), deferred, public :: get_latlon_vec
end type

abstract interface
    pure function get_latlon_vec(this,phi,lambda,h,klev,time) result(uvw)
        import vector_field_latlon_t
        class(vector_field_latlon_t), intent(in) :: this
        real(kind=8),                 intent(in) :: lambda, phi, h, time
        integer(kind=4),              intent(in) :: klev

        real(kind=8) :: uvw(3)
    end function
end interface

contains

subroutine get_cartesian_vec_tile(this,v,mesh,halo_width,time)

    class(vector_field_latlon_t), intent(in)  :: this
    type(tile_mesh_t),            intent(in)  :: mesh
    integer(kind=4),              intent(in)  :: halo_width
    real(kind=8),                 intent(in)  :: time

    real(kind=8),                 intent(out) :: v(1:4,mesh%is-halo_width:mesh%ie+halo_width, &
                                                       mesh%js-halo_width:mesh%je+halo_width, &
                                                       mesh%ks-halo_width:mesh%ke+halo_width)

    integer(kind=4) :: i, j, k
    real(kind=8)    :: lambda, phi, uvw(3)

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                call cart2sph(mesh%rx(i,j,k),mesh%ry(i,j,k),mesh%rz(i,j,k),lambda,phi)
                uvw(1:3) = this%get_latlon_vec(phi,lambda,mesh%h(i,j,k),k,time)
                call sph2cart_vec(lambda,phi,uvw(1),uvw(2),v(1,i,j,k),v(2,i,j,k),v(3,i,j,k))
                v(4,i,j,k) = uvw(3)
            end do
        end do
    end do

end subroutine

end module
