module orography_test_field_mod

use test_fields_3d_mod,   only : scalar_field3d_t, vector_field3d_t
use grid_field_mod,       only : tile_field_t
use mesh_mod,             only : tile_mesh_t
use const_mod,            only : pi
use sph_coords_mod,       only : cart2sph, sph2cart_vec
use const_mod,            only : Earth_radii


implicit none

type, extends(vector_field3d_t) :: orography_test_grad_t
    real(kind=8) :: h0, scale
    contains
    procedure :: get_vector_component_tile
end type orography_test_grad_t

type, extends(scalar_field3d_t) :: orography_test_field_t
    real(kind=8) :: h0
    contains
    procedure :: get_scalar_field_tile
end type orography_test_field_t

contains

subroutine get_vector_component_tile(this,v,mesh,halo_width, &
                                     base_vec, n_comp)

    class(orography_test_grad_t), intent(in)    :: this
    type(tile_field_t),           intent(inout) :: v
    type(tile_mesh_t),            intent(in)    :: mesh
    integer(kind=4),              intent(in)    :: halo_width
    real(kind=8),                 intent(in)    :: base_vec(n_comp, &
                                                     mesh%is-mesh%halo_width:mesh%ie+mesh%halo_width,&
                                                     mesh%js-mesh%halo_width:mesh%je+mesh%halo_width,&
                                                     mesh%ks:mesh%ke)
    integer(kind=4),         intent(in)    :: n_comp

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke
    real(kind=8)    :: velocity(4)
    real(kind=8)    :: lam, phi, u_sph, v_sph

    is = mesh%is-halo_width; ie = mesh%ie+halo_width
    js = mesh%js-halo_width; je = mesh%je+halo_width
    ks = mesh%ks; ke = mesh%ke

    do k=ks,ke
        do j=js,je
            do i=is,ie
                call cart2sph(mesh%rx(i,j,k), mesh%ry(i,j,k), mesh%rz(i,j,k), lam, phi)
                u_sph =-this%h0*2.0_8*sin(2.0_8*lam)*cos(phi) / this%scale
                v_sph =-this%h0*cos(2.0_8*lam)*sin(2.0_8*phi) / this%scale
                call sph2cart_vec(lam, phi, u_sph, v_sph, velocity(1), velocity(2), velocity(3))
                velocity(4) = 0.0_8
                v%p(i,j,k) = sum(velocity(1:n_comp)*base_vec(1:n_comp,i,j,k))
            end do
        end do
    end do
end subroutine get_vector_component_tile

subroutine get_scalar_field_tile(this,f,mesh,halo_width)
    class(orography_test_field_t), intent(in)    :: this
    type(tile_field_t),            intent(inout) :: f
    type(tile_mesh_t),             intent(in)    :: mesh
    integer(kind=4),               intent(in)    :: halo_width

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke
    real(kind=8)    :: lam, phi

    is = mesh%is-halo_width; ie = mesh%ie+halo_width
    js = mesh%js-halo_width; je = mesh%je+halo_width
    ks = mesh%ks; ke = mesh%ke

    do k=ks,ke
        do j=js,je
            do i=is,ie
                call cart2sph(mesh%rx(i,j,k), mesh%ry(i,j,k), mesh%rz(i,j,k), lam, phi)
                f%p(i,j,k) = this%h0*cos(2.0_8*lam)*cos(phi)**2
            end do
        end do
    end do

end subroutine get_scalar_field_tile

end module orography_test_field_mod
