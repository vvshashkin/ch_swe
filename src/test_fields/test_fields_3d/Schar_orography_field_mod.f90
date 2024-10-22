module Schar_orography_field_mod

use test_fields_3d_mod,   only : scalar_field3d_t, vector_field3d_t
use grid_field_mod,       only : tile_field_t
use mesh_mod,             only : tile_mesh_t
use const_mod,            only : pi
use sph_coords_mod,       only : cart2sph, sph2cart_vec
use const_mod,            only : Earth_radii, pi

implicit none

real(kind=8), parameter :: phi0 = 0.0_8, lam0 = 0.25_8*pi
real(kind=8), parameter :: d1 = 5000.0_8, d2 = 4000.0_8

type, extends(vector_field3d_t) :: Schar_orography_grad_t
    real(kind=8) :: h0, scale
    contains
    procedure :: get_vector_component_tile
end type Schar_orography_grad_t

type, extends(scalar_field3d_t) :: Schar_orography_field_t
    real(kind=8) :: h0, scale
    contains
    procedure :: get_scalar_field_tile
end type Schar_orography_field_t

contains

subroutine get_vector_component_tile(this,v,mesh,halo_width, &
                                     base_vec, n_comp)

    class(Schar_orography_grad_t), intent(in)    :: this
    type(tile_field_t),            intent(inout) :: v
    type(tile_mesh_t),             intent(in)    :: mesh
    integer(kind=4),               intent(in)    :: halo_width
    real(kind=8),                  intent(in)    :: base_vec(n_comp, &
                                                     mesh%is-mesh%halo_width:mesh%ie+mesh%halo_width,&
                                                     mesh%js-mesh%halo_width:mesh%je+mesh%halo_width,&
                                                     mesh%ks:mesh%ke)
    integer(kind=4),         intent(in)    :: n_comp

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke
    real(kind=8)    :: grad(4), x, y, z, cos_psi, r, x0, y0, z0, dh_dd, d

    is = mesh%is-halo_width; ie = mesh%ie+halo_width
    js = mesh%js-halo_width; je = mesh%je+halo_width
    ks = mesh%ks; ke = mesh%ke

    x0 = cos(phi0)*cos(lam0)
    y0 = cos(phi0)*sin(lam0)
    z0 = sin(phi0)

    do k=ks,ke
        do j=js,je
            do i=is,ie
                x = mesh%rx(i,j,k); y = mesh%ry(i,j,k); z = mesh%rz(i,j,k)
                r = sqrt(x**2+y**2+z**2)
                x = x/r; y = y/r; z = z/r
                cos_psi = x*x0+y*y0+z*z0
                d = this%scale*acos(cos_psi)
                dh_dd =-2.0_8*this%h0*exp(-d**2/d1**2)*(d/d1**2*cos(pi*d/d2)**2+&
                                                        pi/d2*cos(pi*d/d2)*sin(pi*d/d2))
                r = max(1e-14,sqrt(1.0_8-cos_psi**2))
                grad(1) =-dh_dd / r * (x0-cos_psi*x)
                grad(2) =-dh_dd / r * (y0-cos_psi*y)
                grad(3) =-dh_dd / r * (z0-cos_psi*z)
                grad(4) = 0.0_8
                v%p(i,j,k) = sum(grad(1:n_comp)*base_vec(1:n_comp,i,j,k))
            end do
        end do
    end do
end subroutine get_vector_component_tile

subroutine get_scalar_field_tile(this,f,mesh,halo_width)
    class(Schar_orography_field_t), intent(in)    :: this
    type(tile_field_t),             intent(inout) :: f
    type(tile_mesh_t),              intent(in)    :: mesh
    integer(kind=4),                intent(in)    :: halo_width

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke
    real(kind=8)    :: rx0, ry0, rz0, d

    is = mesh%is-halo_width; ie = mesh%ie+halo_width
    js = mesh%js-halo_width; je = mesh%je+halo_width
    ks = mesh%ks; ke = mesh%ke

    rx0 = cos(phi0)*cos(lam0)
    ry0 = cos(phi0)*sin(lam0)
    rz0 = sin(phi0)

    do k=ks,ke
        do j=js,je
            do i=is,ie
                d = this%scale*acos(mesh%rx(i,j,k)*rx0+mesh%ry(i,j,k)*ry0+mesh%rz(i,j,k)*rz0)
                f%p(i,j,k) = this%h0*exp(-d**2/d1**2)*cos(pi*d/d2)**2
            end do
        end do
    end do

end subroutine get_scalar_field_tile

end module Schar_orography_field_mod
