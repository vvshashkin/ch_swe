module baroclinic_instability_test_orography_mod

use test_fields_3d_mod,   only : scalar_field3d_t, vector_field3d_t
use grid_field_mod,       only : tile_field_t
use mesh_mod,             only : tile_mesh_t
use const_mod,            only : pi
use sph_coords_mod,       only : cart2sph, sph2cart_vec
use const_mod,            only : Earth_radii, pi, Earth_omega, Earth_grav

implicit none

real(kind=8), parameter :: sigma_0 = 0.252_8, u0 = 35.0_8

type, extends(vector_field3d_t) :: baroclinic_instability_test_orography_grad_t
    contains
    procedure :: get_vector_component_tile
end type

type, extends(scalar_field3d_t) :: baroclinic_instability_test_orography_t
    contains
    procedure :: get_scalar_field_tile
end type

contains

subroutine get_vector_component_tile(this,v,mesh,halo_width, &
                                     base_vec, n_comp)

    class(baroclinic_instability_test_orography_grad_t), intent(in)    :: this
    type(tile_field_t),             intent(inout) :: v
    type(tile_mesh_t),              intent(in)    :: mesh
    integer(kind=4),                intent(in)    :: halo_width
    real(kind=8),                   intent(in)    :: base_vec(n_comp, &
                                                     mesh%is-mesh%halo_width:mesh%ie+mesh%halo_width,&
                                                     mesh%js-mesh%halo_width:mesh%je+mesh%halo_width,&
                                                     mesh%ks:mesh%ke)
    integer(kind=4),         intent(in)    :: n_comp

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke
    real(kind=8)    :: grad(4), x, y, z
    real(kind=8)    :: sin_lat, cos_lat, lam, sigma_v, u, dA, dB, dh

    is = mesh%is-halo_width; ie = mesh%ie+halo_width
    js = mesh%js-halo_width; je = mesh%je+halo_width
    ks = mesh%ks; ke = mesh%ke

    do k=ks,ke
        do j=js,je
            do i=is,ie
                lam = atan2(mesh%ry(i,j,k),mesh%rx(i,j,k))
                sin_lat = mesh%rz(i,j,k)
                cos_lat = sqrt(1.0_8-sin_lat**2)
                sigma_v = 0.5_8*pi*(1.0_8-sigma_0)
                u = u0*cos(sigma_v)**(1.5_8)
                dA = (-12.0_8*sin_lat**5*cos_lat*(cos_lat**2+1.0_8/3.0_8)+&
                        4.0_8*sin_lat**7*cos_lat)
                dB = (-24.0_8/5.0_8*cos_lat**2*sin_lat*(sin_lat**2+2.0_8/3.0_8)+&
                       16.0_8/5.0_8*cos_lat**4*sin_lat)
                dh = u*(dA*u+dB*Earth_omega*Earth_radii)/(Earth_radii*Earth_grav)
                grad(1) = -cos(lam)*sin_lat*dh
                grad(2) = -sin(lam)*sin_lat*dh
                grad(3) = cos_lat*dh
                grad(4) = 0.0_8
                v%p(i,j,k) = sum(grad(1:n_comp)*base_vec(1:n_comp,i,j,k))
            end do
        end do
    end do
end subroutine get_vector_component_tile

subroutine get_scalar_field_tile(this,f,mesh,halo_width)
    class(baroclinic_instability_test_orography_t), intent(in)    :: this
    type(tile_field_t),              intent(inout) :: f
    type(tile_mesh_t),               intent(in)    :: mesh
    integer(kind=4),                 intent(in)    :: halo_width

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke
    real(kind=8)    :: sin_lat, cos_lat, sigma_v, u, A, B

    is = mesh%is-halo_width; ie = mesh%ie+halo_width
    js = mesh%js-halo_width; je = mesh%je+halo_width
    ks = mesh%ks; ke = mesh%ke

    do k=ks,ke
        do j=js,je
            do i=is,ie
                sin_lat = mesh%rz(i,j,k)
                cos_lat = sqrt(1.0_8-sin_lat**2)
                sigma_v = 0.5_8*pi*(1.0_8-sigma_0)
                u = u0*cos(sigma_v)**(1.5_8)
                A = (-2.0_8*sin_lat**6*(cos_lat**2+1.0_8/3.0_8)+10.0_8/63.0_8)
                B = (8.0_8/5.0_8*cos_lat**3*(sin_lat**2+2.0_8/3.0_8)-0.25_8*pi)
                f%p(i,j,k) = u*(A*u+B*Earth_radii*Earth_omega)/Earth_grav
            end do
        end do
    end do

end subroutine get_scalar_field_tile

end module baroclinic_instability_test_orography_mod
