module baroclinic_instability_testcase_mod

use test_fields_3d_mod,                         only : vector_field3d_t, scalar_field3d_t
use mesh_mod,                                   only : tile_mesh_t
use grid_field_mod,                             only : tile_field_t
use const_mod,                                  only : pi, Day24h_sec
use baroclinic_instability_solve_H_sigma_mod,   only : solve_h_for_sigma
use baroclinic_instability_test_parameters_mod, only : u0, sigma_0, xc, yc, zc
use baroclinic_instability_test_therm_mod,      only : calc_theta, calc_Pexner, calc_Temp

implicit none

type, extends(vector_field3d_t) :: baroclinic_instability_test_wind_z_t
    real(kind=8) :: u_pert
    contains
    procedure :: get_vector_component_tile
end type baroclinic_instability_test_wind_z_t

type, extends(scalar_field3d_t) :: baroclinic_instability_test_PExner_z_t
    contains
    procedure :: get_scalar_field_tile => get_Pexner
end type

type, extends(scalar_field3d_t) :: baroclinic_instability_test_theta_z_t
    contains
    procedure :: get_scalar_field_tile => get_theta
end type

!Isobaric temperature field, originally introduced to test postprocessing
type, extends(scalar_field3d_t) :: baroclinic_instability_test_isobaric_T_t
    real(kind=8) :: p
    contains
    procedure :: get_scalar_field_tile => get_isobaric_T
end type

contains

subroutine get_vector_component_tile(this,v,mesh,halo_width, &
                                     base_vec, n_comp)

    class(baroclinic_instability_test_wind_z_t),  intent(in)    :: this
    type(tile_field_t),                           intent(inout) :: v
    type(tile_mesh_t),                            intent(in)    :: mesh
    integer(kind=4),                              intent(in)    :: halo_width
    real(kind=8), intent(in)    :: base_vec(n_comp, &
                                            mesh%is-mesh%halo_width:mesh%ie+mesh%halo_width,&
                                            mesh%js-mesh%halo_width:mesh%je+mesh%halo_width,&
                                            mesh%ks:mesh%ke)
    integer(kind=4),         intent(in)    :: n_comp

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke
    real(kind=8)    :: u, vs(3), sigma_v, sinphi, x, y, z, r, xy
    real(kind=8), dimension(mesh%is-halo_width:mesh%ie+halo_width,&
                            mesh%js-halo_width:mesh%je+halo_width,&
                            mesh%ks:mesh%ke) :: sigma, h
    real(kind=8)    :: phi(mesh%is-halo_width:mesh%ie+halo_width,&
                           mesh%js-halo_width:mesh%je+halo_width)

    is = mesh%is-halo_width; ie = mesh%ie+halo_width
    js = mesh%js-halo_width; je = mesh%je+halo_width
    ks = mesh%ks; ke = mesh%ke

    h(is:ie,js:je,ks:ke) = mesh%h(is:ie,js:je,ks:ke)
    phi(is:ie,js:je) = asin(mesh%rz(is:ie,js:je,ks))

    call solve_h_for_sigma(sigma,h,phi)

    do k=ks,ke
        do j=js,je
            do i=is,ie
                x = mesh%rx(i,j,k)
                y = mesh%ry(i,j,k)
                z = mesh%rz(i,j,k)
                sigma_v = (sigma(i,j,k)-sigma_0)*0.5_8*pi
                r = acos(xc*x+yc*y+zc*z)
                u = u0*cos(sigma_v)**1.5_8*sin(2.0_8*phi(i,j))**2 + &
                                       this%u_pert*exp(-100.0_8*r**2)
                xy = max(sqrt(x**2+y**2),1e-10)
                vs(1) = -y*u/xy; vs(2) = x*u/xy; vs(3) = 0.0_8
                v%p(i,j,k) = sum(vs(1:3)*base_vec(1:3,i,j,k))
            end do
        end do
    end do

end subroutine get_vector_component_tile

subroutine get_theta(this,f,mesh,halo_width)
    class(baroclinic_instability_test_theta_z_t), intent(in)    :: this
    type(tile_field_t),                           intent(inout) :: f
    type(tile_mesh_t),                            intent(in)    :: mesh
    integer(kind=4),                              intent(in)    :: halo_width

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke
    real(kind=8)    :: phi(mesh%is-halo_width:mesh%ie+halo_width, &
                           mesh%js-halo_width:mesh%je+halo_width)
    real(kind=8)    :: sigma(mesh%is-halo_width:mesh%ie+halo_width,&
                             mesh%js-halo_width:mesh%je+halo_width, &
                             mesh%ks:mesh%ke)
    real(kind=8)    :: h(mesh%is-halo_width:mesh%ie+halo_width,&
                         mesh%js-halo_width:mesh%je+halo_width, &
                         mesh%ks:mesh%ke)

    is = mesh%is-halo_width; ie = mesh%ie+halo_width
    js = mesh%js-halo_width; je = mesh%je+halo_width
    ks = mesh%ks; ke = mesh%ke

    h(is:ie,js:je,ks:ke) = mesh%h(is:ie,js:je,ks:ke)
    phi(is:ie,js:je) = asin(mesh%rz(is:ie,js:je,ks))
    call solve_h_for_sigma(sigma,h,phi)
    call calc_theta(h,sigma,phi)
    f%p(is:ie,js:je,ks:ke) = h(is:ie,js:je,ks:ke)

end subroutine get_theta

subroutine get_Pexner(this,f,mesh,halo_width)
    class(baroclinic_instability_test_Pexner_z_t), intent(in)    :: this
    type(tile_field_t),                            intent(inout) :: f
    type(tile_mesh_t),                             intent(in)    :: mesh
    integer(kind=4),                               intent(in)    :: halo_width

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke
    real(kind=8)    :: phi(mesh%is-halo_width:mesh%ie+halo_width, &
                           mesh%js-halo_width:mesh%je+halo_width)
    real(kind=8)    :: sigma(mesh%is-halo_width:mesh%ie+halo_width,&
                             mesh%js-halo_width:mesh%je+halo_width, &
                             mesh%ks:mesh%ke)
    real(kind=8)    :: h(mesh%is-halo_width:mesh%ie+halo_width,&
                         mesh%js-halo_width:mesh%je+halo_width, &
                         mesh%ks:mesh%ke)

    is = mesh%is-halo_width; ie = mesh%ie+halo_width
    js = mesh%js-halo_width; je = mesh%je+halo_width
    ks = mesh%ks; ke = mesh%ke

    h(is:ie,js:je,ks:ke) = mesh%h(is:ie,js:je,ks:ke)
    phi(is:ie,js:je) = asin(mesh%rz(is:ie,js:je,ks))
    call solve_h_for_sigma(sigma,h,phi)
    f%p(is:ie,js:je,ks:ke) = calc_Pexner(sigma(is:ie,js:je,ks:ke))

end subroutine get_Pexner

subroutine get_isobaric_T(this,f,mesh,halo_width)
    class(baroclinic_instability_test_isobaric_T_t), intent(in)    :: this
    type(tile_field_t),                              intent(inout) :: f
    type(tile_mesh_t),                               intent(in)    :: mesh
    integer(kind=4),                                 intent(in)    :: halo_width

    real(kind=8)    :: phi(mesh%is:mesh%ie, mesh%js:mesh%je)
    real(kind=8)    :: sigma(mesh%is:mesh%ie, mesh%js:mesh%je,1)
    real(kind=8)    :: temp(mesh%is:mesh%ie,mesh%js:mesh%je,1)

    integer(kind=4) :: k,is,ie,js,je,ks,ke

    is = mesh%is; ie = mesh%ie; js = mesh%js; je = mesh%je; ks = mesh%ks; ke = mesh%ke

    phi(is:ie,js:je)   = asin(mesh%rz(is:ie,js:je,ks))
    sigma(is:ie,js:je,1) = this%p / 1e5_8
    call calc_Temp(temp,sigma,phi)
    do k=ks,ke
        f%p(is:ie,js:je,k) = temp(is:ie,js:je,1)
    end do

end subroutine get_isobaric_T

end module baroclinic_instability_testcase_mod
