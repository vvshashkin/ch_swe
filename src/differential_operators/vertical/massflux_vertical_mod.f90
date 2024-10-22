module massflux_vertical_mod

use abstract_massflux_vertical_mod, only : massflux_vertical_t
use grid_field_mod,                 only : tile_field_t
use mesh_mod,                       only : tile_mesh_t

implicit none

type, extends(massflux_vertical_t) :: massflux_vert_up1_t
    contains
        procedure, public :: calc_vertical_massflux_tile => calc_vertical_massflux_tile_up1
end type

type, extends(massflux_vertical_t) :: massflux_vert_up4_t
    contains
        procedure, public :: calc_vertical_massflux_tile => calc_vertical_massflux_tile_up4
end type

type, extends(massflux_vertical_t) :: massflux_vert_c2_t
    contains
        procedure, public :: calc_vertical_massflux_tile => calc_vertical_massflux_tile_c2
end type

contains

subroutine calc_vertical_massflux_tile_up1(this, massflux, f, eta_dot, mesh)
    class(massflux_vert_up1_t),  intent(in)    :: this
    type(tile_field_t),          intent(in)    :: f, eta_dot
    type(tile_mesh_t),           intent(in)    :: mesh
    !output
    type(tile_field_t),          intent(inout) :: massflux

    integer(kind=4) :: i,j,k
    real(kind=8)    :: ap, am

    if(mesh%ks <=1) massflux%p(mesh%is:mesh%ie,mesh%js:mesh%je,1) = 0.0_8

    call calc_up1_flux(massflux,f,eta_dot,mesh,max(mesh%ks,2),min(mesh%ke,mesh%nz-1))

    if(mesh%ke >= mesh%nz) massflux%p(mesh%is:mesh%ie,mesh%js:mesh%je,mesh%nz) = 0.0_8

end subroutine

subroutine calc_vertical_massflux_tile_up4(this, massflux, f, eta_dot, mesh)
    class(massflux_vert_up4_t),  intent(in)    :: this
    type(tile_field_t),          intent(in)    :: f, eta_dot
    type(tile_mesh_t),           intent(in)    :: mesh
    !output
    type(tile_field_t),          intent(inout) :: massflux

    integer(kind=4) :: i,j,k
    real(kind=8)    :: ap, am

    if(mesh%ks <=1) massflux%p(mesh%is:mesh%ie,mesh%js:mesh%je,1) = 0.0_8

    if(mesh%ks<=2.and.mesh%ke>=2) &
        call calc_up1_flux(massflux,f,eta_dot,mesh,2,2)

    if(mesh%ks<=3.and.mesh%ke>=3) &
        call calc_up3_flux(massflux,f,eta_dot,mesh,3,3)

    call calc_up4_flux(massflux,f,eta_dot,mesh,max(4,mesh%ks),min(mesh%ke,mesh%nz-3))

    if(mesh%ks<=mesh%nz-2.and.mesh%ke>=mesh%nz-2) &
        call calc_up3_flux(massflux,f,eta_dot,mesh,mesh%nz-2,mesh%nz-2)

    if(mesh%ks<=mesh%nz-1.and.mesh%ke>=mesh%nz-1) &
        call calc_up1_flux(massflux,f,eta_dot,mesh,mesh%nz-1,mesh%nz-1)

    if(mesh%ke >= mesh%nz) massflux%p(mesh%is:mesh%ie,mesh%js:mesh%je,mesh%nz) = 0.0_8

end subroutine

subroutine calc_vertical_massflux_tile_c2(this, massflux, f, eta_dot, mesh)
    class(massflux_vert_c2_t),   intent(in)    :: this
    type(tile_field_t),          intent(in)    :: f, eta_dot
    type(tile_mesh_t),           intent(in)    :: mesh
    !output
    type(tile_field_t),          intent(inout) :: massflux

    integer(kind=4) :: i,j,k
    real(kind=8)    :: ap, am

    if(mesh%ks <= 1) massflux%p(mesh%is:mesh%ie,mesh%js:mesh%je,1) = 0.0_8

    call calc_c2_flux(massflux,f,eta_dot,mesh,max(2,mesh%ks),min(mesh%ke,mesh%nz-1))

    if(mesh%ke >= mesh%nz) massflux%p(mesh%is:mesh%ie,mesh%js:mesh%je,mesh%nz) = 0.0_8

end subroutine

subroutine calc_up1_flux(massflux, f, eta_dot, mesh, k1, k2)
    type(tile_field_t),          intent(in)    :: f, eta_dot
    type(tile_mesh_t),           intent(in)    :: mesh
    integer(kind=4),             intent(in)    :: k1, k2
    !output
    type(tile_field_t),          intent(inout) :: massflux

    integer(kind=4) :: i,j,k
    real(kind=8)    :: ap, am

    do k = k1, k2
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                ap = 0.5_8+sign(0.5_8,eta_dot%p(i,j,k))
                am = 1.0_8-ap
                massflux%p(i,j,k) = eta_dot%p(i,j,k)*(ap*f%p(i,j,k-1)+am*f%p(i,j,k))
            end do
        end do
    end do

end subroutine

subroutine calc_c2_flux(massflux, f, eta_dot, mesh, k1, k2)
    type(tile_field_t),          intent(in)    :: f, eta_dot
    type(tile_mesh_t),           intent(in)    :: mesh
    integer(kind=4),             intent(in)    :: k1, k2
    !output
    type(tile_field_t),          intent(inout) :: massflux

    integer(kind=4) :: i,j,k
    real(kind=8)    :: ap, am

    do k = k1, k2
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                massflux%p(i,j,k) = 0.5_8*eta_dot%p(i,j,k)*(f%p(i,j,k-1)+f%p(i,j,k))
            end do
        end do
    end do

end subroutine

subroutine calc_up3_flux(massflux, f, eta_dot, mesh, k1, k2)
    type(tile_field_t),          intent(in)    :: f, eta_dot
    type(tile_mesh_t),           intent(in)    :: mesh
    integer(kind=4),             intent(in)    :: k1, k2
    !output
    type(tile_field_t),          intent(inout) :: massflux

    integer(kind=4) :: i,j,k
    real(kind=8)    :: ap, am, fp, fm

    do k = k1, k2
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                ap = 0.5_8+sign(0.5_8,eta_dot%p(i,j,k))
                am = 1.0_8-ap
                fp = (3._8*f%p(i,j,k)  +6._8*f%p(i,j,k-1)-f%p(i,j,k-2))
                fm = (3._8*f%p(i,j,k-1)+6._8*f%p(i,j,k)  -f%p(i,j,k+1))
                massflux%p(i,j,k) = 0.125_8*eta_dot%p(i,j,k)*(ap*fp+am*fm)
            end do
        end do
    end do

end subroutine

subroutine calc_up4_flux(massflux, f, eta_dot, mesh, k1, k2)
    type(tile_field_t),          intent(in)    :: f, eta_dot
    type(tile_mesh_t),           intent(in)    :: mesh
    integer(kind=4),             intent(in)    :: k1, k2
    !output
    type(tile_field_t),          intent(inout) :: massflux

    integer(kind=4) :: i,j,k
    real(kind=8)    :: ap, am, fp, fm

    do k = k1, k2
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                ap = 0.5_8+sign(0.5_8,eta_dot%p(i,j,k))
                am = 1.0_8-ap
                fp = (5._8*f%p(i,j,k)  +15._8*f%p(i,j,k-1)-5._8*f%p(i,j,k-2)+f%p(i,j,k-3))
                fm = (5._8*f%p(i,j,k-1)+15._8*f%p(i,j,k)  -5._8*f%p(i,j,k+1)+f%p(i,j,k+2))
                massflux%p(i,j,k) = 0.0625_8*eta_dot%p(i,j,k)*(ap*fp+am*fm)
            end do
        end do
    end do

end subroutine

end module
