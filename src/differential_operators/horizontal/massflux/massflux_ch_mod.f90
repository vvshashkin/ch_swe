module massflux_ch_mod

use abstract_massflux_mod,        only : massflux_operator_t
use abstract_interpolators2d_mod, only : interpolator2d_vec2vec_t
use domain_mod,                   only : domain_t
use grid_field_mod,               only : grid_field_t, tile_field_t
use mesh_mod,                     only : tile_mesh_t
use halo_mod,                     only : halo_t

implicit none

type, extends(massflux_operator_t) :: massflux_Ch_t

    class(interpolator2d_vec2vec_t), allocatable :: interp_f2uv_op

    contains

    procedure :: calc_massflux => calc_ch_sbp_massflux

end type

type, extends(massflux_operator_t) :: massflux_Ch_up4_t

    class(halo_t), allocatable :: halo_h

    contains

    procedure :: calc_massflux => calc_ch_up4_massflux

end type massflux_Ch_up4_t

type, extends(massflux_operator_t) :: massflux_Ch_up5_t

    class(halo_t), allocatable :: halo_h

    contains

    procedure :: calc_massflux => calc_ch_up5_massflux

end type massflux_Ch_up5_t

contains

subroutine calc_ch_sbp_massflux(this, fx, fy, f, u, v, domain)

    class(massflux_Ch_t), intent(inout) :: this
    type(domain_t),       intent(in)    :: domain
    type(grid_field_t),   intent(inout) :: f, u, v
    !output:
    type(grid_field_t),   intent(inout) :: fx, fy

    integer(kind=4) :: t

    call this%interp_f2uv_op%interp2d_vec2vec(fx, fy, f, f, domain)

    call fx%assign_prod(1.0_8, fx, u, domain%mesh_y)
    call fy%assign_prod(1.0_8, fy, v, domain%mesh_x)

end subroutine

subroutine calc_ch_up4_massflux(this, fx, fy, f, u, v, domain)

    class(massflux_Ch_up4_t), intent(inout) :: this
    type(domain_t),           intent(in)    :: domain
    type(grid_field_t),       intent(inout) :: f, u, v
    !output:
    type(grid_field_t),       intent(inout) :: fx, fy

    integer(kind=4) :: t

    call this%halo_h%get_halo_scalar(f, domain, halo_width=3)

    do t = domain%mesh_p%ts, domain%mesh_p%te
        call calc_up4_massflux_tile(fx%tile(t), fy%tile(t), &
                                    f%tile(t), u%tile(t), v%tile(t), &
                                    domain%mesh_y%tile(t), domain%mesh_x%tile(t))
    end do

end subroutine

subroutine calc_up4_massflux_tile(fx,fy,f,u,v,mesh_u,mesh_v)

    type(tile_field_t), intent(inout) :: fx, fy
    type(tile_field_t), intent(in)    :: f, u, v
    type(tile_mesh_t),  intent(in)    :: mesh_u, mesh_v

    integer(kind=4) :: i, j, k, is, ie, js, je, ks, ke
    real(kind=8)    :: zl, zr

    ks = mesh_u%ks; ke = mesh_v%ke

    do k=ks,ke
        is = mesh_u%is; ie = mesh_u%ie
        js = mesh_u%js; je = mesh_u%je

        do j=js,je; do i=is,ie
            zl = .5_8+sign(.5_8,u%p(i,j,k))
            zr = 1._8-zl
            fx%p(i,j,k) = u%p(i,j,k)*( &
                          zl*(5._8*f%p(i+1,j,k)+15._8*f%p(i,j,k)-&
                                              5._8*f%p(i-1,j,k)+f%p(i-2,j,k))+&
                          zr*(5._8*f%p(i,j,k)+15._8*f%p(i+1,j,k)-&
                                              5._8*f%p(i+2,j,k)+f%p(i+3,j,k)))/16._8
        end do; end do

        is = mesh_v%is; ie = mesh_v%ie
        js = mesh_v%js; je = mesh_v%je

        do j=js,je; do i=is,ie
            zl = .5_8+sign(.5_8,v%p(i,j,k))
            zr = 1._8-zl
            fy%p(i,j,k) = v%p(i,j,k)*( &
                            zl*(5._8*f%p(i,j+1,k)+15._8*f%p(i,j,k)-&
                                              5._8*f%p(i,j-1,k)+f%p(i,j-2,k))+&
                            zr*(5._8*f%p(i,j,k)+15._8*f%p(i,j+1,k)-&
                                             5._8*f%p(i,j+2,k)+f%p(i,j+3,k)))/16._8
        end do; end do

    end do

end subroutine calc_up4_massflux_tile

subroutine calc_ch_up5_massflux(this, fx, fy, f, u, v, domain)

    class(massflux_Ch_up5_t), intent(inout) :: this
    type(domain_t),           intent(in)    :: domain
    type(grid_field_t),       intent(inout) :: f, u, v
    !output:
    type(grid_field_t),       intent(inout) :: fx, fy

    integer(kind=4) :: t

    call this%halo_h%get_halo_scalar(f, domain, halo_width=3)

    do t = domain%mesh_p%ts, domain%mesh_p%te
        call calc_up5_massflux_tile(fx%tile(t), fy%tile(t), &
                                    f%tile(t), u%tile(t), v%tile(t), &
                                    domain%mesh_y%tile(t), domain%mesh_x%tile(t))
    end do

end subroutine

subroutine calc_up5_massflux_tile(fx,fy,f,u,v,mesh_u,mesh_v)

    type(tile_field_t), intent(inout) :: fx, fy
    type(tile_field_t), intent(in)    :: f, u, v
    type(tile_mesh_t),  intent(in)    :: mesh_u, mesh_v

    integer(kind=4) :: i, j, k, is, ie, js, je, ks, ke
    real(kind=8)    :: zl, zr

    real(kind=8), parameter :: c1 = 3.0_8/128.0_8, c2 =-5.0_8/32.0_8,  &
                               c3 = 45.0_8/64.0_8, c4 = 15.0_8/32.0_8, &
                               c5 =-5.0_8/128.0_8

    ks = mesh_u%ks; ke = mesh_v%ke

    do k=ks,ke
        is = mesh_u%is; ie = mesh_u%ie
        js = mesh_u%js; je = mesh_u%je

        do j=js,je; do i=is,ie
            zl = .5_8+sign(.5_8,u%p(i,j,k))
            zr = 1._8-zl
            fx%p(i,j,k) = u%p(i,j,k)*( &
                          zl*(c5*f%p(i+2,j,k)+c4*f%p(i+1,j,k)+c3*f%p(i,j,k)+&
                              c2*f%p(i-1,j,k)+c1*f%p(i-2,j,k))+&
                          zr*(c5*f%p(i-1,j,k)+c4*f%p(i,j,k)+c3*f%p(i+1,j,k)+&
                              c2*f%p(i+2,j,k)+c1*f%p(i+3,j,k)))
        end do; end do

        is = mesh_v%is; ie = mesh_v%ie
        js = mesh_v%js; je = mesh_v%je

        do j=js,je; do i=is,ie
            zl = .5_8+sign(.5_8,v%p(i,j,k))
            zr = 1._8-zl
            fy%p(i,j,k) = v%p(i,j,k)*( &
                            zl*(c5*f%p(i,j+2,k)+c4*f%p(i,j+1,k)+c3*f%p(i,j,k)+&
                                c2*f%p(i,j-1,k)+c1*f%p(i,j-2,k))+&
                            zr*(c5*f%p(i,j-1,k)+c4*f%p(i,j,k)+c3*f%p(i,j+1,k)+&
                                c2*f%p(i,j+2,k)+c1*f%p(i,j+3,k)))
        end do; end do

    end do

end subroutine calc_up5_massflux_tile

end module massflux_ch_mod