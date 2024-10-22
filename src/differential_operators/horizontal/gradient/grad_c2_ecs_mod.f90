module grad_c2_ecs_mod

use abstract_grad_mod,      only : grad_operator_t
use grid_field_mod,         only : grid_field_t, tile_field_t
use mesh_mod,               only : mesh_t, tile_mesh_t
use halo_mod,               only : halo_t, halo_vec_t
use domain_mod,             only : domain_t
use exchange_abstract_mod,  only : exchange_t

implicit none

type, public, extends(grad_operator_t) :: grad_c2_ecs_t
    class(halo_t), allocatable :: halo_procedure
contains
    procedure, public :: calc_grad => calc_grad_c2_ecs
end type grad_c2_ecs_t

contains

subroutine calc_grad_c2_ecs(this, gx, gy, f, domain)
    class(grad_c2_ecs_t),  intent(inout) :: this
    type(grid_field_t),    intent(inout) :: gx
    type(grid_field_t),    intent(inout) :: gy
    type(grid_field_t),    intent(inout) :: f
    type(domain_t),        intent(in)    :: domain

    integer(kind=4) :: t
    integer(kind=4), parameter :: halo_width=1

    call this%halo_procedure%get_halo_scalar(f,domain,halo_width)

    do t = domain%partition%ts, domain%partition%te
        call calc_grad_on_tile(gx%tile(t), gy%tile(t), f%tile(t),            &
                               domain%mesh_x%tile(t), domain%mesh_y%tile(t), &
                               domain%mesh_o%tile(t),domain%mesh_o%scale)
    end do

end subroutine calc_grad_c2_ecs

subroutine calc_grad_on_tile(gx, gy, f, mesh_x, mesh_y, mesh_o, scale)

    type(tile_field_t),     intent(inout) :: gx, gy
    type(tile_field_t),     intent(in)    :: f
    type(tile_mesh_t),      intent(in)    :: mesh_x, mesh_y, mesh_o
    real(kind=8),           intent(in)    :: scale

    real(kind=8)    :: hx, mult_loc
    integer(kind=4) :: ks, ke
    integer(kind=4) :: jsu, jeu, isu, ieu
    integer(kind=4) :: jsv, jev, isv, iev
    integer(kind=4) :: jsx, jex, isx, iex
    integer(kind=4) :: jsy, jey, isy, iey
    integer(kind=4) :: i, j, k
    real(kind=8), allocatable :: fdx(:,:), fdy(:,:)
    real(kind=8)    :: fdx_at_y, fdy_at_x

    ks = mesh_o%ks; ke = mesh_o%ke

    isu = mesh_x%is; ieu = mesh_x%ie
    jsu = mesh_x%js; jeu = mesh_x%je
    isv = mesh_y%is; iev = mesh_y%ie
    jsv = mesh_y%js; jev = mesh_y%je

    isx = min(isu,isv);   iex = max(ieu,iev+1)
    jsx = min(jsu,jsv-1); jex = max(jeu,jev)
    allocate(fdx(isx:iex,jsx:jex))
    isy = min(isu-1,isv); iey = max(ieu,iev)
    jsy = min(jsu,jsv);   jey = max(jeu+1,jev)
    allocate(fdy(isy:iey,jsy:jey))

    hx = mesh_o%hx

    do k = ks, ke

        do j= jsx, jex
            do i= isx, iex
                gx%p(i,j,k) = (f%p(i,j,k)-f%p(i-1,j,k))/(hx*scale)
            end do
        end do

        do j= jsy, jey
            do i= isy, iey
                gy%p(i,j,k) = (f%p(i,j,k)-f%p(i,j-1,k))/(hx*scale)
            end do
        end do

    end do

end subroutine calc_grad_on_tile

end module grad_c2_ecs_mod
