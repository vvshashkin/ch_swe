module div_c2_mod

use abstract_div_mod,   only : div_operator_t
use domain_mod,         only : domain_t
use mesh_mod,           only : mesh_t, tile_mesh_t
use grid_field_mod,     only : grid_field_t, tile_field_t
use exchange_halo_mod,  only : exchange_t
use halo_mod,           only : halo_vec_t

implicit none

type, public, extends(div_operator_t) :: div_c2_t
    class(halo_vec_t), allocatable :: halo_procedure
contains
    procedure, public :: calc_div => calc_div_c2
end type div_c2_t

contains

subroutine calc_div_c2(this, div, u, v, domain)

    class(div_c2_t),        intent(inout) :: this
    type(domain_t),         intent(in)    :: domain
    type(grid_field_t),     intent(inout) :: u, v
    !out put
    type(grid_field_t),     intent(inout) :: div

    integer(kind=4) :: i, j, k, t
    real(kind=8) hx

    call this%halo_procedure%get_halo_vector(u,v,domain,1)

    do t = domain%partition%ts, domain%partition%te
        call calc_div_on_tile(div%tile(t), u%tile(t), v%tile(t),            &
                              domain%mesh_x%tile(t), domain%mesh_y%tile(t), &
                              domain%mesh_o%tile(t),domain%mesh_o%scale)
    end do

end subroutine calc_div_c2

subroutine calc_div_on_tile(div, u, v, mesh_u, mesh_v, mesh_p, scale)

    type(tile_field_t),  intent(inout) :: div
    type(tile_field_t),  intent(in)    :: u, v
    type(tile_mesh_t),   intent(in)    :: mesh_u, mesh_v, mesh_p
    real(kind=8),        intent(in)    :: scale

    real(kind=8) :: hx
    integer(kind=4) :: ks, ke, js, je, is, ie, i, j, k

    is = mesh_p%is; ie = mesh_p%ie
    js = mesh_p%js; je = mesh_p%je
    ks = mesh_p%ks; ke = mesh_p%ke

    hx = mesh_p%hx

    do k = ks, ke
        do j = js, je
            do i = is, ie
                div%p(i,j,k) = (mesh_u%J(i+1,j,k)*u%p(i+1,j,k)-mesh_u%J(i,j,k)*u%p(i,j,k)  +  &
                                mesh_v%J(i,j+1,k)*v%p(i,j+1,k)-mesh_v%J(i,j,k)*v%p(i,j,k))/  &
                                (mesh_p%J(i,j,k)*hx*scale)
            end do
        end do
    end do

end subroutine calc_div_on_tile

end module div_c2_mod
