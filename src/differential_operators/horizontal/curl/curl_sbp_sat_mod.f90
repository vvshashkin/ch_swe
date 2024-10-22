module curl_sbp_sat_mod

use grid_field_mod,        only : grid_field_t, tile_field_t
use domain_mod,            only : domain_t
use mesh_mod,              only : tile_mesh_t
use abstract_curl_mod,     only : curl_operator_t
use sbp_diff_mod,          only : sbp_diff_t
use exchange_abstract_mod, only : exchange_t
use mesh_mod,              only : mesh_t
use halo_mod,              only : halo_t

implicit none

type, public, extends(curl_operator_t) :: curl_sbp_sat_t
    type(grid_field_t)             :: Dx, Dy
    class(exchange_t), allocatable :: exch_halo
    class(sbp_diff_t), allocatable :: sbp_diff
    class(halo_t),     allocatable :: sync_edges
    type(mesh_t),      pointer     :: mesh_w
contains
    procedure, public :: calc_curl
    procedure, public :: calc_curl_tile
end type curl_sbp_sat_t

contains

subroutine calc_curl(this, curl, u, v, domain)

    class(curl_sbp_sat_t), intent(inout) :: this
    type(domain_t),        intent(in)    :: domain
    type(grid_field_t),    intent(inout) :: u, v !covariant components
    type(grid_field_t),    intent(inout) :: curl

    integer(kind=4) :: t

    call this%exch_halo%do_vec(u, v, domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call this%calc_curl_tile(curl%tile(t), u%tile(t), v%tile(t),   &
                                 this%Dx%tile(t), this%Dy%tile(t),    &
                                 this%mesh_w%tile(t), this%mesh_w%scale)
    end do

    if(allocated(this%sync_edges)) call this%sync_edges%get_halo_scalar(curl, domain, 1)

end subroutine calc_curl

subroutine calc_curl_tile(this, curl, u, v, Dx, Dy, mesh_w, scale)

    class(curl_sbp_sat_t), intent(inout) :: this
    type(tile_field_t),    intent(inout) :: curl
    type(tile_field_t),    intent(in)    :: u, v
    type(tile_field_t),    intent(inout) :: Dx, Dy
    type(tile_mesh_t),     intent(in)    :: mesh_w
    real(kind=8),          intent(in)    :: scale

    integer(kind=4) :: i, j, k

    call this%sbp_diff%apply_tile(Dy, u, mesh_w, "y")
    call this%sbp_diff%add_SAT_correction(Dy, u, mesh_w, "y")

    call this%sbp_diff%apply_tile(Dx, v, mesh_w, "x")
    call this%sbp_diff%add_SAT_correction(Dx, v, mesh_w, "x")

    do k = mesh_w%ks, mesh_w%ke
        do j = mesh_w%js, mesh_w%je
            do i = mesh_w%is, mesh_w%ie
                curl%p(i,j,k) = (Dx%p(i,j,k) / mesh_w%hx - Dy%p(i,j,k)/ mesh_w%hy) &
                             / (mesh_w%J(i,j,k)*scale)
            end do
        end do
    end do

end subroutine calc_curl_tile

end module curl_sbp_sat_mod
