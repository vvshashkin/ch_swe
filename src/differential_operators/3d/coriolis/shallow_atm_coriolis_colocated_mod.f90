module shallow_atm_coriolis_colocated_mod

use abstract_coriolis_3d_mod, only : coriolis3d_operator_t
use grid_field_mod,           only : grid_field_t, tile_field_t
use mesh_mod,                 only : tile_mesh_t
use domain_mod,               only : domain_t
use parcomm_mod,              only : parcomm_global

implicit none

type, public, extends(coriolis3d_operator_t) :: shallow_atm_coriolis_colocated_t
    type(grid_field_t) :: u2u, u2v, v2u, v2v
contains
    procedure :: calc_coriolis
end type shallow_atm_coriolis_colocated_t

contains

subroutine calc_coriolis(this, u_tend, v_tend, w_tend, u, v, w, domain)
    class(shallow_atm_coriolis_colocated_t), intent(inout) :: this
    type(domain_t),                          intent(in)    :: domain
    type(grid_field_t),                      intent(inout) :: u, v, w
    type(grid_field_t),                      intent(inout) :: u_tend, v_tend, w_tend

    integer(kind=4) :: t

    call w_tend%assign(0.0_8,domain%mesh_w)

    do t = domain%mesh_u%ts, domain%mesh_v%te
        call calc_coriolis_tile(u_tend%tile(t), v_tend%tile(t),     &
                                u%tile(t), v%tile(t),               &
                                this%u2u%tile(t), this%u2v%tile(t), &
                                this%v2u%tile(t), this%v2v%tile(t), &
                                domain%mesh_u%tile(t))
    end do

end subroutine calc_coriolis

subroutine calc_coriolis_tile(u_tend, v_tend, u, v, u2u, u2v, v2u, v2v, mesh)
    type(tile_field_t), intent(inout) :: u_tend, v_tend
    type(tile_field_t), intent(in)    :: u, v, u2u, u2v, v2u, v2v
    type(tile_mesh_t),  intent(in)    :: mesh

    integer(kind=4) :: i, j, k

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                u_tend%p(i,j,k) = u2u%p(i,j,k)*u%p(i,j,k)+v2u%p(i,j,k)*v%p(i,j,k)
                v_tend%p(i,j,k) = u2v%p(i,j,k)*u%p(i,j,k)+v2v%p(i,j,k)*v%p(i,j,k)
            end do
        end do
    end do
end subroutine

end module shallow_atm_coriolis_colocated_mod
