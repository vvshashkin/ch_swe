module abstract_coriolis_3d_mod

use grid_field_mod,    only : grid_field_t
use domain_mod,        only : domain_t
use parcomm_mod,       only : parcomm_global

implicit none

type, abstract, public :: coriolis3d_operator_t
contains
    procedure(calc_coriolis_i),         deferred :: calc_coriolis
end type coriolis3d_operator_t

abstract interface
    subroutine calc_coriolis_i(this, u_tend, v_tend, w_tend, u, v, w, domain)
        import coriolis3d_operator_t, grid_field_t, domain_t
        class(coriolis3d_operator_t), intent(inout) :: this
        type(domain_t),               intent(in)    :: domain
        type(grid_field_t),           intent(inout) :: u, v, w
        type(grid_field_t),           intent(inout) :: u_tend, v_tend, w_tend
    end subroutine calc_coriolis_i
end interface

contains

end module abstract_coriolis_3d_mod
