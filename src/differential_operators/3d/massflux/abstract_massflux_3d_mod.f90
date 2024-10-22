module abstract_massflux_3d_mod

use grid_field_mod, only : grid_field_t
use domain_mod,     only : domain_t

implicit none

type, abstract, public :: massflux_3d_operator_t
    contains
        procedure(calc_massflux), public, deferred :: calc_massflux
end type

abstract interface

    subroutine calc_massflux(this, fx, fy, fz, f, u, v, w, domain)
        import massflux_3d_operator_t, grid_field_t, domain_t
        class(massflux_3d_operator_t), intent(inout) :: this
        type(domain_t),                intent(in)    :: domain
        type(grid_field_t),            intent(inout) :: f, u, v, w
        !output:
        type(grid_field_t),            intent(inout) :: fx, fy, fz
    end subroutine

end interface

end module
