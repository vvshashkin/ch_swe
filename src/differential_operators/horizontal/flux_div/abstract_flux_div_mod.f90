module abstract_flux_div_mod

use grid_field_mod,    only : grid_field_t
use domain_mod,        only : domain_t

implicit none

type, abstract :: flux_div_operator_t
    contains
        procedure(calc_flux_div_i), deferred :: calc_flux_div
end type

abstract interface
    subroutine calc_flux_div_i(this, flux_div, h, u, v, domain)

        import grid_field_t, domain_t, flux_div_operator_t

        !input:
        class(flux_div_operator_t), intent(inout) :: this
        type(grid_field_t),         intent(inout) :: h, u, v
        type(domain_t),             intent(in)    :: domain

        !output:
        type(grid_field_t),         intent(inout) :: flux_div

    end subroutine
end interface

end module abstract_flux_div_mod