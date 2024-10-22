module abstract_curl_mod

use grid_field_mod,    only : grid_field_t
use domain_mod,        only : domain_t

implicit none

type, abstract, public :: curl_operator_t
contains
    procedure(calc_curl_i), deferred :: calc_curl
end type curl_operator_t

abstract interface
    subroutine calc_curl_i(this, curl, u, v, domain)
        import curl_operator_t, grid_field_t, domain_t
        class(curl_operator_t),  intent(inout) :: this
        type(domain_t),          intent(in)    :: domain
        type(grid_field_t),      intent(inout) :: u, v
        type(grid_field_t),      intent(inout) :: curl
    end subroutine calc_curl_i
end interface

contains

end module abstract_curl_mod
