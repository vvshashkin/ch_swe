module linear_operator_mod

use domain_mod, only : domain_t
use vector_mod, only : abstract_vector_t

implicit none

private

! Base abstract type for linear operator
! defines an interface for the linear operator application y = A*x
type, public, abstract :: linear_operator_t
    private
contains
    procedure(apply_i), deferred :: apply
end type linear_operator_t

abstract interface
    subroutine apply_i(this, vout, vin, domain)
        import :: linear_operator_t, domain_t, abstract_vector_t
        class(linear_operator_t), intent(inout) :: this
        class(abstract_vector_t), intent(inout) :: vout, vin
        type(domain_t),           intent(in)    :: domain
    end subroutine apply_i
end interface

contains

end module linear_operator_mod
