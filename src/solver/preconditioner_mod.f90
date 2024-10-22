module preconditioner_mod

use domain_mod, only : domain_t
use vector_mod, only : abstract_vector_t

implicit none

private

! Base abstract class for preconditioner
! Defines an abstract class for the application of preconditioner y=P^{-1}x

type, public, abstract :: preconditioner_t
private
contains
    procedure(apply_preconditioner_i), public, deferred :: apply
end type preconditioner_t

abstract interface
    subroutine apply_preconditioner_i(this, y, x, domain)
        import :: preconditioner_t, abstract_vector_t, domain_t
        class(preconditioner_t),  intent(inout) :: this
        class(abstract_vector_t), intent(inout) :: x, y
        type(domain_t),           intent(in)    :: domain
    end subroutine apply_preconditioner_i
end interface

contains

end module preconditioner_mod
