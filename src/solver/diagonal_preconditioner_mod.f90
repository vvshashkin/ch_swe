module diagonal_preconditioner_mod

use domain_mod,         only : domain_t
use preconditioner_mod, only : preconditioner_t
use vector_mod,         only : abstract_vector_t

implicit none

type, public, extends(preconditioner_t) :: diagonal_preconditioner_t
    class(abstract_vector_t), allocatable :: diag
contains
    procedure, public :: apply
end type diagonal_preconditioner_t

contains

subroutine apply(this, y, x, domain)

    class(diagonal_preconditioner_t), intent(inout) :: this
    class(abstract_vector_t),         intent(inout) :: x, y
    type(domain_t),                   intent(in)    :: domain

    call y%assign_prod(this%diag, x, domain)

end subroutine apply

end module diagonal_preconditioner_mod
