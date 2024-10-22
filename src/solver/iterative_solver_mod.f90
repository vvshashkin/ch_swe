module iterative_solver_mod

use domain_mod,          only : domain_t
use linear_operator_mod, only : linear_operator_t
use vector_mod,          only : abstract_vector_t
use preconditioner_mod,  only : preconditioner_t

implicit none

private

! Base abstract class for iterative solver
! Defines an abstract class for iterative solution
! of the the system Ax=b.
type, public, abstract :: iterative_solver_t
    real(kind=8)    :: relative_tolerance
    integer(kind=4) :: max_iter
    logical         :: is_verbose ! print convergence info flag 
contains
    procedure(solve_i), deferred :: solve
end type iterative_solver_t

abstract interface
    subroutine solve_i(this, x, linear_operator, rhs, domain, precond)
        import :: iterative_solver_t, linear_operator_t, abstract_vector_t, domain_t, preconditioner_t
        class(iterative_solver_t), intent(inout) :: this
        class(linear_operator_t),  intent(inout) :: linear_operator
        class(abstract_vector_t),  intent(inout) :: x, rhs
        type(domain_t),            intent(in)    :: domain
        !ATTENTION this is very preliminary approach to preconditionning
        !please feel free to rewrite it as soon as alternative conception is ready
        class(preconditioner_t),   intent(inout), optional :: precond
    end subroutine solve_i
end interface

contains

end module iterative_solver_mod
