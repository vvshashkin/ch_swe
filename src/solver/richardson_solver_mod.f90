module richardson_solver_mod

use domain_mod,           only : domain_t
use iterative_solver_mod, only : iterative_solver_t
use linear_operator_mod,  only : linear_operator_t
use vector_mod,           only : abstract_vector_t
use preconditioner_mod,   only : preconditioner_t

implicit none

type, public, extends(iterative_solver_t) :: richardson_solver_t
    real(kind=8) :: omega
contains
    procedure, public :: solve
end type richardson_solver_t

contains

subroutine solve(this, x, linear_operator, rhs, domain, precond)
    class(richardson_solver_t), intent(inout) :: this
    class(linear_operator_t),   intent(inout) :: linear_operator
    class(abstract_vector_t),   intent(inout) :: x, rhs
    type(domain_t),             intent(in)    :: domain
    class(preconditioner_t),    intent(inout), optional :: precond

    class(abstract_vector_t), allocatable :: r, r_prec

    integer(kind=4) :: iter

    call x%create_similar(r, domain)

    if (present(precond)) call x%create_similar(r_prec, domain)

    do iter = 1, this%max_iter

        call linear_operator%apply(r, x, domain)
        call r%assign(-1.0_8, r, 1.0_8, rhs, domain)

        if (present(precond)) then
            call precond%apply(r_prec, r, domain)
            call x%update(this%omega, r_prec, domain)
        else
            call x%update(this%omega, r, domain)
        end if

    end do

end subroutine solve

end module richardson_solver_mod
