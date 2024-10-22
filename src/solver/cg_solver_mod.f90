module cg_solver_mod

use domain_mod,           only : domain_t
use iterative_solver_mod, only : iterative_solver_t
use linear_operator_mod,  only : linear_operator_t
use vector_mod,           only : abstract_vector_t
use preconditioner_mod,   only : preconditioner_t

implicit none

private

type, public, extends(iterative_solver_t) :: cg_solver_t
private
contains
    procedure, public :: solve
end type cg_solver_t

contains

subroutine solve(this, x, linear_operator, rhs, domain, precond)

    class(cg_solver_t),        intent(inout) :: this
    class(linear_operator_t),  intent(inout) :: linear_operator
    class(abstract_vector_t),  intent(inout) :: x, rhs
    type(domain_t),            intent(in)    :: domain
    !ATTENTION this is very preliminary approach to preconditionning
    !please feel free to rewrite it as soon as alternative conception is
    class(preconditioner_t),   intent(inout), optional :: precond

    class(abstract_vector_t), allocatable :: r, p, Ap

    real(kind=8)       :: rhs_norm, r_norm, r_norm_new, alpha, beta
    integer(kind=4)    :: iter
    logical            :: is_converged

    !make abort to warn user that precond takes no effect here
    if(present(precond)) call domain%parcomm%abort("CG is currently unable to use preconditioner")

    is_converged = .false.

    call x%create_similar(r,  domain)
    call x%create_similar(p,  domain)
    call x%create_similar(Ap, domain)

    call linear_operator%apply(r, x, domain)     ! r = A*x
    call r%assign(-1.0_8, r, 1.0_8, rhs, domain) ! r = b - A*x

    r_norm   =   r%norm(domain)
    rhs_norm = rhs%norm(domain)

    if (r_norm / rhs_norm < this%relative_tolerance) then
        is_converged = .true.
        if (this%is_verbose) call domain%parcomm%print("CG converged in 0 iterations")
        return
    end if

    call p%copy(r, domain) ! r = p

    do iter = 1, this%max_iter

        call linear_operator%apply(Ap, p, domain) ! Ap = A*p

        alpha = r_norm**2 / p%dot(Ap, domain)     ! alpha = (r,r)/(Ap,p)

        call x%update( alpha, p,  domain)         ! x = x + alpha*p
        call r%update(-alpha, Ap, domain)         ! r = r - alpha*Ap

        r_norm_new = r%norm(domain)

        if ( r_norm_new / rhs_norm < this%relative_tolerance ) then
            is_converged = .true.
            exit
        end if

        beta = (r_norm_new / r_norm)**2
        r_norm = r_norm_new

        call p%assign(1.0_8, r, beta, p, domain) ! p = r + beta*p

    end do

    if (domain%parcomm%myid == 0) then
        if (is_converged) then
            if (this%is_verbose ) print*, "CG converged in ", iter, " iterations"
        else
            print*, "CG failed to converge in ", this%max_iter, " iterations"
        end if
    end if

end subroutine solve

end module cg_solver_mod
