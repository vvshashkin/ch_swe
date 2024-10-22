module chebyshev_solver_mod

use domain_mod,           only : domain_t
use iterative_solver_mod, only : iterative_solver_t
use linear_operator_mod,  only : linear_operator_t
use vector_mod,           only : abstract_vector_t
use preconditioner_mod,   only : preconditioner_t

implicit none

type, public, extends(iterative_solver_t) :: chebyshev_solver_t
    real(kind=8) :: lambda_min, lambda_max
contains
    procedure, public :: solve
end type chebyshev_solver_t

contains


subroutine solve(this, x, linear_operator, rhs, domain, precond)

    class(chebyshev_solver_t), intent(inout) :: this
    class(linear_operator_t),  intent(inout) :: linear_operator
    class(abstract_vector_t),  intent(inout) :: x, rhs
    type(domain_t),            intent(in)    :: domain
    !ATTENTION this is very preliminary approach to preconditionning
    !please feel free to rewrite it as soon as alternative conception is
    class(preconditioner_t),   intent(inout), optional :: precond

    class(abstract_vector_t), allocatable :: r, d

    real(kind=8)    :: rho, rho_old, sigma_1, delta, theta, rhs_norm, r_norm
    integer(kind=4) :: iter
    logical         :: is_converged

    is_converged = .false.

    call x%create_similar(r, domain)
    call x%create_similar(d, domain)

    theta = (this%lambda_max + this%lambda_min) / 2.0_8
    delta = (this%lambda_max - this%lambda_min) / 2.0_8

    sigma_1 = theta / delta
    rho_old = 1.0_8 / sigma_1

    call linear_operator%apply(r, x, domain)     ! r = A*x
    call r%assign(-1.0_8, r, 1.0_8, rhs, domain) ! r = b - A*x

    call d%assign(1.0_8 / theta, r, domain)

    rhs_norm = rhs%norm(domain)

    do iter = 1, this%max_iter

        call x%update(1.0_8, d, domain)

        call linear_operator%apply(r, x, domain)     ! r = A*x
        call r%assign(-1.0_8, r, 1.0_8, rhs, domain) ! r = b - A*x

        rho = 1.0_8 / (2.0_8*sigma_1 - rho_old)
        call d%assign(rho*rho_old, d, 2.0_8*rho/delta, r, domain)

        rho_old = rho
        r_norm = r%norm(domain)

        if (domain%parcomm%myid == 0) then
            print*, "iter ", iter, r_norm / rhs_norm
        end if

        if ( r_norm / rhs_norm < this%relative_tolerance ) then
            is_converged = .true.
            exit
        end if

    end do

    if (domain%parcomm%myid == 0) then
        if (is_converged) then
            if (this%is_verbose ) print*, "Chebyshev solver converged in ", iter, " iterations"
        else
            print*, "Chebyshev solver failed to converge in ", this%max_iter, " iterations"
        end if
    end if

end subroutine solve

end module chebyshev_solver_mod
