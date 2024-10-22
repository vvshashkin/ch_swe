module bicgstab_solver_mod

use domain_mod,           only : domain_t
use iterative_solver_mod, only : iterative_solver_t
use linear_operator_mod,  only : linear_operator_t
use vector_mod,           only : abstract_vector_t
use preconditioner_mod,  only : preconditioner_t

implicit none

private

type, public, extends(iterative_solver_t) :: bicgstab_solver_t
private
contains
    procedure, public :: solve
end type bicgstab_solver_t

contains

subroutine solve(this, x, linear_operator, rhs, domain, precond)

    class(bicgstab_solver_t), intent(inout) :: this
    class(linear_operator_t), intent(inout) :: linear_operator
    class(abstract_vector_t), intent(inout) :: x, rhs
    type(domain_t),           intent(in)    :: domain
    !ATTENTION this is very preliminary approach to preconditionning
    !please feel free to rewrite it as soon as alternative conception is
    class(preconditioner_t),   intent(inout), optional :: precond

    class(abstract_vector_t), allocatable :: r, v, p, r0, s, t, pre_buff, pre_rhs

    real(kind=8)    :: rhs_norm, r_norm, r_norm_new
    real(kind=8)    :: alpha, beta, rho, rho_old, omega
    integer(kind=4) :: iter
    logical         :: is_converged

    is_converged = .false.

    if(present(precond)) then
        call x%create_similar(pre_buff, domain)
        call x%create_similar(pre_rhs, domain)
        call precond%apply(pre_rhs,rhs,domain)
    end if

    call x%create_similar(r,  domain)
    if(present(precond)) then
        call linear_operator%apply(pre_buff, x, domain)     ! r = B*A*x
        call precond%apply(r,pre_buff,domain)
        call r%assign(-1.0_8, r, 1.0_8, pre_rhs, domain) ! r = B*b - B*A*x
    else
        call linear_operator%apply(r, x, domain)     ! r = A*x
        call r%assign(-1.0_8, r, 1.0_8, rhs, domain) ! r = b - A*x
    end if

    r_norm   =   r%norm(domain)
    if(present(precond)) then
        rhs_norm = pre_rhs%norm(domain)
    else
        rhs_norm = rhs%norm(domain)
    end if

    if (r_norm / rhs_norm < this%relative_tolerance) then
        is_converged = .true.
        if (this%is_verbose) call domain%parcomm%print("bicgstab converged in 0 iterations")
        return
    end if

    alpha   = 1.0_8
    omega   = 1.0_8
    rho_old = 1.0_8

    call x%create_similar(v, domain)
    call x%create_similar(p, domain)
    call x%create_similar(s, domain)
    call x%create_similar(t, domain)

    call v%set_scalar(0.0_8, domain)  ! v = 0
    call p%set_scalar(0.0_8, domain)  ! p = 0

    call x%create_similar(r0, domain)
    call r0%copy(r, domain)           ! r0 = r

    do iter = 1, this%max_iter

        rho = r%dot(r0, domain) ! rho = (r,r0)

        beta = (rho / rho_old) * (alpha / omega)

        call p%scale(beta, domain)            ! p = beta*p
        call p%update(1.0_8,       r, domain) ! p = p + r
        call p%update(-beta*omega, v, domain) ! p = p - beta*omega*v

        if(present(precond)) then
            call linear_operator%apply(pre_buff, p, domain) ! v = A*p
            call precond%apply(v,pre_buff,domain) ! v = B*v
        else
            call linear_operator%apply(v, p, domain) ! v = A*p
        end if

        alpha = rho / v%dot(r0, domain)

        call s%assign(1.0_8, r, -alpha, v, domain) ! s = r -alpha*v

        if(present(precond)) then
            call linear_operator%apply(pre_buff, s, domain) ! t = A*s
            call precond%apply(t,pre_buff,domain) ! s = B*s
        else
            call linear_operator%apply(t, s, domain) !t = A*s
        end if

        omega = t%dot(s, domain) / (t%norm(domain))**2 ! omega = (t,s)/(t,t)

        call x%update( alpha, p,  domain) ! x = x + alpha*p
        call x%update( omega, s,  domain) ! x = x + omega*s

        call r%assign(1.0_8, s, -omega, t, domain) ! r = s - omega*t

        rho_old = rho

        r_norm = r%norm(domain)

        if ( r_norm / rhs_norm < this%relative_tolerance ) then
            is_converged = .true.
            exit
        end if

    end do

    if (domain%parcomm%myid == 0) then
        if (is_converged) then
            if (this%is_verbose )  print*, "bicgstab converged in ", iter, " iterations"
        else
            print*, "bicgstab failed to converge in ", this%max_iter, " iterations, rel_tol=", r_norm/rhs_norm
        end if
    end if

end subroutine solve


end module bicgstab_solver_mod
