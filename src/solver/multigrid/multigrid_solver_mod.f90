module multigrid_solver_mod

use domain_mod,           only : domain_t
use iterative_solver_mod, only : iterative_solver_t
use linear_operator_mod,  only : linear_operator_t
use vector_mod,           only : abstract_vector_t
use iterative_solver_mod, only : iterative_solver_t
use preconditioner_mod,   only : preconditioner_t

use abstract_restriction_mod,  only : restriction_operator_t
use abstract_prolongation_mod, only : prolongation_operator_t

implicit none

type, public, extends(iterative_solver_t) :: multigrid_solver_t
    integer(kind=4) :: n_levels
    integer(kind=4) :: n_pre_smooth, n_post_smooth, n_coarse_smooth
    type(mg_data_t), allocatable :: mg_data(:)
contains
    procedure, public  :: solve
    procedure, private :: V_cycle
    procedure, private :: smooth
end type multigrid_solver_t

type, public :: mg_data_t
    type(domain_t)                              :: domain
    class(abstract_vector_t),       allocatable :: x, r, b, r_prec
    class(linear_operator_t),       allocatable :: linear_operator
    class(restriction_operator_t),  allocatable :: restriction_op
    class(prolongation_operator_t), allocatable :: prolongation_op
    class(preconditioner_t),        allocatable :: smooth_precond
    class(iterative_solver_t),      allocatable :: coarse_solver
    class(iterative_solver_t),      allocatable :: smoother
contains
end type mg_data_t

contains

subroutine solve(this, x, linear_operator, rhs, domain, precond)

    class(multigrid_solver_t), intent(inout) :: this
    class(linear_operator_t),  intent(inout) :: linear_operator
    class(abstract_vector_t),  intent(inout) :: x, rhs
    type(domain_t),            intent(in)    :: domain
    class(preconditioner_t),   intent(inout), optional :: precond

    integer(kind=4) :: iter
    real(kind=8) :: r_norm, rhs_norm
    logical :: is_converged

    is_converged = .true.

    call this%mg_data(1)%x%copy(x,   domain)
    call this%mg_data(1)%b%copy(rhs, domain)

    call linear_operator%apply(this%mg_data(1)%r, this%mg_data(1)%x, this%mg_data(1)%domain) ! r = A*x
    call this%mg_data(1)%r%assign(-1.0_8, this%mg_data(1)%r, 1.0_8, rhs, domain) ! r = b - A*x

    r_norm   =   this%mg_data(1)%r%norm(domain)
    rhs_norm = rhs%norm(domain)

    if (r_norm / rhs_norm < this%relative_tolerance) then
        is_converged = .true.
        if (this%is_verbose) call domain%parcomm%print("Multigrif converged in 0 iterations")
        return
    end if

    do iter = 1, this%max_iter

        call this%V_cycle()
        call linear_operator%apply(this%mg_data(1)%r, this%mg_data(1)%x, this%mg_data(1)%domain) ! r = A*x
        call this%mg_data(1)%r%assign(-1.0_8, this%mg_data(1)%r, 1.0_8, rhs, domain) ! r = b - A*x

        r_norm   =   this%mg_data(1)%r%norm(domain)
        rhs_norm = rhs%norm(domain)

        if (domain%parcomm%myid == 0) then
                if (this%is_verbose ) print*, iter, r_norm / rhs_norm
        end if

        if ( r_norm / rhs_norm < this%relative_tolerance ) then
            is_converged = .true.
            exit
        end if
    end do

    call x%copy(this%mg_data(1)%x, domain)

    if (domain%parcomm%myid == 0) then
        if (is_converged) then
            if (this%is_verbose ) print*, "Multigrid converged in ", iter, " iterations"
        else
            print*, "Multigrid failed to converge in ", this%max_iter, " iterations"
        end if
    end if

end subroutine solve

subroutine V_cycle(this)

    class(multigrid_solver_t), intent(inout) :: this

    integer(kind=4) :: l

    do l = 1, this%n_levels - 1

        this%mg_data(l)%smoother%max_iter = this%n_pre_smooth

        if (l /= 1) call this%mg_data(l)%x%set_scalar(0.0_8, this%mg_data(l)%domain)

        call this%mg_data(l)%smoother%solve(this%mg_data(l)%x, &
                                         this%mg_data(l)%linear_operator, &
                                         this%mg_data(l)%b, &
                                         this%mg_data(l)%domain, &
                                         this%mg_data(l)%smooth_precond)

        ! call this%smooth(this%n_pre_smooth, 4.0_8/5, l, l /= 1)

        call this%mg_data(l)%linear_operator%apply(&
                this%mg_data(l)%r, this%mg_data(l)%x, this%mg_data(l)%domain)

        call this%mg_data(l)%r%assign(&
                1.0_8, this%mg_data(l)%b, -1.0_8, this%mg_data(l)%r, this%mg_data(l)%domain)

        call this%mg_data(l)%restriction_op%apply(&
               this%mg_data(l+1)%b, this%mg_data(l)%r, this%mg_data(l+1)%domain, this%mg_data(l)%domain)

    end do

    l = this%n_levels

    call this%mg_data(l)%x%set_scalar(0.0_8, this%mg_data(l)%domain)

    call this%mg_data(l)%coarse_solver%solve(this%mg_data(l)%x, this%mg_data(l)%linear_operator, this%mg_data(l)%b, this%mg_data(l)%domain)

    ! call this%smooth(this%n_coarse_smooth, 1.0_8, l, this%n_levels /= 1)

    do l = this%n_levels - 1, 1, -1

        call this%mg_data(l)%prolongation_op%apply(&
                  this%mg_data(l)%r, this%mg_data(l+1)%x, this%mg_data(l)%domain, this%mg_data(l+1)%domain)

        call this%mg_data(l)%x%update(1.0_8, this%mg_data(l)%r, this%mg_data(l)%domain)

        this%mg_data(l)%smoother%max_iter = this%n_post_smooth

        call this%mg_data(l)%smoother%solve(this%mg_data(l)%x, &
                                         this%mg_data(l)%linear_operator, &
                                         this%mg_data(l)%b, &
                                         this%mg_data(l)%domain, &
                                         this%mg_data(l)%smooth_precond)

        ! call this%smooth(this%n_post_smooth, 4.0_8/5, l, .false.)

    end do

end subroutine V_cycle

subroutine smooth(this, n_smooth, w, level, is_zero_initial)

    class(multigrid_solver_t), intent(inout) :: this
    integer(kind=4),           intent(in)    :: n_smooth, level
    real(kind=8),              intent(in)    :: w
    logical,                   intent(in)    :: is_zero_initial

    integer(kind=4) :: n, iter

    iter = 1
    if (is_zero_initial) then
        call this%mg_data(level)%smooth_precond%apply( this%mg_data(level)%r_prec, &
                                            this%mg_data(level)%b,      &
                                            this%mg_data(level)%domain )
        call this%mg_data(level)%x%assign(w, this%mg_data(level)%r_prec, this%mg_data(level)%domain)
        iter = 2
    end if

    do n = iter, n_smooth

        call this%mg_data(level)%linear_operator%apply( this%mg_data(level)%r, &
                                                        this%mg_data(level)%x, &
                                                        this%mg_data(level)%domain )

        call this%mg_data(level)%r%assign( 1.0_8, this%mg_data(level)%b, &
                                  -1.0_8, this%mg_data(level)%r, this%mg_data(level)%domain)

        call this%mg_data(level)%smooth_precond%apply( this%mg_data(level)%r_prec, &
                                            this%mg_data(level)%r,      &
                                            this%mg_data(level)%domain )

        call this%mg_data(level)%x%update(w, this%mg_data(level)%r_prec, this%mg_data(level)%domain)

    end do

end subroutine smooth



end module multigrid_solver_mod
