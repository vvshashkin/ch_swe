module imex_M2_mod

use timescheme_mod, only : timescheme_t
use operator_mod,   only : operator_t, imex_operator_t
use stvec_mod,      only : stvec_t
use domain_mod,     only : domain_t
use parcomm_mod,    only : parcomm_global

implicit none

real(kind=8), parameter :: imex_M2_c(5) = [1.0_8/4.0_8, 1.0_8/6.0_8, 3.0_8/8.0_8, 1.0_8/2.0_8, 1.0_8]

type, extends(timescheme_t), public :: imex_M2_t

    real(kind=8) :: c(5), d(6)
    class(stvec_t),    allocatable :: ye, yi, res
    class(stvec_t),    allocatable :: r,s

contains

    procedure, public :: step => step_imex_M2

end type imex_M2_t

contains

subroutine step_imex_M2(this, v0, operator, domain, dt)

    class(imex_M2_t),       intent(inout) :: this
    class(stvec_t),         intent(inout) :: v0
    class(operator_t),      intent(inout) :: operator
    type(domain_t),         intent(in)    :: domain
    real(kind=8),           intent(in)    :: dt

    integer(kind=4) :: i

    select type(operator)
    class is(imex_operator_t)

    !first stage:
    call this%res%assign(v0, domain)
    call operator%apply_explicit(this%ye,v0,domain)
    if(this%d(1) /= 0.0_8) then
        call operator%apply_implicit(this%yi,v0,domain)
        call this%res%update(this%d(1),this%yi,domain)
    end if

    !inner stages:
    do i = 1, 4
        call this%r%assign(1.0_8,v0,this%c(i)*dt,this%ye,domain)
        call operator%solve_implicit(this%s,this%r,this%c(i)*dt,domain)
        call operator%apply_explicit(this%ye,this%s,domain)
        call this%yi%assign(1.0_8/(this%c(i)*dt),this%s,-1.0_8/(this%c(i)*dt),this%r,domain)
        if(this%d(i+1) /= 0.0_8) &
            call this%res%update(dt*this%d(i+1),this%yi,domain)
    end do

    !last stage
    call this%res%update(this%c(5)*dt,this%ye,domain)
    call operator%solve_implicit(v0,this%res,this%d(6)*dt,domain)

    class default !default operator
        call parcomm_global%abort("imex_M2 can work only with imex operators")
    end select

end subroutine step_imex_M2

end module imex_M2_mod
