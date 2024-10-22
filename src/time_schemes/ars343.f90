module ars343_mod

use timescheme_mod, only : timescheme_t
use operator_mod,   only : operator_t, imex_operator_t
use stvec_mod,      only : stvec_t
use domain_mod,     only : domain_t
use parcomm_mod,    only : parcomm_global

implicit none

real(kind=8), parameter :: ae(4,4) = [ [0.0_8         , 0.0_8         , 0.0_8         , 0.0_8], &
                                       [0.4358665215_8, 0.0_8         , 0.0_8         , 0.0_8], &
                                       [0.3212788860_8, 0.3966543747_8, 0.0_8         , 0.0_8], &
                                       [-0.105858296_8, 0.5529291479_8, 0.5529291479_8, 0.0_8]]
real(kind=8), parameter :: b(4)    =   [0.0_8,          1.208496649_8,  -0.644363171_8, 0.4358665215_8]
real(kind=8), parameter :: ai(4,4) = [ [0.0_8,         0.0_8,          0.0_8,          0.0_8],&
                                       [0.0_8,         0.4358665215_8, 0.0_8,          0.0_8],&
                                       [0.0_8,         0.2820667392_8, 0.4358665215_8, 0.0_8],&
                                       [0.0_8,         1.208496649_8,  -0.644363171_8, 0.4358665215_8]]

type, extends(timescheme_t), public :: ars343_t

    class(stvec_t),    allocatable :: y1,y2,y3,y4
    class(stvec_t),    allocatable :: q2,q3,q4
    class(stvec_t),    allocatable :: r,s

contains

    procedure, public :: step => step_ars343

end type ars343_t

contains

subroutine step_ars343(this, v0, operator, domain, dt)

    class(ars343_t),                    intent(inout) :: this
    class(stvec_t),         intent(inout) :: v0
    class(operator_t),      intent(inout) :: operator
    type(domain_t),         intent(in)    :: domain
    real(kind=8),           intent(in)    :: dt

    select type(operator)
    class is(imex_operator_t) !OK, let's work

    call operator%apply_explicit(this%y1,v0,domain)

    call this%r%assign(1.0_8,v0,ae(1,2)*dt,this%y1,domain)
    call operator%solve_implicit(this%s,this%r,ai(2,2)*dt,domain)
    call operator%apply_explicit(this%y2,this%s,domain)
    call this%q2%assign(1.0_8/(ai(2,2)*dt),this%s,-1.0_8/(ai(2,2)*dt),this%r,domain)

    call this%r%assign(     1.0_8,     v0,ae(1,3)*dt,this%y1,domain)
    call this%r%update(ae(2,3)*dt,this%y2,ai(2,3)*dt,this%q2,domain)
    call operator%solve_implicit(this%s,this%r,ai(3,3)*dt,domain)
    call operator%apply_explicit(this%y3,this%s,domain)
    call this%q3%assign(1.0_8/(ai(3,3)*dt),this%s,-1.0_8/(ai(3,3)*dt),this%r,domain)

    call this%r%assign(     1.0_8,     v0,ae(1,4)*dt,this%y1,domain)
    call this%r%update(ae(2,4)*dt,this%y2,ae(3,4)*dt,this%y3,domain)
    call this%r%update(ai(2,4)*dt,this%q2,ai(3,4)*dt,this%q3,domain)
    call operator%solve_implicit(this%s,this%r,ai(4,4)*dt,domain)
    call operator%apply_explicit(this%y4,this%s,domain)
    call this%q4%assign(1.0_8/(ai(4,4)*dt),this%s,-1.0_8/(ai(4,4)*dt),this%r,domain)

    call v0%update(b(2)*dt,this%y2,dt*b(2),this%q2,domain)
    call v0%update(b(3)*dt,this%y3,dt*b(3),this%q3,domain)
    call v0%update(b(4)*dt,this%y4,dt*b(4),this%q4,domain)

    class default !default operator
        call parcomm_global%abort("ARS343 can work only with imex operators")
    end select

end subroutine step_ars343

end module ars343_mod
