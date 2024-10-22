module crank_nikolson_mod

use stvec_mod,      only : stvec_t
use timescheme_mod, only : timescheme_t
use operator_mod,   only : operator_t
use domain_mod,     only : domain_t

implicit none

private

type, public, extends(timescheme_t) :: crank_nikolson_t
    class(stvec_t), allocatable :: v_new
    contains
    procedure, public :: step => step_crank_nikolson
end type crank_nikolson_t

contains

subroutine step_crank_nikolson(this, v0, operator, domain, dt)

    class(crank_nikolson_t), intent(inout) :: this
    class(stvec_t),         intent(inout) :: v0
    class(operator_t),      intent(inout) :: operator
    type(domain_t),         intent(in)    :: domain
    real(kind=8),           intent(in)    :: dt

    call operator%apply(this%v_new, v0, domain)
    call this%v_new%assign(dt/2, this%v_new, 1.0_8, v0, domain)
    call operator%solve(v0, this%v_new, dt/2, domain)

end subroutine step_crank_nikolson

end module crank_nikolson_mod
