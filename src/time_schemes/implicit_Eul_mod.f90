module implicit_Eul_mod

use stvec_mod,      only : stvec_t
use timescheme_mod, only : timescheme_t
use operator_mod,   only : operator_t
use domain_mod,     only : domain_t

implicit none

private

type, public, extends(timescheme_t) :: implicit_Eul_t
    class(stvec_t), allocatable :: v_new
    contains
    procedure, public :: step => step_implicit_Eul
end type implicit_Eul_t

contains

subroutine step_implicit_Eul(this, v0, operator, domain, dt)

    class(implicit_Eul_t), intent(inout) :: this
    class(stvec_t),         intent(inout) :: v0
    class(operator_t),      intent(inout) :: operator
    type(domain_t),         intent(in)    :: domain
    real(kind=8),           intent(in)    :: dt

    call operator%solve(this%v_new, v0, dt, domain)
    call v0%assign(this%v_new,domain)
end subroutine step_implicit_Eul

end module implicit_Eul_mod
