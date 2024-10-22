module lsrk_mod

!source : J Niegemann 2012 Efficient low-storage Runge-Kutta schemes with optimized stability regions

use stvec_mod,      only: stvec_t
use timescheme_mod, only: timescheme_t
use operator_mod,   only: operator_t
use domain_mod,     only: domain_t

implicit none

!number of stages is 12
integer(kind=4), parameter :: method_size = 12

!arrays of constants for method
real(kind=8), parameter :: A(method_size) =                                    &
      [0.0000000000000000_8,  -0.0923311242368072_8,  -0.9441056581158819_8,   &
      -4.3271273247576394_8,  -2.1557771329026072_8,  -0.9770727190189062_8,   &
      -0.7581835342571139_8,  -1.7977525470825499_8,  -2.6915667972700770_8,   &
      -4.6466798960268143_8,  -0.1539613783825189_8,  -0.5943293901830616_8]
real(kind=8), parameter :: B(method_size) =                                    &
      [0.0650008435125904_8, 0.0161459902249842_8, 0.5758627178358159_8,       &
       0.1649758848361671_8, 0.3934619494248182_8, 0.0443509641602719_8,       &
       0.2074504268408778_8, 0.6914247433015102_8, 0.3766646883450449_8,       &
       0.0757190350155483_8, 0.2027862031054088_8, 0.2167029365631842_8]

private

type, public, extends(timescheme_t) :: lsrk_t
    class(stvec_t), allocatable :: k1, k2, temp
contains
    procedure, public :: step => step_lsrk
end type lsrk_t

contains

subroutine step_lsrk(this, v0, operator, domain, dt)

    class(lsrk_t),     intent(inout) :: this
    class(stvec_t),    intent(inout) :: v0
    class(operator_t), intent(inout) :: operator
    type(domain_t),    intent(in)    :: domain
    real(kind=8),      intent(in)    :: dt

    integer(kind=4)  :: i

    call this%k1%assign(1.0_8, v0, domain)
    call this%k2%assign(1.0_8, v0, domain)

    do i = 1, method_size
        call operator%apply(this%temp, this%k1, domain)
        call this%k2%assign(A(i), this%k2, dt, this%temp, domain)
        call this%k1%update(B(i), this%k2, domain)
    end do
    !answer=k1
    call v0%assign(1.0_8, this%k1, domain)

end subroutine step_lsrk

end module lsrk_mod

