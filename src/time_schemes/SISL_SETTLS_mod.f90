module SISL_SETTLS_mod

!semi implicit semi-Lagrangian time integration scheme as proposed by M.Hortal
!ECMWF technical memorandum 292, doi: 10.21957/ap8qpeghx
!For the state-vector v, operator A = L+N (linearized+non-linear)
!v^{n+1} - dt/2*(1+eps)*L*v^(n+1) =
!    = dt/2*N(v^n) + (v^n + dt/2*(1-eps)*L*v^n +dt/2*N(v)_extrapol)))_dep
!N(v)_extrapol = 2N(v^n)-N(v^(n-1))
!eps is decentering (epsilon in the type given below)

use stvec_mod,      only : stvec_t
use timescheme_mod, only : timescheme_t
use operator_mod,   only : operator_t, SISL_operator_t
use domain_mod,     only : domain_t
use parcomm_mod,    only : parcomm_global

implicit none

type, public, extends(timescheme_t) :: SISL_SETTLS_t
    class(stvec_t), allocatable :: v_dp, rhs
    class(stvec_t), allocatable :: explicit_tend, explicit_tend_old
    class(stvec_t), allocatable :: wind_arr, wind_dp
    logical, private :: first_step = .true.
    real(kind=8) :: epsilon = 0.0_8 ! time off-centering for implicit part, not used currently
    contains
    procedure, public :: step => step_SISL_SETTLS
end type SISL_SETTLS_t

contains

subroutine step_SISL_SETTLS(this, v0, operator, domain, dt)

    class(SISL_SETTLS_t),   intent(inout) :: this
    class(stvec_t),         intent(inout) :: v0
    class(operator_t),      intent(inout) :: operator
    type(domain_t),         intent(in)    :: domain
    real(kind=8),           intent(in)    :: dt

    select type(operator)
    class is (SISL_operator_t)

        call operator%apply_implicit(this%v_dp,v0,domain)
        call operator%apply_explicit(this%explicit_tend,v0,domain)

        if(this%first_step) then
            call this%v_dp%update(1.0_8,this%explicit_tend,domain)

            call this%wind_arr%assign(v0,domain)
            call this%wind_dp%assign(v0,domain)
            this%first_step = .false.
        else
            !use 2 -1 extrapolation
            call this%v_dp%update(2.0_8,this%explicit_tend,-1.0_8,this%explicit_tend_old,domain)
            !the same for wind
            call this%wind_dp%assign(-1.0_8,this%wind_arr,domain)
            call this%wind_arr%assign(v0,domain)
            call this%wind_dp%update(2.0_8,this%wind_arr,domain)
        end if

        call this%explicit_tend_old%assign(this%explicit_tend, domain)

        call this%v_dp%assign(1.0_8,v0,0.5_8*dt,this%v_dp,domain)

        call operator%make_SL_calculations(this%rhs,this%v_dp,this%wind_arr,this%wind_dp, dt, domain)

        call this%rhs%update(0.5_8*dt,this%explicit_tend,domain)
        call operator%solve_implicit(v0,this%rhs,0.5_8*dt,domain)

    class default
        call parcomm_global%abort("SISL_SETTLS scheme works only with SISL_operator_t")
    end select

end subroutine step_SISL_SETTLS

end module