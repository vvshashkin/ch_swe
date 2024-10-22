module SL_UKMO_mod

!UKMO-END_GAME / Canadian GEM semi-Lagrangian  uterative implicit
!time-integration scheme
!10.1515/caim-2016-0020 - ENDGame
!https://doi.org/10.1175/MWR-D-18-0438.1 - GEM
!summarized in Eqs. (1.5-1.6) of https://doi.org/10.1515/rnam-2021-0020

use stvec_mod,      only : stvec_t
use timescheme_mod, only : timescheme_t
use operator_mod,   only : operator_t, SISL_operator_t
use domain_mod,     only : domain_t
use parcomm_mod,    only : parcomm_global

use timer_mod,      only : start_timer, stop_timer

implicit none

type, public, extends(timescheme_t) :: SL_UKMO_t
    class(stvec_t), allocatable :: v_dp, rhs
    class(stvec_t), allocatable :: wind_arr, wind_dp
    integer(kind=4) :: num_iter = 3
    real(kind=8)    :: epsilon = 0.01_8 !time off-centering
    contains
    procedure, public :: step => step_SL_UKMO
end type SL_UKMO_t

contains

subroutine step_SL_UKMO(this, v0, operator, domain, dt)

    class(SL_UKMO_t),       intent(inout) :: this
    class(stvec_t),         intent(inout) :: v0
    class(operator_t),      intent(inout) :: operator
    type(domain_t),         intent(in)    :: domain
    real(kind=8),           intent(in)    :: dt

    integer(kind=4) :: it
    real(kind=8) :: alpha, beta

    alpha = (1+this%epsilon)
    beta  = (1-this%epsilon)

    select type(operator)
    class is (SISL_operator_t)

        call start_timer("explicit_tend")
        call operator%apply(this%v_dp,v0,domain)
        call this%v_dp%assign(0.5_8*beta*dt,this%v_dp,1.0_8,v0,domain)

        call this%wind_dp%assign(v0,domain)
        call stop_timer("explicit_tend")

        do it = 1, this%num_iter
            call start_timer("sl_calc")
            call this%wind_arr%assign(v0,domain)
            call operator%make_SL_calculations(this%rhs,this%v_dp,this%wind_arr,this%wind_dp, dt, domain)
            call stop_timer("sl_calc")

            call start_timer("solver")
            call operator%solve(v0,this%rhs,0.5_8*alpha*dt,domain)
            call stop_timer("solver")
        end do

    class default
        call parcomm_global%abort("SISL_SETTLS scheme works only with SISL_operator_t")
    end select

end subroutine step_SL_UKMO

end module