module operator_mod

use stvec_mod,   only : stvec_t
use domain_mod,  only : domain_t
use parcomm_mod, only : parcomm_global

implicit none

private

type, abstract, public :: operator_t
contains
    procedure(apply_i),  public, deferred :: apply !vout=A*vin
    procedure,           public           :: solve !vout=inverse(I-dt*A)*rhs
    procedure,           public           :: set_boundary_conditions
    procedure,           public           :: get_diagnostics
end type operator_t

type, abstract, extends(operator_t), public :: imex_operator_t
contains
    procedure(apply_imex),  public, deferred :: apply_explicit !vout=A*vin
    procedure(apply_imex),  public, deferred :: apply_implicit !vout=A*vin
    procedure,              public           :: apply_implicit_Jac !non-obligatory, but useful for testing
    procedure(solve_i),     public, deferred :: solve_implicit !vout=inverse(I-dt*A_impl)*rhs
    procedure(solve_i_Jac), public, deferred :: solve_implicit_Jac !vout=inverse(I-dt*J_impl)*rhs
end type imex_operator_t

type, abstract, public, extends(imex_operator_t) :: SISL_operator_t
contains
    procedure(make_SL_calculations), deferred :: make_SL_calculations
    procedure solve_implicit_Jac => solve_implicit_Jac_SISL    
end type

abstract interface
    subroutine apply_i(this, vout, vin, domain)
        import stvec_t, operator_t, domain_t
        class(operator_t), intent(inout) :: this
        class(stvec_t),    intent(inout) :: vout !inout to enable preallocated vectors
        class(stvec_t),    intent(inout) :: vin
        type(domain_t),    intent(in)    :: domain
    end subroutine apply_i
    subroutine apply_imex(this, vout, vin, domain)
        import stvec_t, imex_operator_t, domain_t
        class(imex_operator_t), intent(inout) :: this
        class(stvec_t),         intent(inout) :: vout !inout to enable preallocated vectors
        class(stvec_t),         intent(inout) :: vin
        type(domain_t),         intent(in)    :: domain
    end subroutine apply_imex
    subroutine solve_i(this, vout, rhs, dt, domain)
        import stvec_t, imex_operator_t, domain_t
        class(imex_operator_t), intent(inout) :: this
        class(stvec_t),         intent(inout) :: vout !inout to enable preallocated vectors
        class(stvec_t),         intent(inout) :: rhs
        real(kind=8),           intent(in)    :: dt
        type(domain_t),         intent(in)    :: domain
    end subroutine solve_i
    subroutine solve_i_Jac(this, vout, rhs, v0, dt, domain)
        import stvec_t, imex_operator_t, domain_t
        class(imex_operator_t), intent(inout) :: this
        class(stvec_t),         intent(inout) :: vout
        class(stvec_t),         intent(inout) :: rhs
        class(stvec_t),         intent(inout) :: v0 !Jac = dA/dv with v=v0
        real(kind=8),           intent(in)    :: dt
        type(domain_t),         intent(in)    :: domain
    end subroutine solve_i_Jac
    subroutine make_SL_calculations(this, vout, vin, wind_arr, wind_dp, dt, domain)
        import stvec_t, SISL_operator_t, domain_t
        class(SISL_operator_t), intent(inout) :: this
        class(stvec_t),         intent(inout) :: vout, vin, wind_arr, wind_dp
        real(kind=8),           intent(in)    :: dt
        type(domain_t),         intent(in)    :: domain
    end subroutine
end interface

contains

subroutine solve(this, vout, rhs, dt, domain)
    class(operator_t), intent(inout) :: this
    class(stvec_t),    intent(inout) :: vout !inout to enable preallocated vectors
    class(stvec_t),    intent(inout) :: rhs
    real(kind=8),      intent(in)    :: dt
    type(domain_t),    intent(in)    :: domain

    call parcomm_global%abort("Solve function not implemented for specific operator class")
end subroutine solve

subroutine set_boundary_conditions(this, stvec, domain)
    class(operator_t), intent(inout) :: this
    class(stvec_t),    intent(inout) :: stvec
    type(domain_t),    intent(in)    :: domain
    !by default, do nothing
end subroutine set_boundary_conditions

function get_diagnostics(this, v, domain) result(diagnostics)

    use key_value_mod, only : key_value_r8_t

    class(operator_t),     intent(inout) :: this
    class(stvec_t),        intent(inout) :: v
    type(domain_t),        intent(in)    :: domain

    type(key_value_r8_t)  :: diagnostics

    diagnostics = key_value_r8_t()

    call parcomm_global%abort("get_diagnostics function not implemented for specific operator class")
end function get_diagnostics

subroutine apply_implicit_Jac(this, vout, vin, v0, domain)
    class(imex_operator_t), intent(inout) :: this
    class(stvec_t),         intent(inout) :: vout
    class(stvec_t),         intent(inout) :: vin
    class(stvec_t),         intent(inout) :: v0 !Jac = dA/dv with v=v0
    type(domain_t),         intent(in)    :: domain

    call parcomm_global%abort("Solve function not implemented for specific operator class")
end subroutine apply_implicit_Jac

subroutine solve_implicit_Jac_SISL(this, vout, rhs, v0, dt, domain)

    class(SISL_operator_t), intent(inout) :: this
    class(stvec_t),         intent(inout) :: vout
    class(stvec_t),         intent(inout) :: rhs
    class(stvec_t),         intent(inout) :: v0 !Jac = dA/dv with v=v0
    real(kind=8),           intent(in)    :: dt
    type(domain_t),         intent(in)    :: domain

    call parcomm_global%abort("Solve implicit Jac function not implemented for specific operator class")

end subroutine solve_implicit_Jac_SISL

end module operator_mod
