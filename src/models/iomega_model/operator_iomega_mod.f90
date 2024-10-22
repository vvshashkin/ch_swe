module operator_iomega_mod

use operator_mod,     only: operator_t, imex_operator_t
use stvec_mod,        only: stvec_t
use stvec_iomega_mod, only: stvec_iomega_t
use domain_mod,       only: domain_t
use parcomm_mod,      only: parcomm_global


implicit none

private
public :: init_iomega_operator

type, public, extends(imex_operator_t) :: operator_iomega_t
    integer(kind=4) N
    complex(kind=8), allocatable :: omega(:) !eigen values, imag == oscillation frequency,
                                             !              real == amplification/decay
    complex(kind=8), allocatable :: omega_implicit(:)
    complex(kind=8), allocatable :: omega_explicit(:)
    contains
    procedure, public :: apply => apply_iomega
    procedure, public :: solve => solve_iomega

    procedure, public :: apply_explicit
    procedure, public :: apply_implicit
    procedure, public :: solve_implicit
    procedure, public :: solve_implicit_Jac
end type operator_iomega_t

contains

subroutine init_iomega_operator(operator, omega, omega_implicit)
    type(operator_iomega_t), allocatable :: operator
    complex(kind=8)                      :: omega(:)
    complex(kind=8), optional            :: omega_implicit(:)

    allocate(operator)
    operator%N = size(omega)
    operator%omega = omega

    if(present(omega_implicit)) then
        operator%omega_implicit = omega_implicit
        operator%omega_explicit = omega-omega_implicit
    else
        operator%omega_implicit = 0.0_8*omega(1:size(omega))
        operator%omega_explicit = omega
    end if
end

subroutine apply_iomega(this,vout,vin,domain)
    class(operator_iomega_t),  intent(inout) :: this
    class(stvec_t),            intent(inout) :: vout !inout to enable preallocated vectors
    class(stvec_t),            intent(inout) :: vin
    type(domain_t),            intent(in)    :: domain

    call apply_iomega_all(vout,vin,this%omega)

end subroutine apply_iomega

subroutine apply_explicit(this,vout,vin,domain)
    class(operator_iomega_t),  intent(inout) :: this
    class(stvec_t),            intent(inout) :: vout !inout to enable preallocated vectors
    class(stvec_t),            intent(inout) :: vin
    type(domain_t),            intent(in)    :: domain

    call apply_iomega_all(vout,vin,this%omega_explicit)

end subroutine apply_explicit

subroutine apply_implicit(this,vout,vin,domain)
    class(operator_iomega_t),  intent(inout) :: this
    class(stvec_t),            intent(inout) :: vout
    class(stvec_t),            intent(inout) :: vin
    type(domain_t),            intent(in)    :: domain

    call apply_iomega_all(vout,vin,this%omega_implicit)

end subroutine apply_implicit

subroutine apply_iomega_all(vout,vin,omega)
    class(stvec_t),   intent(inout) :: vout
    class(stvec_t),   intent(inout) :: vin
    complex(kind=8),  intent(in)    :: omega(:)

    integer(kind=4) :: i

    select type (vout)
    class is (stvec_iomega_t)
    select type (vin)
    class is (stvec_iomega_t)
        do i=1,size(omega)
            vout%f(i) = omega(i)*vin%f(i)
        end do
    class default
        call parcomm_global%abort("iomega operator failure: vin of wrong type")
    end select
    class default
        call parcomm_global%abort("iomega operator failure: vout of wrong type")
    end select

end subroutine apply_iomega_all

subroutine solve_iomega(this,vout,rhs,dt,domain)
    class(operator_iomega_t),  intent(inout) :: this
    class(stvec_t),            intent(inout) :: vout !inout to enable preallocated vectors
    class(stvec_t),            intent(inout) :: rhs
    real(kind=8),              intent(in)    :: dt
    type(domain_t),            intent(in)    :: domain

    call solve_iomega_all(vout,rhs,dt,this%omega)

end subroutine solve_iomega

subroutine solve_implicit(this,vout,rhs,dt,domain)
    class(operator_iomega_t),  intent(inout) :: this
    class(stvec_t),            intent(inout) :: vout !inout to enable preallocated vectors
    class(stvec_t),            intent(inout) :: rhs
    real(kind=8),              intent(in)    :: dt
    type(domain_t),            intent(in)    :: domain

    call solve_iomega_all(vout,rhs,dt,this%omega_implicit)

end subroutine solve_implicit

subroutine solve_implicit_Jac(this,vout,rhs,v0,dt,domain)
    class(operator_iomega_t),  intent(inout) :: this
    class(stvec_t),            intent(inout) :: vout !inout to enable preallocated vectors
    class(stvec_t),            intent(inout) :: rhs
    class(stvec_t),            intent(inout) :: v0
    real(kind=8),              intent(in)    :: dt
    type(domain_t),            intent(in)    :: domain

    call solve_iomega_all(vout,rhs,dt,this%omega_implicit)

end subroutine solve_implicit_jac

subroutine solve_iomega_all(vout, rhs, dt, omega)
    class(stvec_t),           intent(inout) :: vout !inout to enable preallocated vectors
    class(stvec_t),           intent(inout) :: rhs
    real(kind=8),             intent(in)    :: dt
    complex(kind=8),          intent(in)    :: omega(:)

    integer(kind=4) i

    select type (vout)
    class is (stvec_iomega_t)
    select type (rhs)
    class is (stvec_iomega_t)
        do i=1,size(omega)
            vout%f(i) = rhs%f(i) / (1.0_8 - dt*omega(i))
        end do
    class default
        call parcomm_global%abort("iomega solve-operator failure: rhs of wrong type")
    end select
    class default
        call parcomm_global%abort("iomega solve-operator failure: vout of wrong type")
    end select
end subroutine solve_iomega_all

end module operator_iomega_mod
