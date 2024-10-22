module swm_forcing_mod

use domain_mod, only : domain_t
use stvec_mod,  only : stvec_t

implicit none

type, abstract :: swm_forcing_t

    contains
        procedure(get_tend_i), deferred :: get_tend
        procedure(apply_i),    deferred :: apply
        
end type

abstract interface

    subroutine get_tend_i(this, tend, state, domain)
        import swm_forcing_t, domain_t, stvec_t

        class(swm_forcing_t), intent(inout) :: this

        class(stvec_t),       intent(inout) :: tend

        class(stvec_t),       intent(inout) :: state
        type(domain_t),       intent(in)    :: domain
    end subroutine

    subroutine apply_i(this, state, dt, domain)
        import swm_forcing_t, domain_t, stvec_t

        class(swm_forcing_t), intent(inout) :: this

        class(stvec_t),       intent(inout) :: state

        real(kind=8),         intent(in)    :: dt
        type(domain_t),       intent(in)    :: domain
    end subroutine

end interface

end module