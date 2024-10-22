module abstract_scorer_mod

use stvec_mod,  only : stvec_t
use domain_mod, only : domain_t

implicit none

type, abstract :: scorer_t
    contains
        procedure(get_scores_i), deferred :: print_scores
end type

abstract interface

    subroutine get_scores_i(this, state, domain, time)

        import scorer_t, stvec_t, domain_t

        class(scorer_t), intent(inout) :: this
        class(stvec_t),  intent(in)    :: state
        type(domain_t),  intent(in)    :: domain
        real(kind=8),    intent(in)    :: time

    end subroutine

end interface

end module