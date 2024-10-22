module abstract_prolongation_mod

use domain_mod, only : domain_t
use vector_mod, only : abstract_vector_t

implicit none

type, public, abstract :: prolongation_operator_t
contains
    procedure(apply_prolongation_i), public, deferred :: apply
end type prolongation_operator_t

abstract interface
    subroutine apply_prolongation_i(this, f_fine, f_coarse, domain_fine, domain_coarse)
        import prolongation_operator_t, abstract_vector_t, domain_t
        class(prolongation_operator_t), intent(inout) :: this
        class(abstract_vector_t),       intent(inout) :: f_coarse
        class(abstract_vector_t),       intent(inout) :: f_fine
        type(domain_t),                 intent(in)    :: domain_coarse
        type(domain_t),                 intent(in)    :: domain_fine
    end subroutine apply_prolongation_i
end interface

contains

end module abstract_prolongation_mod
