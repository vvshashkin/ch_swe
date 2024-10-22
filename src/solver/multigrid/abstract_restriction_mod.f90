module abstract_restriction_mod

use domain_mod, only : domain_t
use vector_mod, only : abstract_vector_t

implicit none

type, public, abstract :: restriction_operator_t
contains
    procedure(apply_restriction_i), public, deferred :: apply
end type restriction_operator_t

abstract interface
    subroutine apply_restriction_i(this, f_coarse, f_fine, domain_coarse, domain_fine)
        import restriction_operator_t, abstract_vector_t, domain_t
        class(restriction_operator_t), intent(inout) :: this
        class(abstract_vector_t),      intent(inout) :: f_coarse
        class(abstract_vector_t),      intent(inout) :: f_fine
        type(domain_t),                intent(in)    :: domain_coarse
        type(domain_t),                intent(in)    :: domain_fine
    end subroutine apply_restriction_i
end interface

contains

end module abstract_restriction_mod
