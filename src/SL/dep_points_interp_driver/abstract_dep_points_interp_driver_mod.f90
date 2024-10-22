module abstract_dep_points_interp_driver_mod

use stvec_flexible_mod, only : stvec_flexible_t
use domain_mod,         only : domain_t

implicit none

type, abstract, public :: dep_points_interp_driver_t
    contains
        procedure(do_i), deferred :: do
end type

abstract interface
    subroutine do_i(this, q_dep, q, dep_points_coords, domain)
        import dep_points_interp_driver_t, stvec_flexible_t, domain_t
        class(dep_points_interp_driver_t), intent(inout) :: this
        class(stvec_flexible_t),           intent(inout) :: q_dep, q, dep_points_coords
        type(domain_t),                    intent(in)    :: domain
    end subroutine
end interface

end module