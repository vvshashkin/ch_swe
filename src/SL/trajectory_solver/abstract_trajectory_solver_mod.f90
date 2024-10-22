module abstract_trajectory_solver_mod

use domain_mod,          only : domain_t
use stvec_flexible_mod,  only : stvec_flexible_t

implicit none

type, abstract :: trajectory_solver_t
    contains
        procedure(find_departure_points), deferred :: find_departure_points
end type

abstract interface

    subroutine find_departure_points(this, dep_points_coords, dt, wind_arr, wind_dep, domain)
        import trajectory_solver_t, stvec_flexible_t, domain_t
        class(trajectory_solver_t), intent(inout) :: this
        class(stvec_flexible_t),    intent(inout) :: dep_points_coords
        real(kind=8),               intent(in)    :: dt
        class(stvec_flexible_t),    intent(inout) :: wind_arr, wind_dep
        type(domain_t),             intent(in)    :: domain
    end subroutine

end interface

contains

end module