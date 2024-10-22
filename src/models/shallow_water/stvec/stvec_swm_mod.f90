module stvec_swm_mod

use stvec_flexible_mod,  only : stvec_flexible_t, stvec_t
use grid_field_mod,      only : grid_field_t
use parcomm_mod,         only : parcomm_global
use domain_mod,          only : domain_t

implicit none

type, extends(stvec_flexible_t) :: stvec_swm_t
    type(grid_field_t), pointer :: h, u, v
contains
    procedure, public :: create_similar
    procedure, public :: associate_pointers
end type stvec_swm_t

contains

subroutine create_similar(this, destination)
    class(stvec_swm_t),          intent(in)    :: this
    class(stvec_t), allocatable, intent(inout) :: destination

    allocate(stvec_swm_t :: destination)
    select type (destination)
    class is (stvec_swm_t)

        call this%stvec_flexible_t%create_similar_allocated(destination)
        call destination%associate_pointers()

    class default
        call parcomm_global%abort("stvec_swm_t%create_similar: types error")
    end select
end subroutine create_similar

subroutine associate_pointers(this)
    class(stvec_swm_t), intent(inout) :: this

    call this%get_field(this%u,"u")
    call this%get_field(this%v,"v")
    call this%get_field(this%h,"h")

end subroutine

end module stvec_swm_mod
