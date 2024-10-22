module grid_field_collection_mod

use linked_list_mod, only : linked_list_t
use grid_field_mod,  only : grid_field_t
use string_mod,      only : string_t

implicit none

type, public :: grid_field_collection_t
    private
    type(linked_list_t) :: grid_field_list
contains
    procedure, public :: add_grid_field
    procedure, public :: add_grid_field_reference
    procedure, public :: get_grid_field
    procedure, public :: get_grid_field_names
end type grid_field_collection_t

type, public :: grid_field_reference_t
    type(grid_field_t), pointer :: grid_field_pointer
contains
end type grid_field_reference_t

contains

subroutine add_grid_field(this, new_grid_field, grid_field_name)

    class(grid_field_collection_t), intent(inout) :: this
    type(grid_field_t),             intent(in)    :: new_grid_field
    character(len=*),               intent(in)    :: grid_field_name

    call this%grid_field_list%append(new_grid_field, grid_field_name)

end subroutine add_grid_field

subroutine add_grid_field_reference(this, grid_field, grid_field_name)

    class(grid_field_collection_t), intent(inout) :: this
    type(grid_field_t), pointer,    intent(in)    :: grid_field
    character(len=*),               intent(in)    :: grid_field_name

    type(grid_field_reference_t), pointer :: grid_field_reference

    allocate(grid_field_reference)
    grid_field_reference%grid_field_pointer => grid_field
    call this%grid_field_list%append(grid_field_reference, grid_field_name)

end subroutine add_grid_field_reference

subroutine get_grid_field(this, grid_field, grid_field_name)

    class(grid_field_collection_t), intent(in)    :: this
    type(grid_field_t), pointer,    intent(inout) :: grid_field
    character(len=*),               intent(in)    :: grid_field_name

    class(*), pointer :: data

    data => null()
    grid_field => null()

    call this%grid_field_list%get_item(data, grid_field_name)

    select type (data)
    type is (grid_field_t)
        grid_field => data
    type is (grid_field_reference_t)
        grid_field => data%grid_field_pointer
    end select

end subroutine get_grid_field

function get_grid_field_names(this) result(names)

    class(grid_field_collection_t), intent(in) :: this

    type(string_t), allocatable :: names(:)

    names = this%grid_field_list%get_item_names()

end function

end module grid_field_collection_mod
