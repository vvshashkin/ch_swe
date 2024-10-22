module orography_mod

use grid_field_mod,  only : grid_field_t
use linked_list_mod, only : linked_list_t

implicit none

type orography_t
    type(grid_field_t) :: h, dh_alpha, dh_beta
end type orography_t

type orography_collection_t
    ! type(orography_1mesh_t) :: o, x, y, xy
    private
    type(linked_list_t) :: orography_list
contains
    procedure, public :: add_orography
    procedure, public :: get_orography
end type orography_collection_t

contains

subroutine add_orography(this, new_orography, orography_name)

    class(orography_collection_t), intent(inout) :: this
    class(orography_t),            intent(in)    :: new_orography
    character(len=*),              intent(in)    :: orography_name

    call this%orography_list%append(new_orography, orography_name)

end subroutine add_orography

subroutine get_orography(this, orography, orography_name)

    class(orography_collection_t),     intent(in)    :: this
    type(orography_t),        pointer, intent(inout) :: orography
    character(len=*),                  intent(in)    :: orography_name

    class(*), pointer :: data

    data => null()
    orography => null()

    call this%orography_list%get_item(data, orography_name)

    if(.not. associated(data)) return

    select type (data)
    type is (orography_t)
        orography => data
    end select

end subroutine get_orography

end module orography_mod
