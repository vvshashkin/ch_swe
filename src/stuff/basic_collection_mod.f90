module basic_collection_mod

use linked_list_mod, only : linked_list_t
use string_mod,      only : string_t
use parcomm_mod,     only : parcomm_global

implicit none

type, extends(linked_list_t) :: basic_collection_t
    contains
        procedure, public :: get_string, get_r8
end type

contains

function get_string(this,name) result(str)
    class(basic_collection_t), intent(in) :: this
    character(len=*),          intent(in) :: name

    character(len=:), allocatable :: str

    class(*), pointer :: p

    call this%get_item(p,name)

    if(.not. associated(p)) call parcomm_global%abort(__FILE__//": get_str, collection item "//name//" not found")
    
    select type(p)
    type is (character(len=*))
        str = p
    type is (string_t)
        str = p%str
    class default
        call parcomm_global%abort(__FILE__//": get_str, linked list item "//name//" is not of character of string type")
    end select
end function

function get_r8(this,name) result(x)
    class(basic_collection_t), intent(in) :: this
    character(len=*),          intent(in) :: name

    real(kind=8) :: x

    class(*), pointer :: p

    call this%get_item(p,name)

    if(.not. associated(p)) call parcomm_global%abort(__FILE__//": get_str, collection item "//name//" not found")

    select type(p)
    type is (real(kind=8))
        x = p
    class default
        call parcomm_global%abort(__FILE__//": get_str, collection item "//name//" is not of real(kind=8) type")
    end select
end function

end module
