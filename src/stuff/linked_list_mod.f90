module linked_list_mod

use string_mod,   only : string_t
use parcomm_mod,  only : parcomm_global

implicit none

type, public :: linked_list_node_t
    type(linked_list_node_t),  pointer     :: next  => null() ! pointer to next list node
    class(*),                  pointer     :: value => null()
    character(len=:),          allocatable :: name            ! names to reference list nodes
contains
end type linked_list_node_t

type, public :: linked_list_t

    integer(kind=4) :: length = 0
    !WORKAROUND #safe_final
    logical         :: safe_final = .true.

    type(linked_list_node_t), pointer :: head => null()
    type(linked_list_node_t), pointer :: tail => null()
contains
    procedure, public  :: append
    procedure, public  :: get_item
    procedure, public  :: get_item_names
    procedure, public  :: change_item
    procedure, public  :: copy_to
    final              :: finalize_linked_list

end type linked_list_t

interface
    subroutine node_func(node)
        import linked_list_node_t
        type(linked_list_node_t), pointer, intent(inout) :: node
    end subroutine node_func
end interface

contains

subroutine append(this, data, data_name)

    class(linked_list_t), intent(inout) :: this
    class(*),             intent(in)    :: data
    character(len=*),     intent(in)    :: data_name

    type(linked_list_node_t), pointer :: new_item

    allocate(new_item)
    allocate(new_item%value, source = data)
    allocate(new_item%name,  source = data_name)
    new_item%next => null()

    this%length = this%length + 1

    if (.not. associated(this%head)) then
        this%head => new_item
        this%tail => new_item
    else
        this%tail%next => new_item
        this%tail      => new_item
    end if

end subroutine append

subroutine get_item(this, out, name)
    class(linked_list_t), intent(in)           :: this
    class(*),             intent(out), pointer :: out
    character(len=*),     intent(in)           :: name

    type(linked_list_node_t), pointer :: current_node

    current_node => this%head
    out => null()

    do
        if (.not. associated(current_node)) exit
        if (current_node%name == name) then
            out => current_node%value
            exit
        end if
        current_node => current_node%next
    end do

end subroutine get_item

function get_item_names(this) result(names)

    class(linked_list_t), intent(in) :: this

    type(string_t), allocatable :: names(:)

    type(linked_list_node_t), pointer :: current
    integer(kind=4)                   :: i

    allocate(names(this%length))

    current => this%head
    do i=1, this%length
        names(i)%str = current%name
        current => current%next
    end do

end function get_item_names

subroutine change_item(this, data, data_name)

    class(linked_list_t), intent(inout) :: this
    class(*),             intent(in)    :: data
    character(len=*),     intent(in)    :: data_name

    class(*), pointer :: data_old => null()

    call this%get_item(data_old,data_name)

    if(.not. associated(data_old)) return

    deallocate(data_old)
    allocate(data_old,source=data)

end subroutine change_item

subroutine copy_to(this, destination)
    class(linked_list_t), intent(in)  :: this
    class(linked_list_t), intent(out) :: destination

    class(*), pointer :: data
    type(linked_list_node_t), pointer :: current_node

    current_node => this%head
    do
        if (.not. associated(current_node)) exit

        call destination%append(current_node%value, current_node%name)

        current_node => current_node%next
    end do

end subroutine copy_to

impure elemental subroutine finalize_linked_list(this)

    type(linked_list_t), intent(inout) :: this

    type(linked_list_node_t), pointer :: tmp

    !WORKAROUND #safe_final
    if(.not. this%safe_final) return

    do
        if (.not. associated(this%head)) exit

        tmp => this%head
        this%head => this%head%next

        deallocate(tmp%value)
        deallocate(tmp%name)
        deallocate(tmp)
    end do

    this%length = 0

    nullify(this%head)
    nullify(this%tail)

end subroutine finalize_linked_list

end module linked_list_mod
