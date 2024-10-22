module error_mod

type, public :: error_t
    logical :: raised = .false.
    character(len=:), allocatable :: message, location
    contains
       procedure, public :: raise
end type

contains

subroutine raise(this,message,location)
    class(error_t), intent(out) :: this
    character(len=*), optional, intent(in) :: message, location

    this%raised = .true.

    this%location = ""
    if(present(location)) this%location = location

    this%message = ""
    if(present(message)) this%message = message

end subroutine

end module