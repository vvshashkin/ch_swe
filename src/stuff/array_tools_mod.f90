module array_tools_mod

implicit none

type array_container_1d_r8_t
    real(kind=8), allocatable :: a(:)
end type

interface shrink_array
    module procedure :: shrink_1d_int4_array
    module procedure :: shrink_1d_real8_array
    !add another options if you need it
end interface

interface enlarge_array
    module procedure :: enlarge_1d_real8_array
    module procedure :: enlarge_1d_int4_array
end interface

interface derealloc
    module procedure :: derealloc_i4
    module procedure :: derealloc_r8
end interface

contains

subroutine shrink_1d_int4_array(a,to_length)
        
    integer(kind=4), allocatable, intent(inout) :: a(:)
    integer(kind=4),              intent(in)    :: to_length

    integer(kind=4), allocatable :: tmp(:)

    tmp = a(1:to_length)
    deallocate(a)
    call move_alloc(tmp,a)

end subroutine

subroutine shrink_1d_real8_array(a,to_length)
        
    real(kind=8), allocatable, intent(inout) :: a(:)
    integer(kind=4),           intent(in)    :: to_length

    real(kind=8), allocatable :: tmp(:)

    tmp = a(1:to_length)
    deallocate(a)
    call move_alloc(tmp,a)

end subroutine

subroutine enlarge_1d_real8_array(a, to_length)

    real(kind=8), allocatable, intent(inout) :: a(:)
    integer(kind=4),           intent(in)    :: to_length

    real(kind=8)    :: tmp(size(a,1))
    integer(kind=4) :: current_size

    current_size = size(a,1)

    tmp(1:current_size) = a(1:current_size)
    deallocate(a)
    allocate(a(to_length))
    a(1:current_size) = tmp(1:current_size)

end subroutine

subroutine enlarge_1d_int4_array(a, to_length)

    integer(kind=4), allocatable, intent(inout) :: a(:)
    integer(kind=4),              intent(in)    :: to_length

    integer(kind=4) :: tmp(size(a,1))
    integer(kind=4) :: current_size

    current_size = size(a,1)

    tmp(1:current_size) = a(1:current_size)
    deallocate(a)
    allocate(a(to_length))
    a(1:current_size) = tmp(1:current_size)

end subroutine

subroutine derealloc_r8(a, new_length)

    real(kind=8), allocatable, intent(inout) :: a(:)
    integer(kind=4),           intent(in)    :: new_length

    deallocate(a)
    allocate(a(new_length))

end subroutine

subroutine derealloc_i4(a, new_length)

    integer(kind=4), allocatable, intent(inout) :: a(:)
    integer(kind=4),              intent(in)    :: new_length

    deallocate(a)
    allocate(a(new_length))

end subroutine

end module