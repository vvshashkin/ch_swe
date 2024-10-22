module smart_array_mod

use array_tools_mod, only : shrink_array, enlarge_array

implicit none

type, public :: smart_array_i4_t

    integer(kind=4) :: size_alloc, size_used
    integer(kind=4), allocatable :: a(:)

    contains
        procedure :: init => init_i4
        procedure :: append_element => append_i4_element
        procedure :: append_array => append_i4_array
        generic   :: append => append_element, append_array
        procedure :: shrink => shrink_i4
        procedure :: reset  => reset_i4
        procedure :: find_val => find_val_i4

end type

type, public :: smart_array_r8_t

    integer(kind=4) :: size_alloc, size_used
    real(kind=8), allocatable :: a(:)

    contains
        procedure :: init   => init_r8
        procedure :: append => append_r8
        procedure :: shrink => shrink_r8
        procedure :: reset  => reset_r8

end type

type, public :: vec_array_r8_t

    integer(kind=4) :: vec_size, size_alloc, size_used
    real(kind=8), allocatable :: a(:,:)

    contains
        procedure :: init   => init_r8_vec
        procedure :: append_element => append_r8_vec_element
        procedure :: append_array => append_r8_vec_array
        generic   :: append => append_array, append_element
        procedure :: shrink => shrink_r8_vec
        procedure :: reset  => reset_r8_vec

end type

!and not so smart container for 2d array
type array_container_i4_2d
    integer(kind=4), allocatable :: a(:,:)
end type

contains

subroutine init_i4(this, size_alloc)

    class(smart_array_i4_t), intent(out) :: this
    integer(kind=4),         intent(in)  :: size_alloc

    this%size_alloc = size_alloc
    this%size_used  = 0

    allocate(this%a(size_alloc))

end subroutine

subroutine append_i4_element(this, a, spare_size)

    class(smart_array_i4_t), intent(inout) :: this
    integer(kind=4),         intent(in)    :: a
    integer(kind=4),         intent(in), optional :: spare_size

    integer(kind=4) :: spare_size_loc

    if(this%size_used >= this%size_alloc) then

        spare_size_loc = 1
        if(present(spare_size)) spare_size_loc = spare_size

        call enlarge_array(this%a,this%size_alloc+spare_size_loc)
        this%size_alloc = this%size_alloc+spare_size_loc

    end if
    
    this%size_used = this%size_used+1
    this%a(this%size_used) = a

end subroutine

subroutine append_i4_array(this, a, spare_size)

    class(smart_array_i4_t), intent(inout) :: this
    integer(kind=4),         intent(in)    :: a(:)
    integer(kind=4),         intent(in), optional :: spare_size

    integer(kind=4) :: spare_size_loc, i

    if(this%size_used + size(a,1) >= this%size_alloc) then

        spare_size_loc = 1
        if(present(spare_size)) spare_size_loc = spare_size

        call enlarge_array(this%a,this%size_alloc+spare_size_loc+size(a,1)-1)
        this%size_alloc = this%size_alloc+spare_size_loc+size(a,1)-1

    end if
    
    do i = 1, size(a,1)
        this%a(this%size_used+i) = a(i)
    end do
    this%size_used = this%size_used+size(a,1)

end subroutine

subroutine shrink_i4(this)

    class(smart_array_i4_t), intent(inout) :: this

    call shrink_array(this%a, this%size_used)
    this%size_alloc = this%size_used

end subroutine

subroutine reset_i4(this)

    class(smart_array_i4_t), intent(inout) :: this

    this%size_used = 0

end subroutine

integer(kind=4) function find_val_i4(this, val) result(ind)

    class(smart_array_i4_t), intent(in) :: this
    integer(kind=4),         intent(in) :: val

    ind = findloc(this%a(1:this%size_used),val,dim=1)

end function

subroutine init_r8(this, size_alloc)

    class(smart_array_r8_t), intent(out) :: this
    integer(kind=4),         intent(in)  :: size_alloc

    this%size_alloc = size_alloc
    this%size_used  = 0

    allocate(this%a(size_alloc))

end subroutine

subroutine append_r8(this, a, spare_size)

    class(smart_array_r8_t), intent(inout) :: this
    real(kind=8),            intent(in)    :: a
    integer(kind=4),         intent(in), optional :: spare_size

    integer(kind=4) :: spare_size_loc

    spare_size_loc = 1
    if(present(spare_size)) spare_size_loc = spare_size

    if(this%size_used >= this%size_alloc) then
        call enlarge_array(this%a,this%size_alloc+spare_size_loc)
        this%size_alloc = this%size_alloc+spare_size_loc
    end if

    this%size_used = this%size_used+1
    this%a(this%size_used) = a

end subroutine

subroutine shrink_r8(this)

    class(smart_array_r8_t), intent(inout) :: this

    call shrink_array(this%a, this%size_used)
    this%size_alloc = this%size_used

end subroutine

subroutine reset_r8(this)

    class(smart_array_r8_t), intent(inout) :: this

    this%size_used = 0

end subroutine

subroutine init_r8_vec(this, vec_size, nvec_alloc)

    class(vec_array_r8_t), intent(out) :: this
    integer(kind=4),       intent(in)  :: vec_size, nvec_alloc

    this%vec_size = vec_size
    this%size_alloc = nvec_alloc
    this%size_used  = 0

    allocate(this%a(vec_size, nvec_alloc))

end subroutine

subroutine append_r8_vec_element(this, a, spare_size)

    class(vec_array_r8_t), intent(inout) :: this
    real(kind=8),          intent(in)    :: a(this%vec_size)
    integer(kind=4),       intent(in), optional :: spare_size

    integer(kind=4) :: spare_size_loc
    real(kind=8), allocatable :: tmp(:,:)

    if(this%size_used >= this%size_alloc) then

        spare_size_loc = 1
        if(present(spare_size)) spare_size_loc = spare_size

        allocate(tmp(this%vec_size,this%size_alloc+spare_size_loc))
        tmp(1:this%vec_size,1:this%size_used) = this%a(1:this%vec_size,1:this%size_used)

        deallocate(this%a)
        call move_alloc(tmp,this%a)

        this%size_alloc = this%size_alloc+spare_size_loc

    end if

    this%size_used = this%size_used+1
    this%a(1:this%vec_size,this%size_used) = a(1:this%vec_size)

end subroutine

subroutine append_r8_vec_array(this, a, spare_size)

    class(vec_array_r8_t), intent(inout) :: this
    real(kind=8),          intent(in)    :: a(:,:)
    integer(kind=4),       intent(in), optional :: spare_size

    integer(kind=4) :: spare_size_loc, n_elements, i
    real(kind=8), allocatable :: tmp(:,:)

    n_elements = size(a,2)

    if(this%size_used + n_elements >= this%size_alloc) then

        spare_size_loc = 1
        if(present(spare_size)) spare_size_loc = spare_size

        allocate(tmp(this%vec_size,this%size_alloc+spare_size_loc+n_elements-1))
        tmp(1:this%vec_size,1:this%size_used) = this%a(1:this%vec_size,1:this%size_used)

        deallocate(this%a)
        call move_alloc(tmp,this%a)

        this%size_alloc = this%size_alloc+spare_size_loc+n_elements-1

    end if

    do i = 1, n_elements
        this%a(1:this%vec_size,this%size_used+i) = a(1:this%vec_size,i)
    end do
    this%size_used = this%size_used+n_elements

end subroutine

subroutine shrink_r8_vec(this)

    class(vec_array_r8_t), intent(inout) :: this

    real(kind=8), allocatable :: tmp(:,:)

    allocate(tmp(this%vec_size,this%size_used))
    tmp(1:this%vec_size,1:this%size_used) = this%a(1:this%vec_size,1:this%size_used)

    deallocate(this%a)
    call move_alloc(tmp, this%a)

    this%size_alloc = this%size_used

end subroutine

subroutine reset_r8_vec(this)

    class(vec_array_r8_t), intent(inout) :: this

    this%size_used = 0

end subroutine

end module