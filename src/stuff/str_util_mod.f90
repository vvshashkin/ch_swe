module str_util_mod

implicit none

character(len=*), parameter :: alpha    = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
character(len=*), parameter :: CAPS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
character(len=*), parameter :: lows = "abcdefghijklmnopqrstuvwxyz"
character(len=*), parameter :: digits   = "0123456789"


contains

function to_lower(str) result(out)
    character(len=*), intent(in)  :: str
    character(len=:), allocatable :: out
    
    integer(kind=4) :: l, i, k

    l = len(str)

    allocate(character(len=l) :: out)

    do i = 1, l
        k = index(CAPS,str(i:i))
        if(k == 0) then
            out(i:i) = str(i:i)
        else
            out(i:i) = lows(k:k)
        end if
    end do
    
end function 

function my_trim(str) result(str1) 
    character(len=*), intent(in)  :: str
    character(len=:), allocatable :: str1

    integer(kind=4) :: i

    do i = len(str),1,-1
        if(.not. is_blanc(str(i:i))) then
            str1 = str(1:i)
            return
        end if
    end do
    str1 = ""
end

function my_adjustl(str) result(str1)
    character(len=*), intent(in)  :: str
    character(len=:), allocatable :: str1

    integer(kind=4) :: i

    do i = 1, len(str)
        if(.not. is_blanc(str(i:i))) then
            str1 = str(i:)
            return
        end if
    end do
    str1 = ""
end

logical function is_blanc(c)
    character(len=1) :: c

    character(len=*), parameter :: blancs = " "//new_line('c')//char(9) !char(9) is \t

    is_blanc = (index(blancs,c) /= 0)
end

logical function all_blanc(str)
    character(len=*), intent(in) :: str

    integer(kind=4) :: i

    all_blanc = .false.
    do i=1, len(str)
        if(.not. is_blanc(str(i:i))) return
    end do
    all_blanc = .true.
end

function integer_to_str(i) result(str)
    integer(kind=4), intent(in) :: i
    character(len=:), allocatable :: str

    character(len=32) :: buff

    write(buff,'(I32)') i

    str = my_trim(my_adjustl(buff))
end function

end module