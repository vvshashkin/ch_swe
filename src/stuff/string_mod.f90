module string_mod

implicit none

type, public :: string_t
    character(:), allocatable :: str
contains
end type string_t

interface find_str_index
    module procedure :: find_str_index_str
    module procedure :: find_str_index_char
end interface

contains


subroutine print_strings(strings)

    type(string_t), intent(in) :: strings(:)

    integer(kind=4) :: i

    do i= 1, size(strings,1)
        print *, strings(i)%str
    end do

end subroutine

function strings(s01, s02, s03, s04, s05, s06, s07, s08, s09, s10, s11, &
                 s12, s13, s14, s15, s16, s17, s18, s19, s20, s21, s22, &
                 s23, s24, s25, s26, s27, s28, s29, s30, s31, s32, s33, &
                 s34, s35, s36, s37, s38, s39, s40, s41, s42, s43, s44, &
                 s45, s46, s47, s48, s49, s50, s51, s52, s53, s54, s55) &
                 result(string_array)
                 !fifty five will be enough? 

    character(len=*), intent(in) :: s01
    character(len=*), intent(in), optional :: &
                 s02, s03, s04, s05, s06, s07, s08, s09, s10, s11,      &
                 s12, s13, s14, s15, s16, s17, s18, s19, s20, s21, s22, &
                 s23, s24, s25, s26, s27, s28, s29, s30, s31, s32, s33, &
                 s34, s35, s36, s37, s38, s39, s40, s41, s42, s43, s44, &
                 s45, s46, s47, s48, s49, s50, s51, s52, s53, s54, s55

    type(string_t), allocatable :: string_array(:)
    type(string_t) :: str_buff(55)

    integer(kind=4) :: array_size

    array_size = 1
    str_buff(1)%str = s01

    if(present(s02)) then; array_size = 02; str_buff(02)%str = s02; end if
    if(present(s03)) then; array_size = 03; str_buff(03)%str = s03; end if
    if(present(s04)) then; array_size = 04; str_buff(04)%str = s04; end if
    if(present(s05)) then; array_size = 05; str_buff(05)%str = s05; end if
    if(present(s06)) then; array_size = 06; str_buff(06)%str = s06; end if
    if(present(s07)) then; array_size = 07; str_buff(07)%str = s07; end if
    if(present(s08)) then; array_size = 08; str_buff(08)%str = s08; end if
    if(present(s09)) then; array_size = 09; str_buff(09)%str = s09; end if
    if(present(s10)) then; array_size = 10; str_buff(10)%str = s10; end if
    if(present(s11)) then; array_size = 11; str_buff(11)%str = s11; end if
    if(present(s12)) then; array_size = 12; str_buff(12)%str = s12; end if
    if(present(s13)) then; array_size = 13; str_buff(13)%str = s13; end if
    if(present(s14)) then; array_size = 14; str_buff(14)%str = s14; end if
    if(present(s15)) then; array_size = 15; str_buff(15)%str = s15; end if
    if(present(s16)) then; array_size = 16; str_buff(16)%str = s16; end if
    if(present(s17)) then; array_size = 17; str_buff(17)%str = s17; end if
    if(present(s18)) then; array_size = 18; str_buff(18)%str = s18; end if
    if(present(s19)) then; array_size = 19; str_buff(19)%str = s19; end if
    if(present(s20)) then; array_size = 20; str_buff(20)%str = s20; end if
    if(present(s21)) then; array_size = 21; str_buff(21)%str = s21; end if
    if(present(s22)) then; array_size = 22; str_buff(22)%str = s22; end if
    if(present(s23)) then; array_size = 23; str_buff(23)%str = s23; end if
    if(present(s24)) then; array_size = 24; str_buff(24)%str = s24; end if
    if(present(s25)) then; array_size = 25; str_buff(25)%str = s25; end if
    if(present(s26)) then; array_size = 26; str_buff(26)%str = s26; end if
    if(present(s27)) then; array_size = 27; str_buff(27)%str = s27; end if
    if(present(s28)) then; array_size = 28; str_buff(28)%str = s28; end if
    if(present(s29)) then; array_size = 29; str_buff(29)%str = s29; end if
    if(present(s30)) then; array_size = 30; str_buff(30)%str = s30; end if
    if(present(s31)) then; array_size = 31; str_buff(31)%str = s31; end if
    if(present(s32)) then; array_size = 32; str_buff(32)%str = s32; end if
    if(present(s33)) then; array_size = 33; str_buff(33)%str = s33; end if
    if(present(s34)) then; array_size = 34; str_buff(34)%str = s34; end if
    if(present(s35)) then; array_size = 35; str_buff(35)%str = s35; end if
    if(present(s36)) then; array_size = 36; str_buff(36)%str = s36; end if
    if(present(s37)) then; array_size = 37; str_buff(37)%str = s37; end if
    if(present(s38)) then; array_size = 38; str_buff(38)%str = s38; end if
    if(present(s39)) then; array_size = 39; str_buff(39)%str = s39; end if
    if(present(s40)) then; array_size = 40; str_buff(40)%str = s40; end if
    if(present(s41)) then; array_size = 41; str_buff(41)%str = s41; end if
    if(present(s42)) then; array_size = 42; str_buff(42)%str = s42; end if
    if(present(s43)) then; array_size = 43; str_buff(43)%str = s43; end if
    if(present(s44)) then; array_size = 44; str_buff(44)%str = s44; end if
    if(present(s45)) then; array_size = 45; str_buff(45)%str = s45; end if
    if(present(s46)) then; array_size = 46; str_buff(46)%str = s46; end if
    if(present(s47)) then; array_size = 47; str_buff(47)%str = s47; end if
    if(present(s48)) then; array_size = 48; str_buff(48)%str = s48; end if
    if(present(s49)) then; array_size = 49; str_buff(49)%str = s49; end if
    if(present(s50)) then; array_size = 50; str_buff(50)%str = s50; end if
    if(present(s51)) then; array_size = 51; str_buff(51)%str = s51; end if
    if(present(s52)) then; array_size = 52; str_buff(52)%str = s52; end if
    if(present(s53)) then; array_size = 53; str_buff(53)%str = s53; end if
    if(present(s54)) then; array_size = 54; str_buff(54)%str = s54; end if
    if(present(s55)) then; array_size = 55; str_buff(55)%str = s55; end if

    string_array = str_buff(1:array_size)

end function

logical function is_in(string, strings)

    type(string_t), intent(in) :: string
    type(string_t), intent(in) :: strings(:)

    integer(kind=4) :: i

    is_in = .false.
    do i = 1, size(strings,1)
        if(string%str == strings(i)%str) then
            is_in = .true.
            return
        end if
    end do
end function

!find index of string in the array
integer(kind=4) function find_str_index_str(string, strings) result(ind)

    type(string_t), intent(in) :: string
    type(string_t), intent(in) :: strings(:)

    integer(kind=4) :: i

    ind = 0
    do i = 1, size(strings,1)
        if(string%str == strings(i)%str) then
            ind = i
            return
        end if
    end do
end function

integer(kind=4) function find_str_index_char(string, strings) result(ind)

    character(len=*), intent(in) :: string
    type(string_t),   intent(in) :: strings(:)

    integer(kind=4) :: i

    ind = 0
    do i = 1, size(strings,1)
        if(string == strings(i)%str) then
            ind = i
            return
        end if
    end do
end function

end module
