module time_util_mod

implicit none

integer(kind=4), parameter :: Day_sec = 86400, Hour_sec = 3600, Min_sec = 60

contains

function get_day_hour_min_sec_str(seconds) result(str)
    real(kind=8), intent(in) :: seconds
    character(len=:), allocatable :: str

    character(len=128) :: buff1
    character(len=256) :: buff2
    integer(kind=4) :: days, hours, minutes, seconds1

    buff2 = ""
    seconds1 = int(seconds)
    days = int(seconds1 / Day_sec)
    seconds1 = mod(seconds1, Day_sec)
    if(days > 0) then
        write(buff1,"(I8,A)") days, "d"
        buff2 = trim(buff2)//trim(buff1)
    end if

    hours = int(seconds1 / Hour_sec)
    seconds1 = mod(seconds1, Hour_sec)
    if(hours > 0) then
        write(buff1,"(I3,A)") hours,"hr"
        buff2 = trim(buff2)//trim(buff1)
    end if

    minutes = int(seconds1 / Min_sec)
    seconds1 = mod(seconds1, Min_sec)
    if(minutes > 0) then
        write(buff1,"(I3,A)") minutes,"min"
        buff2 = trim(buff2)//trim(buff1)
    end if

    write(buff1,"(F6.2,A)") seconds1+(seconds-int(seconds)),"sec"
    str = trim(buff2)//trim(buff1)
end function

end module
