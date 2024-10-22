module geosci_config_parser_mod

use generic_config_mod, only : generic_config_t
use string_mod,         only : string_t
use str_util_mod,       only : my_adjustl, my_trim, is_blanc, all_blanc, to_lower, &
                               integer_to_str
use const_mod,          only : pi, Earth_grav, Earth_omega, Earth_radii, &
                               Day24h_sec

implicit none

character(len=*), parameter :: digits = "0123456789"
character(len=*), parameter :: exp_char = "eEdD"
character(len=*), parameter :: plus_minus = "+-"
character(len=*), parameter :: comment_char = "#"
character(len=*), parameter :: valid_token_end = ","//new_line("c")//comment_char

class(generic_config_t), allocatable :: constants

contains

subroutine parse_geosci_config(config,err_msg,str,path,root_config)
    class(generic_config_t), intent(inout), target :: config
    character(len=:),        intent(out), allocatable :: err_msg
    character(len=*),        intent(in)    :: str
    character(len=*),        intent(in),    optional :: path
    class(generic_config_t), intent(inout), optional, target :: root_config

    character(len=:), allocatable :: clean_str
    integer(kind=4) :: pos, pos1
    class(generic_config_t), pointer :: root_config_loc

    if( .not. allocated(constants)) then
        constants = config%get_empty_subconfig()
        call constants%set("pi",pi)
        call constants%set("Earth_grav",Earth_grav)
        call constants%set("Earth_omega",Earth_omega)
        call constants%set("Earth_radii",Earth_radii)
        call constants%set("day",Day24h_sec)
        call constants%set("min",60)
        call constants%set("sec",1)
        call constants%set("hour",3600)
    end if

    root_config_loc => config
    if(present(root_config)) root_config_loc => root_config

    clean_str = my_trim(my_adjustl(str))

    pos = 0
    do while (pos < len(clean_str))
        call parse_var_value(config,pos1,err_msg,clean_str(pos+1:),root_config_loc,path)
        if(allocated(err_msg)) return
        pos = pos+pos1
        if(pos <= len(clean_str)) then
            if(clean_str(pos:pos) == comment_char) pos = max(pos-1,1)
        end if
    end do
end subroutine

subroutine parse_var_value(config, end_pos, err_msg, str, root_config, path)
    character(len=*),        intent(in)               :: str
    class(generic_config_t), intent(inout)            :: root_config
    character(len=*),        intent(in),     optional :: path
    class(generic_config_t), intent(inout)            :: config
    integer(kind=4),         intent(out)              :: end_pos
    character(len=:),        intent(out), allocatable :: err_msg

    character(len=:), allocatable :: var_name, test_comment
    integer(kind=4) :: ind_eq
    logical :: success

    ind_eq    = index(str,"=")

    test_comment = my_adjustl(str)
    if(test_comment(1:1) == comment_char) then
        end_pos = index(str,new_line('c'))
        if(end_pos == 0) end_pos = len(str)+1
        return
    end if

    var_name = my_trim(my_adjustl(str(1:ind_eq-1)))
    if(var_name == "") then
        err_msg = "empty variable name, '"//str//"'"
        return
    end if
    if(index(var_name,new_line("c")) /= 0 .or. &
       index(var_name," ") /= 0) then
        err_msg = "incorrect variable name, '"//str//"'"
        return
    end if

    if(present(path)) var_name = path//"%"//var_name

    call parse_value(config,end_pos, err_msg, var_name, str(ind_eq+1:), root_config, path)
    end_pos = ind_eq+end_pos

end subroutine

subroutine parse_value(config, end_pos, err_msg, var_name, str, root_config, path)
    character(len=*),        intent(in)               :: str
    class(generic_config_t), intent(inout)            :: root_config
    character(len=*),        intent(in)               :: var_name
    character(len=*),        intent(in),     optional :: path
    class(generic_config_t), intent(inout)            :: config
    integer(kind=4),         intent(out)              :: end_pos
    character(len=:),        intent(out), allocatable :: err_msg

    logical :: success
  
    call parse_int(config, end_pos, success, str, var_name)
    if(success) return

    call parse_real(config, end_pos, success, str, var_name)
    if(success) return

    call parse_logical(config, end_pos, success, str, var_name)
    if(success) return

    call parse_str(config, end_pos, success, err_msg, str, var_name)
    if(allocated(err_msg)) return
    if(success) return

    call parse_subconfig(config, end_pos, success, err_msg, str, var_name,root_config)
    if(allocated(err_msg)) return
    if(success) return

    call parse_array(config, end_pos, success, err_msg, str, var_name, root_config)
    if(allocated(err_msg) .or. success) return

    call parse_arithmetic_expression(config, end_pos, success, err_msg, str, var_name, root_config)
    if(allocated(err_msg)) return
    if(success) return

    err_msg = "wrong config entry: '"//str//"'"

end subroutine

subroutine parse_int(config, end_pos, success, str, var_name)
    class(generic_config_t), intent(inout) :: config
    integer(kind=4),         intent(out)   :: end_pos
    logical,                 intent(out)   :: success
    character(len=*),        intent(in)    :: str, var_name

    integer(kind=4) :: i
    character(len=:), allocatable :: str1

    end_pos = find_token_end(str)
    str1 = my_trim(my_adjustl(str(1:end_pos-1)))

    success = .false.

    do i = 1, len(str1)
        if(index(digits,str1(i:i))==0) return
    end do

    success = .true.
    read(str1,*) i
    call config%set(var_name,i)

end subroutine

subroutine parse_real(config, end_pos, success, str, var_name)
    class(generic_config_t), intent(inout) :: config
    integer(kind=4),         intent(out)   :: end_pos
    logical,                 intent(out)   :: success
    character(len=*),        intent(in)    :: str, var_name

    integer(kind=4) :: i, n_dots, n_exp, exp_pos
    character(len=:), allocatable :: str1
    real(kind=8) :: p

    end_pos = find_token_end(str)
    str1 = my_trim(my_adjustl(str(1:end_pos-1)))

    success = .false.
    n_exp = 0
    n_dots = 0
    exp_pos = -1000

    do i = 1, len(str1)
        if(index(digits,str1(i:i)) /=0) then
            !accept and move to next
        else if(index(exp_char,str1(i:i)) /= 0 .and. n_exp == 0) then
            n_exp = 1
            exp_pos = i
        else if(i /= 1 .and. index(plus_minus,str1(i:i)) /= 0 .and. &
                exp_pos == i-1) then
                !continue
        else if(str1(i:i)=="." .and. n_exp == 0 .and. n_dots == 0) then
            n_dots = 1
        else
            return
        end if
    end do

    success = .true.
    read(str1,*) p
    call config%set(var_name,p)

end subroutine

subroutine parse_logical(config, end_pos, success, str, var_name)
    class(generic_config_t), intent(inout) :: config
    integer(kind=4),         intent(out)   :: end_pos
    logical,                 intent(out)   :: success
    character(len=*),        intent(in)    :: str, var_name

    integer(kind=4) :: i
    logical         :: logical_val
    character(len=:), allocatable :: str1

    end_pos = find_token_end(str)
    str1 = to_lower(my_trim(my_adjustl(str(1:end_pos-1))))

    success = .false.

    if(str1 == ".true.") then
        success = .true.
        call config%set(var_name,.true.)
    else if(str1 == ".false.") then
        success = .true.
        call config%set(var_name,.false.)
    end if

end subroutine

subroutine parse_str(config, end_pos, success, err_msg, str, var_name)
    class(generic_config_t), intent(inout) :: config
    integer(kind=4),         intent(out)   :: end_pos
    logical,                 intent(out)   :: success
    character(len=:),        intent(out), allocatable :: err_msg
    character(len=*),        intent(in)    :: str, var_name

    integer(kind=4) :: i, j, k
    logical         :: logical_val
    character(len=len(str)) :: str1
    character(len=1) :: openning

    success = .false.

    do i = 1, len(str)
        if(str(i:i) == "'" .or. str(i:i) == '"') exit
        if(.not. is_blanc(str(i:i))) return
    end do

    openning = str(i:i)

    j = 0
    do i = i+1, len(str)
        if(str(i:i) == openning) then
            success = .true.
            exit
        end if
        j = j+1
        str1(j:j) = str(i:i)
    end do

    if(.not. success) return

    do i = i+1, len(str)
        if(index(valid_token_end,str(i:i)) /= 0) exit

        if(.not. is_blanc(str(i:i))) then
            success = .false.
            err_msg = "trailing characters after "//openning// " : "//str
            return
        end if
    end do

    end_pos = i
    call config%set(var_name, string_t(str1(1:j)))

end subroutine

subroutine parse_subconfig(config, end_pos, success, err_msg, str, var_name, root_config)
    class(generic_config_t), intent(inout) :: config
    integer(kind=4),         intent(out)   :: end_pos
    logical,                 intent(out)   :: success
    character(len=:),        intent(out), allocatable :: err_msg
    character(len=*),        intent(in)    :: str, var_name
    class(generic_config_t), intent(inout) :: root_config

    integer(kind=4) :: i, n_brackets, first_bracket, last_bracket
    class(generic_config_t), allocatable :: subconfig
    character(len=:), allocatable :: str1

    success = .false.

    do i = 1, len(str)
        if(str(i:i) == "{") exit
        if(.not. is_blanc(str(i:i))) return
    end do

    first_bracket = i
    n_brackets = 1
    do i = i+1,len(str)
        if(str(i:i) == "{") n_brackets = n_brackets+1
        if(str(i:i) == "}") n_brackets = n_brackets-1

        if(n_brackets == 0) exit
    end do

    if(n_brackets /= 0) then
        err_msg = "unbalanced {}: '"//str//"'"
        return
    end if

    last_bracket = i

    do i = i+1, len(str)
        if(index(valid_token_end,str(i:i)) /= 0) exit

        if(.not. is_blanc(str(i:i))) then
            err_msg = "trailing characters after '}':  '"//str//"'"
            return
        end if
    end do

    success = .true.
    end_pos = i

    call parse_geosci_config(config,err_msg,str(first_bracket+1:last_bracket-1),path = var_name,root_config = root_config)

end subroutine

subroutine parse_array(config, end_pos, success, err_msg, str, var_name, root_config)
    class(generic_config_t), intent(inout) :: config
    integer(kind=4),         intent(out)   :: end_pos
    logical,                 intent(out)   :: success
    character(len=:),        intent(out), allocatable :: err_msg
    character(len=*),        intent(in)    :: str, var_name
    class(generic_config_t), intent(inout) :: root_config

    integer(kind=4) :: i, first_bracket, last_bracket, pos, pos1
    character(len=:), allocatable :: str1

    success = .false.

    do i = 1, len(str)
        if(str(i:i) == "[") exit
        if(.not. is_blanc(str(i:i))) return
    end do

    first_bracket = i
    do i = i+1,len(str)
        if(str(i:i) == "[") then
            err_msg = "parse geosci config error, multi-level lists not allowed "// str
            return
        end if
        if(str(i:i) == "]") exit
    end do

    if(str(i:i) /= "]") then
        err_msg = "unbalanced []: '"//str//"'"
        return
    end if

    last_bracket = i

    do i = i+1, len(str)
        if(index(valid_token_end,str(i:i)) /= 0) exit

        if(.not. is_blanc(str(i:i))) then
            err_msg = "trailing characters after ']':  '"//str//"'"
            return
        end if
    end do

    end_pos = i

    str1 = my_trim(my_adjustl(str(first_bracket+1:last_bracket-1)))
    i  = 1
    pos = 0
    do while(pos < len(str1))
        call parse_value(config,pos1, err_msg, &
                         var_name//"%"//integer_to_str(i), &
                         str1(pos+1:), root_config)
        i = i+1
        pos = pos+pos1
        if(allocated(err_msg)) return
    end do

    call config%set(var_name//"%length",i-1)

    success = .true.

end subroutine

subroutine parse_arithmetic_expression(config, end_pos, success, err_msg, str, var_name, root_config)
 
    use arithmetic_parser_mod,     only : evaluate_expression
    use arithmetic_parser_val_mod, only : value_t, int_value_t, real_value_t
    use error_mod,                 only : error_t
    
    class(generic_config_t), intent(inout) :: config
    integer(kind=4),         intent(out)   :: end_pos
    logical,                 intent(out)   :: success
    character(len=:),        intent(out), allocatable :: err_msg
    character(len=*),        intent(in)    :: str, var_name
    class(generic_config_t), intent(inout) :: root_config

    character(len=:), allocatable :: str1
    class(value_t),   allocatable :: val
    type(error_t) :: err

    end_pos = find_token_end(str)
    str1 = my_trim(my_adjustl(str(1:end_pos-1)))

    call evaluate_expression(val,err,str1,tmp_vars= config, vars = root_config, global_vars = constants)
    if(err%raised) then
        err_msg = "Arithmetic parser error: " // err%message
        success = .false.
        return
    end if

    success = .true.

    select type (val)
    type is (int_value_t)
        call config%set(var_name,val%val)
    type is (real_value_t)
        call config%set(var_name,val%val)
    class default
        success = .false.
        err_msg = "Arithmetic parser returned unsupported type of value"
    end select
end subroutine

integer(kind=4) function find_token_end(str) result(end_pos)
    character(len=*), intent(in) :: str

    integer(kind=4) :: ind_test, i, j

    !exclude leading blancks:
    do j = 1, len(str)
        if(.not. is_blanc(str(j:j))) exit
    end do
    if(j > len(str)) then
        end_pos = len(str)
        return
    end if

    end_pos = len(str(j:))+1
    do i=1, len(valid_token_end)
        ind_test = index(str(j:),valid_token_end(i:i))

        if(ind_test /= 0) end_pos = min(end_pos, ind_test)
    end do

    end_pos = end_pos + j-1
end

end module
