module arithmetic_parser_mod

use arithmetic_parser_tools_mod
use arithmetic_parser_val_mod
use str_util_mod,                only : alpha, digits, to_lower
use error_mod,                   only : error_t
use generic_config_mod,          only : generic_config_t

implicit none

contains

subroutine evaluate_expression(value, err, str, vars, global_vars, tmp_vars)
    type(error_t),    intent(out) :: err
    character(len=*), intent(in)  :: str
    class(generic_config_t), intent(inout), optional :: vars, global_vars, tmp_vars

    class(value_t),   intent(out), allocatable :: value

    integer(kind=4) :: pos, pos1
    integer(kind=4) :: mode ! 1 - unary op or value expected
                            ! 2 - binary op expected
    logical :: success
    class(operation_t), allocatable :: oper, oper1

    type(value_stack_t) :: values
    type(operation_stack_t) :: operations

    err = error_t(raised=.false.)

    pos  = 1
    mode = 1
    do while (pos <= len(str))

        if(index(" "//char(9)//new_line('c'),str(pos:pos)) /=0) then
            pos = pos+1
            cycle
        end if

        if(mode == 1) then

            if(str(pos:pos) == '(') then
                pos = pos+1
                call operations%push(openning_bracket_t())
                cycle
            end if

            call try_to_get_value(value,success,pos1,str(pos:))
            if(success) then
                pos = pos+pos1
                mode = 2
                call values%push(value)
                cycle
            end if

            call try_to_get_unary_op(oper,success,pos1,str(pos:))
            if(success) then
                pos = pos+pos1
                mode = 1
                call operations%push(oper)
                cycle
            end if

            call try_to_get_variable(value,success,pos1,str(pos:), &
                                     vars,global_vars,tmp_vars)
            if(success) then
                pos = pos+pos1
                mode = 2
                call values%push(value)
                cycle
            end if
        else

            if(str(pos:pos) == ")") then
                do while(operations%top_priority() > 0)
                    oper1 = operations%pop_operation()
                    value = oper1%apply(err,values)
                    if(err%raised) return

                    call values%push(value)
                    if(operations%is_empty()) then
                        call err%raise(message="syntax error, unbalanced parentheses"//str(pos:))
                        return
                    end if
                end do
                pos = pos+1 !shift ')'
                oper1 = operations%pop_operation() !pop oppening bracket
                cycle
            end if

            call try_to_get_binary_op(oper,success,pos1,str(pos:))
            if(success) then
                pos = pos+pos1
                mode = 1
                do while(.not. operations%is_empty())
                    if(operations%top_priority()<oper%priority()) exit
                    oper1 = operations%pop_operation()
                    value = oper1%apply(err,values)
                    if(err%raised) return
                    call values%push(value)
                end do
                call operations%push(oper)
                cycle
            end if
        end if

        call err%raise(message="Error, unexpected syntax: "// str(pos:))
        return
    end do

    do while(.not. operations%is_empty())
        oper1 = operations%pop_operation()
        value = oper1%apply(err,values)

        if(err%raised) return

        call values%push(value)
    end do

    value = values%pop_value()

end subroutine

subroutine try_to_get_value(val, success,end_pos,str)

    class(value_t), allocatable, intent(out) :: val
    logical,                     intent(out) :: success
    integer(kind=4),             intent(out) :: end_pos
    character(len=*),            intent(in)  :: str

    integer(kind=4) :: ival
    real(kind=8)    :: rval

    call try_to_get_float(rval, success, end_pos, str)
    if(success) then
        val = real_value_t(rval)
        return
    end if

    call try_to_get_int(ival, success, end_pos, str)
    if(success) then
        val = int_value_t(ival)
        return
    end if

end subroutine

subroutine try_to_get_varname(name, success,end_pos,str)

    character(len=:), allocatable, intent(out) :: name
    logical,                       intent(out) :: success
    integer(kind=4),               intent(out) :: end_pos
    character(len=*),              intent(in)  :: str

    integer(kind=4) :: i

    success = .false.
    if(index(alpha,str(1:1)) == 0 ) return

    success = .true.
    end_pos = 1

    do i = 2,len(str)
        if(index(alpha//digits//"_",str(i:i)) == 0) exit
        end_pos = i
    end do

    name = str(1:end_pos)
end subroutine

subroutine try_to_get_float(p, success,end_pos,str)

    real(kind=8),     intent(out) :: p
    logical,          intent(out) :: success
    integer(kind=4),  intent(out) :: end_pos
    character(len=*), intent(in)  :: str

    integer(kind=4) :: i, n_dots, n_exp, n_digits, exp_pos

    success = .false.
    n_exp = 0
    n_dots = 0
    n_digits = 0
    exp_pos = -123 !undef

    do i = 1, len(str)
        if(index("0123456789",str(i:i)) /=0) then
            n_digits = n_digits+1
        else if(i/=1.and.index("eEdD",str(i:i)) /= 0 .and. n_exp == 0) then
            n_exp = 1
            exp_pos = i
        else if(index("+-",str(i:i)) /= 0 .and. exp_pos == i-1) then
                !continue
        else if(str(i:i)=="." .and. n_exp == 0 .and. n_dots == 0) then
            n_dots = 1
        else
            exit
        end if
    end do

    if(n_digits > 0 .and. (n_dots == 1 .or. n_exp == 1)) then
        if(index(".0123456789",str(i-1:i-1)) == 0) return
        success = .true.
        end_pos = min(i-1,len(str))
        read(str(1:end_pos),*) p
    end if

end subroutine

subroutine try_to_get_int(ival, success,end_pos,str)

    integer(kind=4),  intent(out) :: ival
    logical,          intent(out) :: success
    integer(kind=4),  intent(out) :: end_pos
    character(len=*), intent(in)  :: str

    integer(kind=4) :: i, n_digits

    success = .false.
    n_digits = 0

    do i = 1, len(str)
        if(index("0123456789",str(i:i)) /=0) then
            n_digits = n_digits+1
        else
            exit
        end if
    end do


    if(n_digits > 0) then
        success = .true.
        end_pos = min(i-1,len(str))
        read(str(1:end_pos),*) ival
    end if

end subroutine

subroutine try_to_get_unary_op(oper,success,end_pos,str)

    class(operation_t), intent(out), allocatable :: oper
    logical,            intent(out) :: success
    integer(kind=4),    intent(out) :: end_pos
    character(len=*),   intent(in)  :: str

    character(len=:), allocatable :: name

    if(str(1:1) == "-") then
        end_pos = 1
        oper = unary_minus_t()
        success = .true.
        return
    end if

    call try_to_get_varname(name,success,end_pos,str)

    if(success) then
        name = to_lower(name)

        select case (name)
        case("sin","cos","tan","acos","asin","atan","exp","log","int","sqrt")
            oper = function_t(name)
        case default
            success = .false.
            end_pos = 0
            return
        end select

        success = .true.
        return
    end if

    success = .false.
    return
end subroutine

subroutine try_to_get_variable(val,success,end_pos,str,vars,global_vars,tmp_vars)

    class(value_t), allocatable, intent(out)    :: val
    logical,                     intent(out)    :: success
    integer(kind=4),             intent(out)    :: end_pos
    character(len=*),            intent(in)     :: str
    class(generic_config_t),     intent(inout), optional  :: vars, global_vars, tmp_vars

    character(len=:), allocatable :: name
    logical :: success_name

    call try_to_get_varname(name,success_name,end_pos,str)

    if(success_name) then

        success = .false.

        if(present(tmp_vars)) call get_var_value(val,success,name,tmp_vars)
        if(success) return

        if(present(vars)) call get_var_value(val,success,name,vars)
        if(success) return

        if(present(global_vars)) call get_var_value(val,success,name,global_vars)
        if(success) return
    end if

end subroutine

subroutine get_var_value(val, success, name, vars)

    class(value_t), allocatable, intent(out)    :: val
    logical,                     intent(out)    :: success
    character(len=*),            intent(in)     :: name
    class(generic_config_t),     intent(inout)  :: vars

    integer(kind=4) :: ival, ierr
    real(kind=8)    :: rval

    success = .false.

    call vars%get(ival,name,error_code=ierr)
    if(ierr == 0) then
        val = int_value_t(ival)
        success = .true.
        return
    end if

    call vars%get(rval,name,error_code=ierr)
    if(ierr == 0) then
        val = real_value_t(rval)
        success = .true.
        return
    end if

end subroutine

subroutine try_to_get_binary_op(oper,success,end_pos,str)

    class(operation_t), allocatable, intent(out) :: oper
    logical,          intent(out) :: success
    integer(kind=4),  intent(out) :: end_pos
    character(len=*), intent(in)  :: str

    if(str(1:1) == "-") then
        oper = subs_t()
        end_pos = 1
    else if(str(1:1) == "+") then
        oper = add_t()
        end_pos = 1
    else if(str(1:min(2,len(str))) == "**") then
        oper = pow_t()
        end_pos = 2
    else if(str(1:1) == "^") then
        oper = pow_t()
        end_pos = 1
    else if(str(1:1) == "*") then
        oper = mul_t()
        end_pos = 1
    else if(str(1:1) == "/") then
        oper = div_t()
        end_pos = 1
    else
        success = .false.
        return
    end if

    success = .true.
    return

end subroutine

end module
