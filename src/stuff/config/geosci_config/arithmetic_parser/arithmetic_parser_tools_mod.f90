module arithmetic_parser_tools_mod

use arithmetic_parser_val_mod
use error_mod, only : error_t

implicit none

type, abstract :: operation_t
    contains
        procedure(apply),    deferred :: apply
        procedure(priority), nopass, deferred :: priority
end type

type, abstract, extends(operation_t) :: unary_operation_t
    contains
        procedure, nopass :: priority => unary_priority
end type

type, private :: stack_frame_t
    class(*), allocatable :: data
    type(stack_frame_t), pointer :: previous, next
end type

type stack_t
    type(stack_frame_t), pointer :: last
    contains
        procedure :: pop
        procedure :: push
        procedure :: is_empty
end type

type, extends(stack_t) :: value_stack_t
    contains
        procedure :: push => push_value
        procedure :: pop_value
end type

type, extends(stack_t) :: operation_stack_t
    contains
        procedure :: push => push_operation
        procedure :: pop_operation
        procedure :: top_priority => get_top_operation_priority
end type

type, extends(operation_t) :: add_t
    contains
        procedure :: apply => add
        procedure, nopass :: priority => priority_add
end type

type, extends(operation_t) :: subs_t
    contains
        procedure :: apply => subs
        procedure, nopass :: priority => priority_subs
end type

type, extends(operation_t) :: mul_t
    contains
        procedure :: apply => mul
        procedure, nopass :: priority => priority_mul
end type

type, extends(operation_t) :: div_t
    contains
        procedure :: apply => div
        procedure, nopass :: priority => priority_div
end type

type, extends(operation_t) :: pow_t
    contains
        procedure :: apply => pow
        procedure, nopass :: priority => priority_pow
end type

type, extends(unary_operation_t) :: unary_minus_t
    contains
        procedure :: apply => unary_minus
end type

type, extends(unary_operation_t) :: function_t
    character(len=:), allocatable :: name
    contains
        procedure :: apply => calc_function
end type

type, extends(operation_t) :: openning_bracket_t
    contains
        procedure :: apply => do_nothing
        procedure, nopass :: priority => priority_bracket
end type

abstract interface

    function apply(this,err,values) result(res)
        import operation_t, value_stack_t, value_t, error_t
        class(operation_t),  intent(in)    :: this
        type(error_t),       intent(inout) :: err
        type(value_stack_t), intent(inout) :: values
        class(value_t), allocatable :: res
    end function
    integer(kind=4) function priority() result(pr)
    end function

end interface

contains

function pop(this) result(data)
    class(stack_t), intent(inout) :: this
    class(*),       allocatable   :: data

    type(stack_frame_t), pointer :: new_last

    data = this%last%data

    new_last => this%last%previous
    deallocate(this%last)
    this%last => new_last

end function

function pop_value(this) result(data)
    class(value_stack_t), intent(inout) :: this
    class(value_t),       allocatable   :: data

    class(*), allocatable :: item

    item = this%stack_t%pop()

    select type(item)
    class is (value_t)
        data = item
    end select

end function

function pop_operation(this) result(data)
    class(operation_stack_t), intent(inout) :: this
    class(operation_t),       allocatable   :: data

    class(*), allocatable :: item

    item = this%stack_t%pop()

    select type(item)
    class is (operation_t)
        data = item
    end select

end function

subroutine push(this,data)
    class(stack_t), intent(inout) :: this
    class(*),       intent(in)    :: data

    if(associated(this%last)) then
        allocate(this%last%next)
        this%last%next%previous => this%last
        this%last => this%last%next
    else
        allocate(this%last)
    end if
    this%last%data = data

end subroutine

subroutine push_value(this,data)
    class(value_stack_t), intent(inout) :: this
    class(*),             intent(in)    :: data

    select type (data)
    class is (value_t)
        call this%stack_t%push(data)
    type is (real(kind=8))
        call this%stack_t%push(real_value_t(data))
    type is (integer(kind=4))
        call this%stack_t%push(int_value_t(data))
    class default
        print *, "error"
    end select

end subroutine

subroutine push_operation(this,data)
    class(operation_stack_t), intent(inout) :: this
    class(*),                 intent(in)    :: data

    select type (data)
    class is (operation_t)
        call this%stack_t%push(data)
    class default
        print *, "error"
    end select

end subroutine

integer(kind=4) function get_top_operation_priority(this) result(pr)
    class(operation_stack_t), intent(inout) :: this

    class(*), pointer :: item

    item => this%last%data

    select type(item)
    class is (operation_t)
        pr = item%priority()
    end select

end function

logical function is_empty(this)
    class(stack_t), intent(inout) :: this
    is_empty = .not. associated(this%last)
end function

subroutine checkout_value(ival,rval,type_id,err,values)
    integer(kind=4),     intent(out)   :: ival
    real(kind=8),        intent(out)   :: rval
    integer(kind=4),     intent(out)   :: type_id
    type(error_t),       intent(inout) :: err
    type(value_stack_t), intent(inout) :: values
    
    class(value_t), allocatable :: value

    if(values%is_empty()) then
        call err%raise("Error: empty value stack!")
        return
    end if

    value = values%pop_value()
    select type (value)
    type is (int_value_t)
        ival = value%val
        type_id = INT_TYPE_ID
    type is (real_value_t)
        rval = value%val
        type_id = REAL_TYPE_ID
    class default
        type_id = UNKNOWN_TYPE_ID
    end select

end subroutine

integer(kind=4) function priority_add() result(pr)
    pr = 1
end function

function add(this,err,values) result(res)
    class(add_t),        intent(in)    :: this
    type(error_t),       intent(inout) :: err
    type(value_stack_t), intent(inout) :: values
    class(value_t), allocatable :: res

    class(value_t), allocatable :: val1, val2
    integer(kind=4) :: val1_i, val2_i, val1_type, val2_type
    real(kind=8)    :: val1_r, val2_r

    call checkout_value(val1_i, val1_r, val1_type, err, values)
    if(err%raised) return
    call checkout_value(val2_i, val2_r, val2_type, err, values)
    if(err%raised) return

    if(val1_type == INT_TYPE_ID .and. val2_type == INT_TYPE_ID) then
        res = int_value_t(val1_i+val2_i)
    else if(val1_type == REAL_TYPE_ID .and. val2_type == INT_TYPE_ID) then
        res = real_value_t(val1_r+val2_i)
    else if(val1_type == INT_TYPE_ID .and. val2_type == REAL_TYPE_ID) then
        res = real_value_t(val1_i+val2_r)
    else if(val1_type == REAL_TYPE_ID .and. val2_type == REAL_TYPE_ID) then
        res = real_value_t(val1_r+val2_r)
    else
        print *, "add type error"
    end if
end function

integer(kind=4) function priority_subs() result(pr)
    pr = 1
end function

function subs(this,err,values) result(res)
    class(subs_t),       intent(in)    :: this
    type(error_t),       intent(inout) :: err
    type(value_stack_t), intent(inout) :: values
    class(value_t), allocatable :: res

    class(value_t), allocatable :: val1, val2
    integer(kind=4) :: val1_i, val2_i, val1_type, val2_type
    real(kind=8)    :: val1_r, val2_r

    call checkout_value(val2_i, val2_r, val2_type, err, values)
    if(err%raised) return
    call checkout_value(val1_i, val1_r, val1_type, err, values)
    if(err%raised) return

    if(val1_type == INT_TYPE_ID .and. val2_type == INT_TYPE_ID) then
        res = int_value_t(val1_i-val2_i)
    else if(val1_type == REAL_TYPE_ID .and. val2_type == INT_TYPE_ID) then
        res = real_value_t(val1_r-val2_i)
    else if(val1_type == INT_TYPE_ID .and. val2_type == REAL_TYPE_ID) then
        res = real_value_t(val1_i-val2_r)
    else if(val1_type == REAL_TYPE_ID .and. val2_type == REAL_TYPE_ID) then
        res = real_value_t(val1_r-val2_r)
    else
        print *, "subs type error"
    end if
end function

integer(kind=4) function priority_mul() result(pr)
    pr = 2
end function

function mul(this,err,values) result(res)
    class(mul_t),        intent(in)    :: this
    type(error_t),       intent(inout) :: err
    type(value_stack_t), intent(inout) :: values
    class(value_t), allocatable :: res

    class(value_t), allocatable :: val1, val2
    integer(kind=4) :: val1_i, val2_i, val1_type, val2_type
    real(kind=8)    :: val1_r, val2_r

    call checkout_value(val1_i, val1_r, val1_type, err, values)
    if(err%raised) return
    call checkout_value(val2_i, val2_r, val2_type, err, values)
    if(err%raised) return

    if(val1_type == INT_TYPE_ID .and. val2_type == INT_TYPE_ID) then
        res = int_value_t(val1_i*val2_i)
    else if(val1_type == REAL_TYPE_ID .and. val2_type == INT_TYPE_ID) then
        res = real_value_t(val1_r*val2_i)
    else if(val1_type == INT_TYPE_ID .and. val2_type == REAL_TYPE_ID) then
        res = real_value_t(val1_i*val2_r)
    else if(val1_type == REAL_TYPE_ID .and. val2_type == REAL_TYPE_ID) then
        res = real_value_t(val1_r*val2_r)
    else
        print *, "mul type error"
    end if
end function

integer(kind=4) function priority_div() result(pr)
    pr = 2
end function

function div(this,err,values) result(res)
    class(div_t),        intent(in)    :: this
    type(error_t),       intent(inout) :: err
    type(value_stack_t), intent(inout) :: values
    class(value_t), allocatable :: res

    class(value_t), allocatable :: val1, val2
    integer(kind=4) :: val1_i, val2_i, val1_type, val2_type
    real(kind=8)    :: val1_r, val2_r

    call checkout_value(val2_i, val2_r, val2_type, err, values)
    if(err%raised) return
    call checkout_value(val1_i, val1_r, val1_type, err, values)
    if(err%raised) return

    if(val1_type == INT_TYPE_ID .and. val2_type == INT_TYPE_ID) then
        if(mod(val1_i,val2_i) == 0) then
            res = int_value_t(val1_i / val2_i)
        else
            res = real_value_t(real(val1_i,8) / real(val2_i,8))
        end if
    else if(val1_type == REAL_TYPE_ID .and. val2_type == INT_TYPE_ID) then
        res = real_value_t(val1_r / val2_i)
    else if(val1_type == INT_TYPE_ID .and. val2_type == REAL_TYPE_ID) then
        res = real_value_t(val1_i / val2_r)
    else if(val1_type == REAL_TYPE_ID .and. val2_type == REAL_TYPE_ID) then
        res = real_value_t(val1_r / val2_r)
    else
        print *, "mul type error"
    end if
end function

integer(kind=4) function priority_pow() result(pr)
    pr = 3
end function

function pow(this,err,values) result(res)
    class(pow_t),        intent(in)    :: this
    type(error_t),       intent(inout) :: err
    type(value_stack_t), intent(inout) :: values
    class(value_t), allocatable :: res

    class(value_t), allocatable :: val1, val2
    integer(kind=4) :: val1_i, val2_i, val1_type, val2_type
    real(kind=8)    :: val1_r, val2_r

    call checkout_value(val2_i, val2_r, val2_type, err, values)
    if(err%raised) return
    call checkout_value(val1_i, val1_r, val1_type, err, values)
    if(err%raised) return

    if(val1_type == INT_TYPE_ID .and. val2_type == INT_TYPE_ID) then
        res = int_value_t(val1_i**val2_i)
    else if(val1_type == REAL_TYPE_ID .and. val2_type == INT_TYPE_ID) then
        res = real_value_t(val1_r**val2_i)
    else if(val1_type == INT_TYPE_ID .and. val2_type == REAL_TYPE_ID) then
        res = real_value_t(val1_i**val2_r)
    else if(val1_type == REAL_TYPE_ID .and. val2_type == REAL_TYPE_ID) then
        res = real_value_t(val1_r**val2_r)
    else
        print *, "pow type error"
    end if
end function

integer(kind=4) function unary_priority() result(pr)
    pr = 10
end function

function unary_minus(this,err,values) result(res)
    class(unary_minus_t), intent(in)    :: this
    type(error_t),        intent(inout) :: err
    type(value_stack_t),  intent(inout) :: values
    class(value_t), allocatable :: res

    class(value_t), allocatable :: val1
    integer(kind=4) :: val1_i, val1_type
    real(kind=8)    :: val1_r

    call checkout_value(val1_i, val1_r, val1_type, err, values)
    if(err%raised) return

    if(val1_type == INT_TYPE_ID) then
        res = int_value_t(-val1_i)
    else if(val1_type == REAL_TYPE_ID) then
        res = real_value_t(-val1_r)
    else
        print *, "unary minus type error"
    end if
end function

function calc_function(this,err,values) result(res)
    class(function_t),   intent(in)    :: this
    type(error_t),       intent(inout) :: err
    type(value_stack_t), intent(inout) :: values
    class(value_t), allocatable :: res

    class(value_t), allocatable :: val1
    integer(kind=4) :: val1_i, val1_type
    real(kind=8)    :: val1_r

    call checkout_value(val1_i, val1_r, val1_type, err, values)
    if(err%raised) return

    if(val1_type /= INT_TYPE_ID .and. val1_type /= REAL_TYPE_ID) then
        print *, "calc function "//this%name//" type error"
        return
    end if

    if(val1_type == INT_TYPE_ID) then
        if(this%name == "int") then
            res = val1
            return
        end if
        val1_r = real(val1_i,8)
    end if

    select case(this%name)
    case("sin")
        res = real_value_t(sin(val1_r))
    case("cos")
        res = real_value_t(cos(val1_r))
    case("tan")
        res = real_value_t(tan(val1_r))
    case("asin")
        res = real_value_t(asin(val1_r))
    case("acos")
        res = real_value_t(acos(val1_r))
    case("atan")
        res = real_value_t(atan(val1_r))
    case("exp")
        res = real_value_t(exp(val1_r))
    case("log")
        res = real_value_t(log(val1_r))
    case("sqrt")
        res = real_value_t(sqrt(val1_r))
    case("int")
        res = int_value_t(val1_r)
    case default
        print *, "unknown function "// this%name
        return
    end select
end function

integer(kind=4) function priority_bracket() result(pr)
    pr = -10
end function

function do_nothing(this,err,values) result(res)
    class(openning_bracket_t), intent(in)    :: this
    type(error_t),             intent(inout) :: err
    type(value_stack_t),       intent(inout) :: values
    class(value_t), allocatable :: res

    call err%raise("unmatched openning bracket")

end function

end module
