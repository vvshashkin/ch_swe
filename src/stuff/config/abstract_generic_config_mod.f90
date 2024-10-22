module generic_config_mod

use string_mod,  only : string_t, is_in
use parcomm_mod, only : parcomm_global

implicit none

type, abstract :: generic_config_t
    contains
    procedure(parse), deferred :: parse
    procedure(get_var_names), deferred :: get_var_names

    procedure(get_i4),            deferred :: get_i4
    procedure(get_r8),            deferred :: get_r8
    procedure(get_string),        deferred :: get_string
    procedure(get_chars),         deferred :: get_chars
    procedure(get_logical),       deferred :: get_logical
    procedure(get_subconfig),     deferred :: get_subconfig
    procedure(get_i4_array),      deferred :: get_i4_array
    procedure(get_r8_array),      deferred :: get_r8_array
    procedure(get_string_array),  deferred :: get_string_array
    procedure(get_logical_array), deferred :: get_logical_array
    procedure                              :: get_subconfig_array
    procedure(get_empty_config),  deferred :: get_empty_subconfig

    generic :: get => get_i4, get_r8, get_logical, get_string, get_chars, &
                      get_subconfig,                                      &
                      get_i4_array, get_r8_array, get_string_array, get_logical_array, get_subconfig_array

    procedure :: check_no_unexpected

    procedure :: set_var, set_array
    generic   :: set => set_var, set_array
end type

abstract interface

    function get_var_names(this) result(var_names)
        import generic_config_t, string_t
        class(generic_config_t), intent(inout) :: this
        type(string_t), allocatable :: var_names(:)
    end function

    subroutine parse(this,str)
        import generic_config_t
        class(generic_config_t), intent(inout) :: this
        character(len=*),        intent(in)    :: str
    end subroutine

    function get_empty_config(this) result(empty)
        import generic_config_t
        class(generic_config_t), intent(in) :: this
        class(generic_config_t), allocatable :: empty
    end function

!get interfaces
#define VAR_TYPE integer(kind=4)
#define INTERFACE_NAME get_i4
#include "generic_config_get_interface.fpp"

#define VAR_TYPE real(kind=8)
#define INTERFACE_NAME get_r8
#include "generic_config_get_interface.fpp"

#define VAR_TYPE type(string_t)
#define INTERFACE_NAME get_string
#include "generic_config_get_interface.fpp"

#define VAR_TYPE character(len=:), allocatable
#define INTERFACE_NAME get_chars
#define CHAR
#include "generic_config_get_interface.fpp"
#undef CHAR

#define VAR_TYPE logical
#define INTERFACE_NAME get_logical
#include "generic_config_get_interface.fpp"

#define VAR_TYPE integer(kind=4)
#define INTERFACE_NAME get_i4_array
#include "generic_config_get_array_interface.fpp"

#define VAR_TYPE real(kind=8)
#define INTERFACE_NAME get_r8_array
#include "generic_config_get_array_interface.fpp"

#define VAR_TYPE logical
#define INTERFACE_NAME get_logical_array
#include "generic_config_get_array_interface.fpp"

#define VAR_TYPE type(string_t)
#define INTERFACE_NAME get_string_array
#include "generic_config_get_array_interface.fpp"

subroutine get_subconfig(this,q,var_name,default,found,error_code)
    import string_t, generic_config_t
    class(generic_config_t), intent(inout)         :: this
    character(len=*),        intent(in)            :: var_name
    logical,                 intent(out), optional :: found
    integer(kind=4),         intent(out), optional :: error_code

    class(generic_config_t), allocatable, intent(out)           :: q
    class(generic_config_t),              intent(in),  optional :: default
end subroutine

end interface

contains

logical function check_no_unexpected(this,expected_varnames, &
                                     unexpected_varnames) result(check)
    class(generic_config_t), intent(inout) :: this
    type(string_t),          intent(in)    :: expected_varnames(:)
    type(string_t), optional,intent(out), allocatable :: unexpected_varnames(:)

    type(string_t), allocatable :: var_names(:)
    integer(kind=4) :: i, n_unexpected
    integer(kind=4), allocatable :: i_unexpected(:)

    var_names = this%get_var_names()
    allocate(i_unexpected(size(var_names,1)))

    n_unexpected = 0
    check = .true.
    do i = 1, size(var_names,1)
        if(.not. is_in(var_names(i),expected_varnames)) then
            check = .false.
            n_unexpected = n_unexpected+1
            i_unexpected(n_unexpected) = i
        end if
    end do

    if(n_unexpected >0 .and. (present(unexpected_varnames))) then
        unexpected_varnames = var_names(i_unexpected(1:n_unexpected))
    end if
end function

subroutine set_var(this, var_name, q, overwrite, error_code)
    class(generic_config_t), intent(inout) :: this
    character(len=*),        intent(in)    :: var_name
    class(*),                intent(in)    :: q
    logical,                 intent(in),  optional :: overwrite
    integer(kind=4),         intent(out), optional :: error_code

    call parcomm_global%abort(__FILE__//": set_var subroutine is not implemented for this config")
end subroutine

subroutine set_array(this, var_name, q, overwrite, error_code)
    class(generic_config_t),   intent(inout) :: this
    character(len=*),          intent(in)    :: var_name
    class(*),                  intent(in)    :: q(:)
    logical,         optional, intent(in)    :: overwrite
    integer(kind=4), optional, intent(out)   :: error_code

    call parcomm_global%abort(__FILE__//": set_array subroutine is not implemented for this config")
end subroutine

subroutine get_subconfig_array(this,q,var_name,default,found,error_code)
    class(generic_config_t), intent(inout)            :: this
    class(generic_config_t), intent(out), allocatable :: q(:)
    character(len=*),        intent(in)               :: var_name
    class(generic_config_t), intent(in),  optional    :: default(:)
    logical,                 intent(out), optional    :: found
    integer(kind=4),         intent(out), optional    :: error_code

    call parcomm_global%abort(__FILE__//": get_subconfig_array subroutine is not implemented for this config")
end subroutine

end module generic_config_mod
