module geosci_config_mod

use generic_config_mod,       only : generic_config_t
use string_mod,               only : string_t
use str_util_mod,             only : integer_to_str
use linked_list_mod,          only : linked_list_t, linked_list_node_t
use parcomm_mod,              only : parcomm_global
use geosci_config_parser_mod, only : parse_geosci_config

implicit none

type, extends(generic_config_t) :: geosci_config_t
    !WORKAROUND #safe_final
    type(linked_list_t) :: data = linked_list_t(safe_final=.false.)
contains
    procedure, public  :: parse
    procedure, public  :: get_var_names
    procedure, private :: get_i4
    procedure, private :: get_r8
    procedure, private :: get_logical
    procedure, private :: get_string
    procedure, private :: get_chars
    procedure, private :: get_subconfig
    procedure, private :: get_i4_array
    procedure, private :: get_r8_array
    procedure, private :: get_logical_array
    procedure, private :: get_string_array
    procedure, private :: get_subconfig_array
    procedure, public  :: get_empty_subconfig

    procedure, public  :: set_var, set_array
end type

integer(kind=4), parameter :: NO_ERROR = 0, NOT_FOUND_ERROR = 1, TYPE_ERROR = 2, &
                              IMPLICIT_OVERWRITE_ERROR = 3

contains

function get_empty_subconfig(this) result(empty)
    class(geosci_config_t), intent(in) :: this
    class(generic_config_t), allocatable :: empty

    empty = geosci_config_t()
end function

subroutine parse(this,str)
    class(geosci_config_t), intent(inout) :: this
    character(len=*),       intent(in)    :: str

    character(len=:), allocatable :: err_msg

    call parse_geosci_config(this, err_msg, str)

    if(allocated(err_msg)) then
        call parcomm_global%abort("geosci_config parse error: "//err_msg)
    end if

end subroutine

function get_var_names(this) result(var_names)
    class(geosci_config_t), intent(inout) :: this
    type(string_t), allocatable :: var_names(:)

    var_names = this%data%get_item_names()

end function


#define VAR_TYPE integer(kind=4)
#define SUBROUTINE_NAME get_i4
#include "geosci_config_get_var.fpp"

#define VAR_TYPE real(kind=8)
#define SUBROUTINE_NAME get_r8
#include "geosci_config_get_var.fpp"

#define VAR_TYPE logical
#define SUBROUTINE_NAME get_logical
#include "geosci_config_get_var.fpp"

#define VAR_TYPE type(string_t)
#define STRING_T
#define SUBROUTINE_NAME get_string
#include "geosci_config_get_var.fpp"
#undef STRING_T

#define VAR_TYPE character(len=:), allocatable
#define CHAR
#define SUBROUTINE_NAME get_chars
#include "geosci_config_get_var.fpp"
#undef CHAR

#define VAR_TYPE class(generic_config_t), allocatable
#define SUBCONFIG
#define SUBROUTINE_NAME get_subconfig
#include "geosci_config_get_var.fpp"
#undef SUBCONFIG

#define VAR_TYPE integer(kind=4)
#define SUBROUTINE_NAME get_i4_array
#include "geosci_config_get_array.fpp"

#define VAR_TYPE real(kind=8)
#define SUBROUTINE_NAME get_r8_array
#include "geosci_config_get_array.fpp"

#define VAR_TYPE logical
#define SUBROUTINE_NAME get_logical_array
#include "geosci_config_get_array.fpp"

#define VAR_TYPE type(string_t)
#define SUBROUTINE_NAME get_string_array
#include "geosci_config_get_array.fpp"

#define VAR_TYPE class(generic_config_t)
#define SUBCONFIG
#define SUBROUTINE_NAME get_subconfig_array
#include "geosci_config_get_array.fpp"
#undef SUBCONFIG

subroutine set_var(this,var_name,q,overwrite,error_code)
    class(geosci_config_t),  intent(inout) :: this
    character(len=*),        intent(in)    :: var_name
    class(*),                intent(in)    :: q
    logical,                 intent(in),  optional :: overwrite
    integer(kind=4),         intent(out), optional :: error_code

    class(*), pointer :: value
    logical :: overwrite_loc

    logical :: crash
    integer(kind=4) :: ind
    class(*), pointer :: subconfig

    ind = index(var_name,"%")
    if(ind /= 0) then

        call this%data%get_item(subconfig,var_name(1:ind-1))
        if(.not. associated(subconfig)) then
            call this%data%append(geosci_config_t(),var_name(1:ind-1))
            call this%data%get_item(subconfig,var_name(1:ind-1))
        end if

        select type (subconfig)
        class is (generic_config_t)
            call subconfig%set(var_name(ind+1:),q,overwrite,error_code)
        class default
            if(present(error_code)) then
                error_code = TYPE_ERROR
                return
            else
                call parcomm_global%abort(__FILE__//": error while trying to set "//var_name//", "//var_name(1:ind-1)//" is not subconfig")
            end if
        end select
        return
    end if

    overwrite_loc = .false.
    if(present(overwrite)) overwrite_loc = overwrite
    if(present(error_code)) error_code = NO_ERROR

    call this%data%get_item(value,var_name)

    if(.not. associated(value)) then
        call this%data%append(q,var_name)
    else if(overwrite_loc) then
        call this%data%change_item(q,var_name)
    else if(present(error_code)) then
        error_code = IMPLICIT_OVERWRITE_ERROR
    else
        call parcomm_global%abort(__FILE__//": trying to redefine existing variable, need explicit overwrite flag to be true")
    end if
end subroutine

subroutine set_array(this, var_name, q, overwrite, error_code)
    class(geosci_config_t),    intent(inout) :: this
    character(len=*),          intent(in)    :: var_name
    class(*),                  intent(in)    :: q(:)
    logical,         optional, intent(in)    :: overwrite
    integer(kind=4), optional, intent(out)   :: error_code

    integer(kind=4) :: i
    class(generic_config_t), allocatable :: array_conf

    array_conf = this%get_empty_subconfig()
    call this%set(var_name, array_conf, overwrite, error_code)
    if(present(error_code)) then
        if(error_code /= 0) return
    end if

    call this%set(var_name//"%length",size(q,1),overwrite,error_code)
    if(present(error_code)) then
        if(error_code /= 0) return
    end if

    do i = 1, size(q,1)
        call this%set(var_name//"%"//integer_to_str(i),q(i),overwrite,error_code)
        if(present(error_code)) then
            if(error_code /= 0) return
        end if
    end do
end subroutine

end module
