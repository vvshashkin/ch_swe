module config_tools_mod

use generic_config_mod,          only : generic_config_t
use geosci_config_mod,           only : geosci_config_t
use namelist_read_mod,           only : read_namelist_as_str
use parcomm_mod,                 only : parcomm_t

implicit none

contains

function parse_config_file(file_name,parcomm) result(config)

    character(len=*),         intent(in) :: file_name
    type(parcomm_t),          intent(in) :: parcomm
    class(generic_config_t), allocatable :: config

    character(len=:), allocatable :: config_str

    call read_namelist_as_str(config_str,file_name,parcomm%myid)

    config = get_config_type_by_extension(file_name)

    call config%parse(config_str)

end function

function get_config_type_by_extension(file_name) result(config)

    character(len=*), intent(in) :: file_name

    class(generic_config_t), allocatable :: config

    integer(kind=4) :: l

    l = len(file_name)

    if(file_name(l-3:l) == ".cfg") then
        config = geosci_config_t()
    else
        config = geosci_config_t()
    end if

end function

end module
