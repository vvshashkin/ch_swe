module config_orography_mod

use config_mod, only : config_t

implicit none

type, extends(config_t) :: config_test_orography_t
    real(kind=8) :: h = 1.0_8
    contains
        procedure :: parse
end type

contains

subroutine parse(this, config_string)
    class(config_test_orography_t), intent(inout) :: this
    character(len=*),               intent(in)    :: config_string

    real(kind=8) :: h
    namelist /orography/ h

    read(config_string,orography)

    this%h = h

end subroutine parse

end module config_orography_mod
