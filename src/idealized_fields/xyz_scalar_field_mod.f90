module xyz_scalar_field_mod

use scalar_field_cartesian_mod, only : scalar_field_cartesian_t

implicit none

type, extends(scalar_field_cartesian_t) :: xyz_scalar_field_t
    contains
        procedure :: get_scalar_cartesian
end type

contains

pure elemental real(kind=8) function get_scalar_cartesian(this,x,y,z,h,klev,time) &
                                                                        result(f)

    class(xyz_scalar_field_t), intent(in) :: this
    real(kind=8),              intent(in) :: x,y,z,h
    integer(kind=4),           intent(in) :: klev
    real(kind=8),              intent(in) :: time

    integer(kind=4) :: k1

    k1 = mod(klev-1,3)+1

    select case(k1)
    case(1)
        f = x
    case(2)
        f = y
    case(3)
        f=z
    end select

end function

end module
