module simple_latlon_scalar_field_mod

use scalar_field_latlon_mod, only : scalar_field_latlon_t

implicit none

type, extends(scalar_field_latlon_t) :: simple_latlon_scalar_field_t
    contains
        procedure :: get_scalar_latlon
end type

contains

pure elemental function get_scalar_latlon(this,phi,lambda,h,klev,time) result(f)

    class(simple_latlon_scalar_field_t), intent(in) :: this
    real(kind=8),                        intent(in) :: phi,lambda,h
    integer(kind=4),                     intent(in) :: klev
    real(kind=8),                        intent(in) :: time

    real(kind=8) :: f

    integer(kind=4) :: k1

    k1 = mod(klev-1,3)+1

    select case(k1)
    case(1)
        f = cos(phi)*cos(lambda)
    case(2)
        f = cos(phi)*sin(lambda)
    case(3)
        f = sin(phi)
    end select

end function

end module
