module cubsphere_scalar_field_mod

use scalar_field_native_mod, only : scalar_field_native_t

implicit none

type, extends(scalar_field_native_t) :: cubsphere_scalar_field_t
    contains
        procedure :: get_scalar_native
end type

contains

pure elemental real(kind=8) function get_scalar_native(this,alpha,beta,eta,panel_ind,time) result(f)

    class(cubsphere_scalar_field_t), intent(in) :: this
    real(kind=8),                    intent(in) :: alpha,beta,eta,time
    integer(kind=4),                 intent(in) :: panel_ind

    real(kind=8) :: ta, tb, sigm
    ta = tan(alpha);   tb = tan(beta)
    sigm = sqrt(1._8+ta**2+tb**2)

    f = panel_ind*(1._8+ta**2)**2*(1._8+tb**2)/sigm**4

end function

end module
