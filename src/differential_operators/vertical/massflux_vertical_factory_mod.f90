module massflux_vertical_factory_mod

use abstract_massflux_vertical_mod, only : massflux_vertical_t
use massflux_vertical_mod,          only : massflux_vert_up1_t, massflux_vert_up4_t, &
                                           massflux_vert_c2_t
use parcomm_mod,                    only : parcomm_global

implicit none

contains

subroutine create_vertical_massflux(massflux_vertical_op,massflux_vertical_op_name)

    class(massflux_vertical_t), allocatable, intent(out) :: massflux_vertical_op
    character(len=*), intent(in) :: massflux_vertical_op_name

    select case(massflux_vertical_op_name)
    case("up1")
        massflux_vertical_op = massflux_vert_up1_t()
    case("up4")
        massflux_vertical_op = massflux_vert_up4_t()
    case("c2")
        massflux_vertical_op = massflux_vert_c2_t()
    case default
        call parcomm_global%abort(__FILE__//": unknown vertical massflux operator name: "// massflux_vertical_op_name)
    end select

end subroutine

end module
