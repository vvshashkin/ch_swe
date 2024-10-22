module massflux_3d_mod

use abstract_massflux_3d_mod,       only : massflux_3d_operator_t
use grid_field_mod,                 only : grid_field_t
use domain_mod,                     only : domain_t
use abstract_massflux_mod,          only : massflux_operator_t !2d horizontal massflux
use abstract_massflux_vertical_mod, only : massflux_vertical_t

implicit none


type, extends(massflux_3d_operator_t) :: massflux_horvert_t
    class(massflux_operator_t), allocatable :: massflux_hor
    class(massflux_vertical_t), allocatable :: massflux_vert
    contains
    procedure calc_massflux
end type

contains

subroutine calc_massflux(this, fx, fy, fz, f, u, v, w, domain)

    class(massflux_horvert_t), intent(inout) :: this
    type(domain_t),            intent(in)    :: domain
    type(grid_field_t),        intent(inout) :: f, u, v, w
    !output:
    type(grid_field_t),        intent(inout) :: fx, fy, fz

    call this%massflux_hor%calc_massflux(fx, fy, f, u, v, domain)
    call this%massflux_vert%calc_vertical_massflux(fz,f,w,domain%mesh_w)

end subroutine

end module
