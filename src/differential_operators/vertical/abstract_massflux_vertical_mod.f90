module abstract_massflux_vertical_mod

use grid_field_mod, only : grid_field_t, tile_field_t
use mesh_mod,       only : mesh_t, tile_mesh_t

implicit none

type, abstract :: massflux_vertical_t
    contains
    procedure calc_vertical_massflux
    procedure(calc_vertical_massflux_tile), deferred :: calc_vertical_massflux_tile
end type

abstract interface
    subroutine calc_vertical_massflux_tile(this, massflux, f, eta_dot, mesh)
        import massflux_vertical_t, tile_field_t, tile_mesh_t
        class(massflux_vertical_t),  intent(in)    :: this
        type(tile_field_t),          intent(in)    :: f, eta_dot
        type(tile_mesh_t),           intent(in)    :: mesh
        !output
        type(tile_field_t),          intent(inout) :: massflux
    end subroutine calc_vertical_massflux_tile
end interface

contains

subroutine calc_vertical_massflux(this, massflux, f, eta_dot, mesh)
    class(massflux_vertical_t), intent(in)    :: this
    type(grid_field_t),         intent(in)    :: f, eta_dot
    type(mesh_t),               intent(in)    :: mesh
    !output
    type(grid_field_t),         intent(inout) :: massflux

    integer(kind=4) :: t

    do t=mesh%ts, mesh%te
        call this%calc_vertical_massflux_tile(massflux%tile(t), f%tile(t), eta_dot%tile(t), mesh%tile(t))
    end do
end subroutine calc_vertical_massflux

end module
