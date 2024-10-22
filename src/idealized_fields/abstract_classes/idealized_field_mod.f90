module idealized_field_mod

!abstract interface for setting field, either scalar or vector component

use grid_field_mod, only : grid_field_t
use mesh_mod,       only : mesh_t

implicit none

type, abstract :: idealized_field_t
    real(kind=8) :: time_default = 0.0_8
    contains
        procedure(get_field), public, deferred :: get_field
end type

abstract interface
    subroutine get_field(this, field, mesh, halo_width, fill_value, time)

        import idealized_field_t, grid_field_t, mesh_t

        class(idealized_field_t), intent(in)    :: this
        type(grid_field_t),       intent(inout) :: field
        type(mesh_t),             intent(in)    :: mesh
        integer(kind=4),          intent(in), optional :: halo_width
        real(kind=8),             intent(in), optional :: fill_value, time

    end subroutine
end interface

end module
