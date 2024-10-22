module sbp_diff_mod

use grid_field_mod, only : grid_field_t, tile_field_t
use mesh_mod,       only : mesh_t, tile_mesh_t

implicit none

type, public, abstract :: sbp_diff_t
contains
    procedure, public :: apply
    procedure(apply_tile_i), public, deferred :: apply_tile
    procedure(apply_tile_i), public, deferred :: add_SAT_correction
end type sbp_diff_t

abstract interface
    subroutine apply_tile_i(this, f_out, f_in, mesh, direction)
        import :: sbp_diff_t, tile_field_t, tile_mesh_t
        class(sbp_diff_t),  intent(in)    :: this
        type(tile_field_t), intent(inout) :: f_out
        type(tile_field_t), intent(in)    :: f_in
        type(tile_mesh_t),  intent(in)    :: mesh
        character(len=*),   intent(in)    :: direction
    end subroutine apply_tile_i
end interface

contains

subroutine apply(this, f_out, f_in, mesh, direction)
    class(sbp_diff_t),  intent(in)    :: this
    type(grid_field_t), intent(inout) :: f_out
    type(grid_field_t), intent(in)    :: f_in
    character(len=*),   intent(in)    :: direction
    type(mesh_t),       intent(in)    :: mesh

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        call this%apply_tile(f_out%tile(t), f_in%tile(t), mesh%tile(t), direction)
    end do
end subroutine apply

end module sbp_diff_mod
