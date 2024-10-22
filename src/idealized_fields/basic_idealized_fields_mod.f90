module basic_idealized_fields_mod

use idealized_field_mod, only : idealized_field_t
use grid_field_mod,      only : grid_field_t
use mesh_mod,            only : mesh_t

implicit none

type, extends(idealized_field_t) :: zero_idealized_field_t
contains
    procedure, public :: get_field
end type

contains

subroutine get_field(this, field, mesh, halo_width, fill_value, time)

    class(zero_idealized_field_t), intent(in)    :: this
    type(grid_field_t),            intent(inout) :: field
    type(mesh_t),                  intent(in)    :: mesh
    integer(kind=4),               intent(in), optional :: halo_width
    real(kind=8),                  intent(in), optional :: fill_value, time

    integer(kind=4) :: t, i, j, k, hw

    if(present(fill_value)) call field%assign(fill_value, mesh)

    hw = 0
    if(present(halo_width)) hw = halo_width

    do t = mesh%ts, mesh%te
        do k = mesh%tile(t)%ks, mesh%tile(t)%ke
            do j = mesh%tile(t)%js, mesh%tile(t)%je
                do i = mesh%tile(t)%is, mesh%tile(t)%ie
                    field%tile(t)%p(i,j,k) = 0.0_8
                end do
            end do
        end do
    end do

end subroutine

end module
