module restriction_4point_cell_average_mod

use abstract_restriction_mod, only : restriction_operator_t
use domain_mod,               only : domain_t
use grid_field_mod,           only : grid_field_t, tile_field_t
use mesh_mod,                 only : tile_mesh_t
use exchange_abstract_mod,    only : exchange_t
use vector_mod,               only : abstract_vector_t

use grid_field_based_vector_mod, only : grid_field_based_vector_t

implicit none

! This type implements 4-point cell average restriction
! for the cell centred coarse and fine grids node locations (c-grid)

type, public, extends(restriction_operator_t) :: restriction_4point_cell_average_t
    class(exchange_t), allocatable :: halo_exchange
contains
    procedure, public :: apply
    procedure, public :: apply_to_grid_fields
end type restriction_4point_cell_average_t

contains

subroutine apply(this, f_coarse, f_fine, domain_coarse, domain_fine)

    class(restriction_4point_cell_average_t), intent(inout) :: this
    class(abstract_vector_t),                 intent(inout) :: f_coarse
    class(abstract_vector_t),                 intent(inout) :: f_fine
    type(domain_t),                           intent(in)    :: domain_coarse
    type(domain_t),                           intent(in)    :: domain_fine

    select type(f_coarse)
    class is (grid_field_based_vector_t)
        select type(f_fine)
        class is (grid_field_based_vector_t)

            call this%apply_to_grid_fields( f_coarse%grid_field, &
                                            f_fine%grid_field,   &
                                            domain_coarse,       &
                                            domain_fine)

        class default
            call domain_fine%parcomm%abort("Wrong type of f_fine!")
        end select
    class default
        call domain_coarse%parcomm%abort("Wrong type of f_coarse!")
    end select

end subroutine apply

subroutine apply_to_grid_fields(this, f_coarse, f_fine, domain_coarse, domain_fine)

    class(restriction_4point_cell_average_t), intent(inout) :: this
    type(grid_field_t),                       intent(inout) :: f_coarse
    type(grid_field_t),                       intent(inout) :: f_fine
    type(domain_t),                           intent(in)    :: domain_coarse
    type(domain_t),                           intent(in)    :: domain_fine

    integer(kind=4) :: t

    call this%halo_exchange%do(f_fine, domain_fine%parcomm)

    do t = domain_coarse%partition%ts, domain_coarse%partition%te

        call apply_tile(f_coarse%tile(t), domain_coarse%mesh_o%tile(t), &
                        f_fine%tile(t),   domain_fine%mesh_o%tile(t))

    end do

end subroutine apply_to_grid_fields

subroutine apply_tile(f_c, mesh_c, f_f, mesh_f)

    type(tile_field_t), intent(inout) :: f_c
    type(tile_field_t), intent(inout) :: f_f
    type(tile_mesh_t),  intent(in)    :: mesh_c
    type(tile_mesh_t),  intent(in)    :: mesh_f

    integer(kind=4) :: k, i, j, j_f, i_f
    real(kind=8)    :: volume_ratio

    volume_ratio = mesh_f%hx * mesh_f%hy / mesh_c%hx / mesh_c%hy

    do k = mesh_c%ks, mesh_c%ke
        do j = mesh_c%js, mesh_c%je
            j_f = 2*j
            do i = mesh_c%is, mesh_c%ie
                i_f = 2*i
                f_c%p(i,j,k) = ( &
                                 f_f%p(i_f,   j_f,   k)*mesh_f%J(i_f,   j_f,   k)  &
                               + f_f%p(i_f-1, j_f,   k)*mesh_f%J(i_f-1, j_f,   k)  &
                               + f_f%p(i_f,   j_f-1, k)*mesh_f%J(i_f,   j_f-1, k)  &
                               + f_f%p(i_f-1, j_f-1, k)*mesh_f%J(i_f-1, j_f-1, k)) &
                               * volume_ratio / mesh_c%J(i, j, k)
            end do
        end do
    end do

end subroutine apply_tile


end module restriction_4point_cell_average_mod
