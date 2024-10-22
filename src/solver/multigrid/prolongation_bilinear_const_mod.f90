module prolongation_bilinear_const_mod

use abstract_prolongation_mod,   only : prolongation_operator_t
use domain_mod,                  only : domain_t
use grid_field_mod,              only : grid_field_t, tile_field_t
use mesh_mod,                    only : tile_mesh_t
use exchange_abstract_mod,       only : exchange_t
use vector_mod,                  only : abstract_vector_t
use grid_field_based_vector_mod, only : grid_field_based_vector_t

implicit none

! Combination of bilinear interpolation at the interior of the panels
! and piecewise constant interpolation (nearest neighbour) near the
! panels edges. Grids points are assumed to be cell centers (c-grid),
! so mesh_o is used in computations

type, public, extends(prolongation_operator_t) :: prolongation_bilinear_const_t
    class(exchange_t), allocatable :: halo_exchange
contains
    procedure, public :: apply
    procedure, public :: apply_to_grid_fields
end type prolongation_bilinear_const_t

contains

subroutine apply(this, f_fine, f_coarse, domain_fine, domain_coarse)

    class(prolongation_bilinear_const_t), intent(inout) :: this
    class(abstract_vector_t),             intent(inout) :: f_fine
    class(abstract_vector_t),             intent(inout) :: f_coarse
    type(domain_t),                       intent(in)    :: domain_fine
    type(domain_t),                       intent(in)    :: domain_coarse

    select type(f_coarse)
    class is (grid_field_based_vector_t)
        select type(f_fine)
        class is (grid_field_based_vector_t)

            call this%apply_to_grid_fields( f_fine%grid_field, &
                                            f_coarse%grid_field, &
                                            domain_fine, &
                                            domain_coarse)

        class default
            call domain_fine%parcomm%abort("Wrong type of f_fine!")
        end select
    class default
        call domain_coarse%parcomm%abort("Wrong type of f_coarse!")
    end select

end subroutine apply

subroutine apply_to_grid_fields(this, f_fine, f_coarse, domain_fine, domain_coarse)

    class(prolongation_bilinear_const_t), intent(inout) :: this
    type(grid_field_t),                   intent(inout) :: f_fine
    type(grid_field_t),                   intent(inout) :: f_coarse
    type(domain_t),                       intent(in)    :: domain_fine
    type(domain_t),                       intent(in)    :: domain_coarse

    integer(kind=4) :: t

    call this%halo_exchange%do(f_coarse, domain_coarse%parcomm)

    do t = domain_fine%partition%ts, domain_fine%partition%te

        call apply_tile(f_fine%tile(t),   domain_fine%mesh_o%tile(t),  &
                        f_coarse%tile(t), domain_coarse%mesh_o%tile(t))

    end do

end subroutine apply_to_grid_fields

subroutine apply_tile(f_f, mesh_f, f_c, mesh_c)

    type(tile_field_t), intent(inout) :: f_f
    type(tile_field_t), intent(in)    :: f_c
    type(tile_mesh_t),  intent(inout) :: mesh_f
    type(tile_mesh_t),  intent(in)    :: mesh_c

    real(kind=4)    :: coef_i, coef_j
    integer(kind=4) :: k, j, i, j_c, i_c

    do k = mesh_f%ks, mesh_f%ke
        do j = max(2, mesh_f%js), min(mesh_f%ny-1, mesh_f%je)
            j_c = j / 2
            coef_j = (1.0_8 + 2*mod(j+1, 2)) / 4.0_8
            do i = max(2, mesh_f%is), min(mesh_f%nx-1, mesh_f%ie)
                i_c = i / 2
                coef_i = (1.0_8 + 2*mod(i+1,2)) / 4.0_8
                f_f%p(i, j, k) = (      coef_i)*(      coef_j)*f_c%p(i_c,   j_c,   k) &
                               + (1.0_8-coef_i)*(      coef_j)*f_c%p(i_c+1, j_c,   k) &
                               + (      coef_i)*(1.0_8-coef_j)*f_c%p(i_c,   j_c+1, k) &
                               + (1.0_8-coef_i)*(1.0_8-coef_j)*f_c%p(i_c+1, j_c+1, k)
            end do
        end do

        if (mesh_f%is<=1 .and. mesh_f%ie>=1) then
            i   = 1
            i_c = 1
            do j = mesh_f%js, min(mesh_f%je, mesh_f%Ny-1)
                j_c = j/2+1
                f_f%p(i, j, k) = f_c%p(i_c, j_c, k)
            end do
        end if

        if (mesh_f%js<=mesh_f%ny .and. mesh_f%je >= mesh_f%ny) then
            j   = mesh_f%je
            j_c = mesh_c%je
            do i = mesh_f%is, min(mesh_f%ie, mesh_f%nx-1)
                i_c = i/2+1
                f_f%p(i, j, k) = f_c%p(i_c, j_c, k)
            end do
        end if

        if (mesh_f%is<=mesh_f%nx .and. mesh_f%ie>=mesh_f%nx) then
            i   = mesh_f%ie
            i_c = mesh_c%ie
            do j = max(2, mesh_f%js), mesh_f%je
                j_c = j/2
                f_f%p(i, j, k) = f_c%p(i_c, j_c, k)
            end do
        end if

        if (mesh_f%js<=1 .and. mesh_f%je>=1) then
            j   = mesh_f%js
            j_c = mesh_c%js
            do i = max(2, mesh_f%is), mesh_f%ie
                i_c = i/2
                f_f%p(i, j, k) = f_c%p(i_c, j_c, k)
            end do
        end if

    end do

end subroutine apply_tile

end module prolongation_bilinear_const_mod
