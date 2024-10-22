module scalar_field_latlon_mod

use idealized_field_mod, only : idealized_field_t
use grid_field_mod,      only : grid_field_t
use mesh_mod,            only : mesh_t
use sph_coords_mod,      only : cart2sph

implicit none

type, abstract, extends(idealized_field_t) :: scalar_field_latlon_t
    contains
        procedure, public :: get_field
        procedure(get_scalar_latlon), deferred :: get_scalar_latlon
end type

abstract interface

    pure elemental function get_scalar_latlon(this,phi,lambda,h,klev,time) result(f)

        import scalar_field_latlon_t

        class(scalar_field_latlon_t), intent(in) :: this
        real(kind=8),                 intent(in) :: phi,lambda,h
        integer(kind=4),              intent(in) :: klev
        real(kind=8),                 intent(in) :: time

        real(kind=8) :: f

    end function

end interface

contains

subroutine get_field(this, field, mesh, halo_width, fill_value, time)

    class(scalar_field_latlon_t), intent(in)    :: this
    type(grid_field_t),              intent(inout) :: field
    type(mesh_t),                    intent(in)    :: mesh
    integer(kind=4),                 intent(in), optional :: halo_width
    real(kind=8),                    intent(in), optional :: fill_value, time

    integer(kind=4) :: t, i, j, k, hw, isv, iev, jsv, jev, ksv, kev
    real(kind=8)    :: phi, lambda, time_loc

    hw = 0
    if(present(halo_width)) hw = halo_width

    time_loc = this%time_default
    if(present(time)) time_loc = time

    do t = mesh%ts, mesh%te

        if(present(fill_value)) then
            isv = field%tile(t)%is; iev = field%tile(t)%ie
            jsv = field%tile(t)%js; jev = field%tile(t)%je
            ksv = field%tile(t)%ks; kev = field%tile(t)%ke
            field%tile(t)%p(isv:iev,jsv:jev,ksv:kev) = fill_value
        end if

        do k = mesh%tile(t)%ks, mesh%tile(t)%ke
            do j = mesh%tile(t)%js-hw, mesh%tile(t)%je+hw
                do i = mesh%tile(t)%is-hw, mesh%tile(t)%ie+hw
                    call cart2sph(mesh%tile(t)%rx(i,j,k),mesh%tile(t)%ry(i,j,k), &
                                  mesh%tile(t)%rz(i,j,k),lambda,phi)
                    field%tile(t)%p(i,j,k) = &
                        this%get_scalar_latlon(phi,lambda,mesh%tile(t)%h(i,j,k),k,time_loc)
                end do
            end do
        end do
    end do

end subroutine

end module
