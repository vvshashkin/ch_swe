module abstract_dp_interpolator_mod

use mesh_mod,       only : tile_mesh_t
use grid_field_mod, only : grid_field_t, tile_field_t
use domain_mod,     only : domain_t

implicit none

type, abstract :: dp_interpolator_t
    integer(kind=4) :: min_alloc_size = 8, extension_overhead = 4
    contains
        procedure(init_i),         deferred :: init
        procedure(calc_weights_i), deferred :: calc_weights
        procedure(interpolate_i),  deferred :: interpolate
        procedure :: ext_halo
end type

type :: dp_interpolator_container_t !for arrays of various interpolators
    class(dp_interpolator_t), allocatable :: interpolator
end type

abstract interface

    subroutine init_i(this, mesh_name, domain)
        import dp_interpolator_t, domain_t
        class(dp_interpolator_t), intent(inout) :: this
        character(len=*),         intent(in)    :: mesh_name
        type(domain_t),           intent(in)    :: domain
    end subroutine

    subroutine calc_weights_i(this, alpha, beta, eta, tile_ind)
        import dp_interpolator_t, tile_mesh_t
        class(dp_interpolator_t), intent(inout) :: this
        real(kind=8),             intent(in)    :: alpha(:), beta(:), eta(:)
        integer(kind=4),          intent(in)    :: tile_ind
    end subroutine

    subroutine interpolate_i(this, qout, qin)
        import dp_interpolator_t, tile_field_t
        class(dp_interpolator_t), intent(in)    :: this
        real(kind=8),             intent(inout) :: qout(:)
        type(tile_field_t),       intent(in)    :: qin
    end subroutine
end interface

contains

subroutine ext_halo(this,q,domain)
    class(dp_interpolator_t), intent(inout) :: this
    type(grid_field_t),       intent(inout) :: q
    type(domain_t),           intent(in)    :: domain

    !just do nothing by default
    !if specific interpolation need some halo-procedure, this subroutine
    !can be redefined
end subroutine

end module