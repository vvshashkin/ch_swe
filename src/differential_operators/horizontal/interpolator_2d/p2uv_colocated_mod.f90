module p2uv_colocated_mod

use domain_mod,                    only : domain_t
use grid_field_mod,                only : grid_field_t
use abstract_interpolators2d_mod,  only : interpolator2d_scalar2vec_t

implicit none

!simply assigns u=p, v=p for colocated grids
type, public, extends(interpolator2d_scalar2vec_t) :: p2uv_colocated_t
contains
    procedure, public :: interp2d_scalar2vec
end type

contains

subroutine interp2d_scalar2vec(this, u, v, p, domain)

    class(p2uv_colocated_t), intent(inout) :: this
    type(grid_field_t),      intent(inout) :: p, u, v
    type(domain_t),          intent(in)    :: domain

    call u%assign(p,domain%mesh_p)
    call v%assign(p,domain%mesh_p)

end subroutine

end module
