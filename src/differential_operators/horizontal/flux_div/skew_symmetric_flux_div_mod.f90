module skew_symmetric_flux_div_mod

use domain_mod,                   only : domain_t
use grid_field_mod,               only : grid_field_t
use abstract_div_mod,             only : div_operator_t
use abstract_grad_mod,            only : grad_operator_t
use abstract_massflux_mod,        only : massflux_operator_t
use halo_mod,                     only : halo_t
use parcomm_mod,                  only : parcomm_global
use vec_math_mod,                 only : multiply_by_J_self, divide_by_J_self
use abstract_flux_div_mod,        only : flux_div_operator_t
use abstract_interpolators2d_mod, only : interpolator2d_vec2vec_t

implicit none

type, extends(flux_div_operator_t) :: skew_symmetric_flux_div_t

    class(div_operator_t),           allocatable :: div_op
    class(grad_operator_t),          allocatable :: grad_op
    class(massflux_operator_t),      allocatable :: flux_op
    class(interpolator2d_vec2vec_t), allocatable :: interp_uv2p_op
    class(halo_t),                   allocatable :: result_projection
    type(grid_field_t)                           :: hu, hv, grad_x, grad_y, div, up, vp

    contains
    procedure :: calc_flux_div

end type skew_symmetric_flux_div_t

contains

subroutine calc_flux_div(this, flux_div, h, u, v, domain)

    !input:
    class(skew_symmetric_flux_div_t), intent(inout) :: this
    type(grid_field_t),               intent(inout) :: h, u, v
    type(domain_t),                   intent(in)    :: domain

    !output:
    type(grid_field_t),               intent(inout) :: flux_div

    call this%flux_op%calc_massflux(this%hu, this%hv, h, u, v, domain)
    call this%div_op%calc_div(flux_div, this%hu, this%hv, domain)
    call flux_div%assign(0.5_8, flux_div, domain%mesh_p)

    call this%grad_op%calc_grad(this%grad_x, this%grad_y, h, domain)
    call this%grad_x%assign_prod(1.0_8, this%grad_x, u, domain%mesh_u)
    call this%grad_y%assign_prod(1.0_8, this%grad_y, v, domain%mesh_v)

    call multiply_by_J_self(this%grad_x, domain%mesh_u)
    call multiply_by_J_self(this%grad_y, domain%mesh_v)

    call this%interp_uv2p_op%interp2d_vec2vec(this%up, this%vp, this%grad_x, this%grad_y, domain)

    call this%up%update(1.0_8, this%vp, domain%mesh_p)
    call divide_by_J_self(this%up, domain%mesh_p)

    call flux_div%update(0.5_8, this%up, domain%mesh_p)

    call this%div_op%calc_div(this%div, u, v, domain)
    call this%div%assign_prod(1.0_8, this%div, h, domain%mesh_p)

    call flux_div%update(0.5_8, this%div, domain%mesh_p)

    !Enforce continuity etc if needed:
    if(allocated(this%result_projection)) &
        call this%result_projection%get_halo_scalar(flux_div, domain, 0) 

end subroutine

end module skew_symmetric_flux_div_mod