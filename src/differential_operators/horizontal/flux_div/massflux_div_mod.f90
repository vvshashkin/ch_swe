module massflux_div_mod

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
    
    type, extends(flux_div_operator_t) :: massflux_div_t
    
        class(div_operator_t),           allocatable :: div_op
        class(massflux_operator_t),      allocatable :: flux_op
        type(grid_field_t)                           :: hu, hv
    
        contains
        procedure :: calc_flux_div
    
    end type massflux_div_t

contains

subroutine calc_flux_div(this, flux_div, h, u, v, domain)

    !input:
    class(massflux_div_t), intent(inout) :: this
    type(grid_field_t),    intent(inout) :: h, u, v
    type(domain_t),        intent(in)    :: domain

    !output:
    type(grid_field_t),    intent(inout) :: flux_div

    call this%flux_op%calc_massflux(this%hu, this%hv, h, u, v, domain)
    call this%div_op%calc_div(flux_div, this%hu, this%hv, domain)

end subroutine

end module massflux_div_mod