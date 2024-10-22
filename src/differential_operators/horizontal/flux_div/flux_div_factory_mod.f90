module flux_div_factory_mod

use abstract_flux_div_mod,       only : flux_div_operator_t
use skew_symmetric_flux_div_mod, only : skew_symmetric_flux_div_t
use massflux_div_mod,            only : massflux_div_t

use grad_factory_mod,            only : create_grad_operator
use div_factory_mod,             only : create_div_operator
use massflux_factory_mod,        only : create_massflux_operator
use interpolator2d_factory_mod,  only : create_vec2vec_interpolator2d
use halo_factory_mod,            only : create_halo_procedure

use grid_field_factory_mod,      only : create_grid_field
use domain_mod,                  only : domain_t

implicit none

contains

subroutine create_flux_div_operator(flux_div_op, flux_div_op_name, domain)

    class(flux_div_operator_t), allocatable, intent(out) :: flux_div_op
    character(len=*),                        intent(in)  :: flux_div_op_name
    type(domain_t),                          intent(in)  :: domain

    select case(flux_div_op_name)
    case("skew_symmetric_flux_div_Ch_sbp21")
        call create_skew_symmetric_flux_div(flux_div_op, "div_ch_sbp_sat_21",           &
                                            "grad_ch_sbp_sat_21", "massflux_ch_sbp21", &
                                            "interp2d_uv2pvec_Ch_sbp21", domain)
    case("skew_symmetric_flux_div_Ch_sbp21_proj")
        call create_skew_symmetric_flux_div(flux_div_op, "div_ch_sbp_sat_21_proj",           &
                                            "grad_ch_sbp_sat_21", "massflux_ch_sbp21", &
                                            "interp2d_uv2pvec_Ch_sbp21", domain, "Ah_scalar_sync")
    case("skew_symmetric_flux_div_Ch_sbp42")
        call create_skew_symmetric_flux_div(flux_div_op, "div_ch_sbp_sat_42",           &
                                            "grad_ch_sbp_sat_42", "massflux_ch_sbp42", &
                                            "interp2d_uv2pvec_Ch_sbp42", domain)
    case("skew_symmetric_flux_div_Ch_sbp42_proj")
        call create_skew_symmetric_flux_div(flux_div_op, "div_ch_sbp_sat_42_proj",           &
                                            "grad_ch_sbp_sat_42", "massflux_ch_sbp42", &
                                            "interp2d_uv2pvec_Ch_sbp42", domain, "Ah_scalar_sync")
    case("skew_symmetric_flux_div_C_sbp21")
        call create_skew_symmetric_flux_div(flux_div_op, "div_c_sbp_sat_21",           &
                                            "grad_c_sbp_sat_21", "massflux_c_sbp21", &
                                            "interp2d_uv2pvec_C_sbp21", domain)
    case("skew_symmetric_flux_div_C_sbp42")
        call create_skew_symmetric_flux_div(flux_div_op, "div_c_sbp_sat_42",           &
                                            "grad_c_sbp_sat_42", "massflux_c_sbp42", &
                                            "interp2d_uv2pvec_C_sbp42", domain)

    case("massflux_div_C_sbp21")
        call create_massflux_div(flux_div_op, "div_c_sbp_sat_21", "massflux_c_sbp21", domain)
    case("massflux_div_C_sbp42")
        call create_massflux_div(flux_div_op, "div_c_sbp_sat_42", "massflux_c_sbp42", domain)
    case default
        call domain%parcomm%abort("create_flux_div_operator, error - unknown flux_div operator name: "// flux_div_op_name)
    end select

end subroutine

subroutine create_massflux_div(flux_div_op, divergence_name, massflux_name, domain)

    class(flux_div_operator_t), allocatable, intent(out) :: flux_div_op

    character(len=*), intent(in) :: divergence_name, massflux_name
    type(domain_t),   intent(in) :: domain

    type(massflux_div_t), allocatable :: massflux_div_op
    !WORKAROUND:
    integer(kind=4), parameter :: halo_width = 5

    allocate(massflux_div_op)

    massflux_div_op%div_op  = create_div_operator(domain, divergence_name)
    massflux_div_op%flux_op = create_massflux_operator(domain, massflux_name)

    call create_grid_field(massflux_div_op%hu, halo_width, 0, domain%mesh_u)
    call create_grid_field(massflux_div_op%hv, halo_width, 0, domain%mesh_v)
   
    call move_alloc(massflux_div_op, flux_div_op)

end subroutine

subroutine create_skew_symmetric_flux_div(flux_div_op, divergence_name, gradient_name, &
                                          massflux_name, uv2p_name, domain, result_proj_name)

    class(flux_div_operator_t), allocatable, intent(out) :: flux_div_op

    character(len=*), intent(in) :: divergence_name, gradient_name, massflux_name, uv2p_name
    type(domain_t),   intent(in) :: domain
    character(len=*), intent(in), optional :: result_proj_name

    type(skew_symmetric_flux_div_t), allocatable :: skew_flux_div_op
    !WORKAROUND:
    integer(kind=4), parameter :: halo_width = 5

    allocate(skew_flux_div_op)

    skew_flux_div_op%grad_op = create_grad_operator(domain, gradient_name)
    skew_flux_div_op%div_op  = create_div_operator(domain, divergence_name)
    skew_flux_div_op%flux_op = create_massflux_operator(domain, massflux_name)
    call create_vec2vec_interpolator2d(skew_flux_div_op%interp_uv2p_op, uv2p_name, domain)

    if(present(result_proj_name)) &
        call create_halo_procedure(skew_flux_div_op%result_projection, domain, 0, result_proj_name)

    call create_grid_field(skew_flux_div_op%hu, halo_width, 0, domain%mesh_u)
    call create_grid_field(skew_flux_div_op%hv, halo_width, 0, domain%mesh_v)
    call create_grid_field(skew_flux_div_op%grad_x, halo_width, 0, domain%mesh_u)
    call create_grid_field(skew_flux_div_op%grad_y, halo_width, 0, domain%mesh_v)
    call create_grid_field(skew_flux_div_op%div, halo_width, 0, domain%mesh_p)
    call create_grid_field(skew_flux_div_op%up, halo_width, 0, domain%mesh_u)
    call create_grid_field(skew_flux_div_op%vp, halo_width, 0, domain%mesh_v)

    call move_alloc(skew_flux_div_op, flux_div_op)

end subroutine

end module flux_div_factory_mod