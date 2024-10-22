module grad_perp_factory_mod

use domain_mod,              only : domain_t
use abstract_grad_perp_mod,  only : grad_perp_operator_t
use parcomm_mod,             only : parcomm_global

implicit none

contains

subroutine create_grad_perp_operator(grad_perp, grad_perp_operator_name, domain)

    class(grad_perp_operator_t), allocatable, intent(out) :: grad_perp
    character(len=*),                         intent(in)  :: grad_perp_operator_name
    type(domain_t),                           intent(in)  :: domain

    select case(grad_perp_operator_name)
    case("grad_perp_c_sbp_sat_21", "grad_perp_c_sbp_sat_42", "grad_perp_c_sbp_sat_63")
        grad_perp = create_grad_perp_c_sbp_sat_operator(grad_perp_operator_name, domain)

    case default
        call parcomm_global%abort("unknown grad_perp operator name: "// grad_perp_operator_name)

    end select
end subroutine create_grad_perp_operator
function create_grad_perp_c_sbp_sat_operator(grad_perp_operator_name, domain) result(grad_perp)

    use grad_perp_sbp_sat_mod,  only : grad_perp_sbp_sat_t
    use exchange_factory_mod,   only : create_xy_points_halo_exchange
    use sbp_diff_i2c_21_mod,    only : sbp_diff_i2c_21_t
    use sbp_diff_i2c_42_mod,    only : sbp_diff_i2c_42_t
    use sbp_diff_i2c_63_mod,    only : sbp_diff_i2c_63_t

    character(len=*), intent(in)  :: grad_perp_operator_name
    type(domain_t),   intent(in)  :: domain
    type(grad_perp_sbp_sat_t)     :: grad_perp

    integer(kind=4) :: halo_width

    select case(grad_perp_operator_name)
    case("grad_perp_c_sbp_sat_21")
        halo_width = 2
        grad_perp%sbp_diff = sbp_diff_i2c_21_t()
    case("grad_perp_c_sbp_sat_42")
        halo_width = 4
        grad_perp%sbp_diff = sbp_diff_i2c_42_t()
    case("grad_perp_c_sbp_sat_63")
        halo_width = 6
        grad_perp%sbp_diff = sbp_diff_i2c_63_t()
    case default
        call parcomm_global%abort("grad_perp_factory_mod, create_grad_perp_c_sbp_sat_operator"// &
                                  " - unknown grad_perp operator: "//grad_perp_operator_name)
    end select

    call domain%get_mesh(grad_perp%mesh_u, "x" )
    call domain%get_mesh(grad_perp%mesh_v, "y" )

    grad_perp%exch_halo = create_xy_points_halo_exchange( &
                    domain%partition, domain%parcomm, domain%topology,  halo_width, 'full')
end function create_grad_perp_c_sbp_sat_operator

end module grad_perp_factory_mod
