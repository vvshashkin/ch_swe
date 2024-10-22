module laplace_true_hor_factory_mod

use abstract_laplace_mod,     only : laplace_operator_t
use laplace_true_hor_mod,     only : laplace_true_hor_t
use domain_mod,               only : domain_t
use grid_field_factory_mod,   only : create_grid_field
use grad_3d_factory_mod,      only : create_grad_3d_operator
use div_3d_factory_mod,       only : create_div_3d_operator
use co2contra_3d_factory_mod, only : create_co2contra_3d_operator
use parcomm_mod,              only : parcomm_global

implicit none

contains

subroutine create_laplace_true_hor(laplace_operator, laplace_operator_name, domain)
    character(len=*), intent(in) :: laplace_operator_name
    type(domain_t),   intent(in) :: domain
    !output:
    class(laplace_operator_t), allocatable, intent(out) :: laplace_operator

    type(laplace_true_hor_t), allocatable :: laplace
    character(len=:), allocatable :: hor_grad_name, ver_grad_name, hor_div_name, &
                                     ver_div_name, co2contra_name
    integer(kind=4) :: halo_width, t, i, j, k

    allocate(laplace)

    select case(laplace_operator_name)
    case("laplace_3d_true_hor_c_sbp21_vsbp21")
        hor_grad_name  = "grad_c_sbp_sat_21"
        ver_grad_name  = "eta_diff_p2w_sbp21"
        hor_div_name   = "div_c_sbp_sat_21"
        ver_div_name   = "eta_diff_w2p_sbp21"
        co2contra_name = "co2contra_3d_Cgrid_h_sbp21_v_sbp21"
        halo_width = 2
    case("laplace_3d_true_hor_c_sbp42_vsbp21")
        hor_grad_name  = "grad_c_sbp_sat_42"
        ver_grad_name  = "eta_diff_p2w_sbp21"
        hor_div_name   = "div_c_sbp_sat_42"
        ver_div_name   = "eta_diff_w2p_sbp21"
        co2contra_name = "co2contra_3d_Cgrid_h_sbp42_v_sbp21"
        halo_width = 4
    case("laplace_3d_true_hor_c_sbp21_vsbp42")
        hor_grad_name  = "grad_c_sbp_sat_21"
        ver_grad_name  = "eta_diff_p2w_sbp42"
        hor_div_name   = "div_c_sbp_sat_21"
        ver_div_name   = "eta_diff_w2p_sbp42"
        co2contra_name = "co2contra_3d_Cgrid_h_sbp21_v_sbp42"
        halo_width = 2
    case("laplace_3d_true_hor_c_sbp42_vsbp42")
        hor_grad_name  = "grad_c_sbp_sat_42"
        ver_grad_name  = "eta_diff_p2w_sbp42"
        hor_div_name   = "div_c_sbp_sat_42"
        ver_div_name   = "eta_diff_w2p_sbp42"
        co2contra_name = "co2contra_3d_Cgrid_h_sbp42_v_sbp42"
        halo_width = 4
    case default
        call parcomm_global%abort(__FILE__//" unknown laplace_true_hor operator_name: "// laplace_operator_name)
    end select

    call create_grad_3d_operator(laplace%grad,domain,hor_grad_name,ver_grad_name)
    call create_div_3d_operator(laplace%div,domain,hor_div_name,ver_div_name)
    call create_co2contra_3d_operator(laplace%co2contra,domain,co2contra_name)

    call create_grid_field(laplace%fx, halo_width,0,domain%mesh_x)
    call create_grid_field(laplace%fy, halo_width,0,domain%mesh_y)
    call create_grid_field(laplace%fz, halo_width,0,domain%mesh_z)
    call create_grid_field(laplace%fxt,halo_width,0,domain%mesh_x)
    call create_grid_field(laplace%fyt,halo_width,0,domain%mesh_y)
    call create_grid_field(laplace%fzt,halo_width,0,domain%mesh_z)
    call create_grid_field(laplace%hor_factor,0,0,domain%mesh_w)

    do t = domain%mesh_w%ts, domain%mesh_w%te
        do k = domain%mesh_w%tile(t)%ks, domain%mesh_w%tile(t)%ke
            do j = domain%mesh_w%tile(t)%js, domain%mesh_w%tile(t)%je
                do i = domain%mesh_w%tile(t)%is, domain%mesh_w%tile(t)%ie
                    laplace%hor_factor%tile(t)%p(i,j,k) = 1.0_8/domain%mesh_w%tile(t)%a3(4,i,j,k)**2
                end do
            end do
        end do
    end do

    call move_alloc(laplace,laplace_operator)
end subroutine

end module
