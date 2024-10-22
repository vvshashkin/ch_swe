module grad_factory_mod

use abstract_grad_mod, only : grad_operator_t
use domain_mod,        only : domain_t
use parcomm_mod,       only : parcomm_global

implicit none

contains

function create_grad_operator(domain, grad_operator_name) result(grad)
    type(domain_t),    intent(in)  :: domain
    character(len=*),  intent(in)  :: grad_operator_name

    class(grad_operator_t), allocatable :: grad

    select case(grad_operator_name)
    case ("grad_c_sbp_sat_21", "grad_c_sbp_sat_42", "grad_c_sbp_sat_63")
        grad =  create_grad_c_sbp_sat_operator(grad_operator_name, domain, is_z = .false.)

    case ("grad_c_sbp_sat_21_z", "grad_c_sbp_sat_42_z", "grad_c_sbp_sat_63_z")
        grad = create_grad_c_sbp_sat_operator(grad_operator_name, domain, is_z = .true.)

    case ("grad_ch_sbp_sat_21", "grad_ch_sbp_sat_42", "grad_ch_sbp_sat_63")
        grad =  create_grad_ch_sbp_sat_operator(grad_operator_name, domain, is_z = .false.)

    case("grad_ah_sbp_proj_21", "grad_ah_sbp_proj_42", "grad_ah_sbp_proj_43", "grad_ah_sbp_proj_63")
            grad = create_grad_ah_sbp_proj_operator(grad_operator_name, domain, is_z = .false.)

    case("gradient_c2_ecs")
        grad = create_grad_c2_ecs_operator(domain)

    case("gradient_a2_ecs", "gradient_a2_cons")
        grad = create_grad_a2_operator(domain,grad_operator_name)

    case("gradient_ch_ecs_halo2", "gradient_ch_ecs_halo4")
        grad = create_grad_ch_halo_operator(domain, grad_operator_name)

    case default
        call parcomm_global%abort("unknown gradient operator: "//grad_operator_name)
    end select

end function create_grad_operator

function create_grad_c2_ecs_operator(domain) result(grad)

    use grad_c2_ecs_mod, only : grad_c2_ecs_t
    use halo_factory_mod,       only : create_halo_procedure

    type(domain_t),   intent(in)  :: domain
    type(grad_c2_ecs_t)           :: grad

    integer(kind=4), parameter :: ecs_halo_width=2

    grad = grad_c2_ecs_t()
    call create_halo_procedure(grad%halo_procedure,domain,ecs_halo_width,"ECS_O")

end function create_grad_c2_ecs_operator

function create_grad_ah_sbp_proj_operator(grad_operator_name, domain, is_z) result(grad)

    use grad_sbp_proj_mod,    only : grad_sbp_proj_t
    use exchange_factory_mod, only : create_xy_points_halo_exchange
    use halo_factory_mod,     only : create_vector_halo_procedure
    use sbp_diff_21_mod,      only : sbp_diff_21_t
    use sbp_diff_42_mod,      only : sbp_diff_42_t
    use sbp_diff_43_mod,      only : sbp_diff_43_t
    use sbp_diff_63_mod,      only : sbp_diff_63_t

    character(len=*),       intent(in) :: grad_operator_name
    type(domain_t),         intent(in) :: domain
    logical,                intent(in) :: is_z
    type(grad_sbp_proj_t)              :: grad

    integer(kind=4) :: halo_width

    select case(grad_operator_name)
    case ("grad_ah_sbp_proj_21")
        grad%sbp_diff = sbp_diff_21_t()
        halo_width = 1
    case ("grad_ah_sbp_proj_42")
        grad%sbp_diff = sbp_diff_42_t()
        halo_width = 3
    case ("grad_ah_sbp_proj_43")
        grad%sbp_diff = sbp_diff_43_t()
        halo_width = 5
    case ("grad_ah_sbp_proj_63")
        grad%sbp_diff = sbp_diff_63_t()
        halo_width = 5
    case default
        call parcomm_global%abort("unknown ah_sbp_proj gradient_operator_name"// grad_operator_name)
    end select

    if(.not. is_z) then
        call create_vector_halo_procedure(grad%proj_op, domain, 0, "ecs_Ah_vec_sync_covariant")
        grad%exch_halo = create_xy_points_halo_exchange(domain%partition, domain%parcomm, &
                                         domain%topology,  halo_width, 'full')

        call domain%get_mesh(grad%mesh_p, "xy")
        call domain%get_mesh(grad%mesh_u, "xy")
        call domain%get_mesh(grad%mesh_v, "xy")
    else
        ! This branch is not implemented yet
    end if

end function create_grad_ah_sbp_proj_operator

function create_grad_c_sbp_sat_operator(grad_operator_name, domain, is_z) result(grad)

    use grad_sbp_sat_mod,     only : grad_sbp_sat_t
    use sbp_diff_c2i_21_mod,  only : sbp_diff_c2i_21_t
    use sbp_diff_c2i_42_mod,  only : sbp_diff_c2i_42_t
    use sbp_diff_c2i_63_mod,  only : sbp_diff_c2i_63_t
    use exchange_factory_mod, only : create_o_points_halo_exchange, &
                                     create_z_points_halo_exchange

    character(len=*),       intent(in) :: grad_operator_name
    type(domain_t),         intent(in) :: domain
    logical,                intent(in) :: is_z
    type(grad_sbp_sat_t)               :: grad

    integer(kind=4) :: halo_width

    select case(grad_operator_name)
    case("grad_c_sbp_sat_21","grad_c_sbp_sat_21_z")
        halo_width = 2
        grad%sbp_diff = sbp_diff_c2i_21_t()
    case("grad_c_sbp_sat_42","grad_c_sbp_sat_42_z")
        halo_width = 3
        grad%sbp_diff = sbp_diff_c2i_42_t()
    case("grad_c_sbp_sat_63", "grad_c_sbp_sat_63_z")
        halo_width = 4
        grad%sbp_diff = sbp_diff_c2i_63_t()
    case default
        call parcomm_global%abort("unknown c_sbp gradient_operator_name"// grad_operator_name)
    end select

    if(.not. is_z) then
        grad%exch_halo = create_o_points_halo_exchange( domain%partition,    &
                          domain%parcomm, domain%topology,  halo_width, 'full')

        call domain%get_mesh(grad%mesh_p, "o")
        call domain%get_mesh(grad%mesh_u, "x")
        call domain%get_mesh(grad%mesh_v, "y")
    else
        grad%exch_halo = create_z_points_halo_exchange( domain%partition,    &
                          domain%parcomm, domain%topology,  halo_width, 'full')
        call domain%get_mesh(grad%mesh_p, "z" )
        call domain%get_mesh(grad%mesh_u, "xz")
        call domain%get_mesh(grad%mesh_v, "yz")
    end if

end function create_grad_c_sbp_sat_operator


function create_grad_ch_sbp_sat_operator(grad_operator_name, domain, is_z) result(grad)

    use grad_sbp_sat_mod,     only : grad_sbp_sat_t
    use exchange_factory_mod, only : create_xy_points_halo_exchange
    use sbp_diff_i2c_21_mod,  only : sbp_diff_i2c_21_t
    use sbp_diff_i2c_42_mod,  only : sbp_diff_i2c_42_t
    use sbp_diff_i2c_63_mod,  only : sbp_diff_i2c_63_t

    character(len=*),       intent(in) :: grad_operator_name
    type(domain_t),         intent(in) :: domain
    logical,                intent(in) :: is_z
    type(grad_sbp_sat_t)               :: grad

    integer(kind=4) :: halo_width

    select case(grad_operator_name)
    case("grad_ch_sbp_sat_21")
        halo_width = 1
        grad%sbp_diff = sbp_diff_i2c_21_t()
    case("grad_ch_sbp_sat_42")
        halo_width = 3
        grad%sbp_diff = sbp_diff_i2c_42_t()
    case("grad_ch_sbp_sat_63")
        halo_width = 6
        grad%sbp_diff = sbp_diff_i2c_63_t()
    case default
        call parcomm_global%abort("unknown ch_sbp gradient_operator_name"// grad_operator_name)
    end select

    if(.not. is_z) then
        grad%exch_halo = create_xy_points_halo_exchange(domain%partition, domain%parcomm, &
                                         domain%topology,  halo_width, 'full')

        call domain%get_mesh(grad%mesh_p, "xy")
        call domain%get_mesh(grad%mesh_u, "y")
        call domain%get_mesh(grad%mesh_v, "x")
    else
        ! This branch is not implemented yet
        ! grad%exch_halo = create_z_points_halo_exchange( domain%partition,    &
        !                   domain%parcomm, domain%topology,  halo_width, 'full')
        ! call domain%get_mesh(grad%mesh_p, "z" )
        ! call domain%get_mesh(grad%mesh_u, "xz")
        ! call domain%get_mesh(grad%mesh_v, "yz")
    end if

end function create_grad_ch_sbp_sat_operator

function create_grad_a2_operator(domain, grad_operator_name) result(grad)

    use grad_a2_mod,        only : grad_a2_t
    use halo_factory_mod,   only : create_halo_procedure

    type(domain_t),   intent(in)  :: domain
    character(len=*), intent(in)  :: grad_operator_name
    type(grad_a2_t)               :: grad

    integer(kind=4), parameter :: ecs_halo_width=2, default_halo_width=1

    grad = grad_a2_t()
    if(grad_operator_name=="gradient_a2_ecs") then
        call create_halo_procedure(grad%halo_procedure,domain,ecs_halo_width,"ECS_O")
    else if(grad_operator_name=="gradient_a2_cons") then
        call create_halo_procedure(grad%halo_procedure,domain,default_halo_width,"A_default")
    else
        call parcomm_global%abort("grad_factory_mod, create_grad_a2_operator "//&
                                  "unknown gradient_a2 subtype: "// grad_operator_name)
    end if

end function create_grad_a2_operator

function create_grad_ch_halo_operator(domain, grad_operator_name) result(grad)

    use grad_ch_halo_mod,     only : grad_ch_halo_t
    use halo_factory_mod,     only : create_halo_procedure

    type(domain_t),   intent(in)  :: domain
    character(len=*), intent(in)  :: grad_operator_name
    type(grad_ch_halo_t)          :: grad

    integer(kind=4)               :: halo_width_interior

    select case (grad_operator_name)
    case("gradient_ch_ecs_halo2")
        grad%order             = 2
        grad%input_halo_width  = 0
        grad%output_halo_width = 0
        call create_halo_procedure(grad%halo,domain,grad%input_halo_width,"ECS_xy")
    case("gradient_ch_ecs_halo4")
        grad%order             = 4
        grad%input_halo_width  = 1
        grad%output_halo_width = 0
        call create_halo_procedure(grad%halo,domain,grad%input_halo_width,"ECS_xy")
    case default
        call parcomm_global%abort("create_grad_ch_halo, unknown operator: "//&
                                  grad_operator_name)
    end select

end function create_grad_ch_halo_operator

end module grad_factory_mod
