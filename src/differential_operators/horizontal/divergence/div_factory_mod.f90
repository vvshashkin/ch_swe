module div_factory_mod

use domain_mod,        only : domain_t
use mesh_mod,          only : mesh_t
use abstract_div_mod,  only : div_operator_t
use parcomm_mod,       only : parcomm_global

implicit none

contains

function create_div_operator(domain, div_operator_name) result(div)
    type(domain_t),    intent(in)    :: domain
    character(len=*),  intent(in)    :: div_operator_name

    class(div_operator_t), allocatable :: div

    select case(div_operator_name)
    case("div_ah_sbp_proj_21", "div_ah_sbp_proj_42", "div_ah_sbp_proj_43", "div_ah_sbp_proj_63")
        div = create_div_ah_sbp_proj_operator(div_operator_name, domain)

    case("div_c_sbp_sat_21", "div_c_sbp_sat_42", "div_c_sbp_sat_63")
        div = create_div_c_sbp_sat_operator(div_operator_name, domain, is_z = .false.)

    case("div_c_sbp_sat_21_z", "div_c_sbp_sat_42_z", "div_c_sbp_sat_63_z")
        div = create_div_c_sbp_sat_operator(div_operator_name, domain, is_z = .true.)

    case("div_ch_sbp_sat_21", "div_ch_sbp_sat_42", "div_ch_sbp_sat_63", &
         "div_ch_sbp_sat_21_proj", "div_ch_sbp_sat_42_proj", "div_ch_sbp_sat_63_proj")
        div = create_div_ch_sbp_sat_operator(div_operator_name, domain, is_z = .false.)

    case("divergence_c2")
        div = create_div_c2_operator(domain)

    case("divergence_a2_ecs","divergence_a2_cons", "divergence_a2_fv")
        div = create_div_a2_operator(domain,div_operator_name)

    case default
        call parcomm_global%abort("unknown divergence operator: "//div_operator_name)
    end select
end function create_div_operator

function create_div_c2_operator(domain) result(div)

    use div_c2_mod,             only : div_c2_t
    use halo_factory_mod,       only : create_vector_halo_procedure

    type(domain_t),   intent(in)  :: domain
    type(div_c2_t)                :: div

    !div%exch_halo = create_symmetric_halo_vec_exchange_C(domain%partition, domain%parcomm, &
    !                                                     domain%topology, halo_width, 'full')
    call create_vector_halo_procedure(div%halo_procedure, domain, 1, "C_vec_default")
end function create_div_c2_operator

function create_div_ah_sbp_proj_operator(div_operator_name, domain) result(div)

    use div_sbp_proj_mod,       only : div_sbp_proj_t
    use exchange_factory_mod,   only : create_xy_points_halo_exchange
    use halo_factory_mod,       only : create_halo_procedure
    use sbp_diff_21_mod,        only : sbp_diff_21_t
    use sbp_diff_42_mod,        only : sbp_diff_42_t
    use sbp_diff_43_mod,        only : sbp_diff_43_t
    use sbp_diff_63_mod,        only : sbp_diff_63_t
    use grid_field_factory_mod, only : create_grid_field

    character(len=*), intent(in)  :: div_operator_name
    type(domain_t),   intent(in)  :: domain
    type(div_sbp_proj_t)          :: div

    integer(kind=4)            :: halo_width

    select case(div_operator_name)
    case ("div_ah_sbp_proj_21")
        halo_width = 1
        div%sbp_diff = sbp_diff_21_t()
    case ("div_ah_sbp_proj_42")
        halo_width = 3
        div%sbp_diff = sbp_diff_42_t()
    case ("div_ah_sbp_proj_43")
        halo_width = 5
        div%sbp_diff = sbp_diff_43_t()
    case ("div_ah_sbp_proj_63")
        halo_width = 5
        div%sbp_diff = sbp_diff_63_t()
    case default
        call parcomm_global%abort("div_factory_mod, create_div_ah_sbp_proj_operator"// &
                                  " - unknown div operator: "//div_operator_name)
    end select

    call domain%get_mesh(div%mesh_u, "xy")
    call domain%get_mesh(div%mesh_v, "xy")
    call domain%get_mesh(div%mesh_p, "xy")

    call create_grid_field(div%Ju, halo_width, 0, div%mesh_u)
    call create_grid_field(div%Jv, halo_width, 0, div%mesh_v)

    call create_grid_field(div%Dx, 0, 0, div%mesh_p)
    call create_grid_field(div%Dy, 0, 0, div%mesh_p)

    div%exch_halo =  create_xy_points_halo_exchange(domain%partition, domain%parcomm, &
                                             domain%topology,  halo_width, 'full')
    call create_halo_procedure(div%proj_op, domain, 1, "Ah_scalar_sync")

end function create_div_ah_sbp_proj_operator

function create_div_c_sbp_sat_operator(div_operator_name, domain, is_z) result(div)

    use div_sbp_SAT_mod,        only : div_sbp_SAT_t
    use sbp_diff_i2c_21_mod,    only : sbp_diff_i2c_21_t
    use sbp_diff_i2c_42_mod,    only : sbp_diff_i2c_42_t
    use sbp_diff_i2c_63_mod,    only : sbp_diff_i2c_63_t
    use grid_field_factory_mod, only : create_grid_field
    use exchange_factory_mod,   only : create_symmetric_halo_vec_exchange_C, &
                                       create_symmetric_halo_vec_exchange_Cz

    character(len=*),   intent(in) :: div_operator_name
    type(domain_t),     intent(in) :: domain
    logical,            intent(in) :: is_z
    type(div_sbp_SAT_t)            :: div

    integer(kind=4) :: halo_width

    select case(div_operator_name)
    case("div_c_sbp_sat_21", "div_c_sbp_sat_21_z")
        halo_width = 1
        div%sbp_diff = sbp_diff_i2c_21_t()
    case("div_c_sbp_sat_42", "div_c_sbp_sat_42_z")
        halo_width = 3
        div%sbp_diff  = sbp_diff_i2c_42_t()
    case("div_c_sbp_sat_63", "div_c_sbp_sat_63_z")
        halo_width = 5
        div%sbp_diff = sbp_diff_i2c_63_t()
    case default
        call parcomm_global%abort("div_factory_mod, create_div_sbp_sat_operator"// &
                                  " - unknown div operator: "//div_operator_name)
    end select

    if (is_z) then
        call domain%get_mesh(div%mesh_u, "xz")
        call domain%get_mesh(div%mesh_v, "yz")
        call domain%get_mesh(div%mesh_p, "z")
        div%exch_halo = create_symmetric_halo_vec_exchange_Cz(domain%partition, domain%parcomm, &
                                                             domain%topology, halo_width, 'full')
    else
        call domain%get_mesh(div%mesh_u, "x")
        call domain%get_mesh(div%mesh_v, "y")
        call domain%get_mesh(div%mesh_p, "o")
        div%exch_halo = create_symmetric_halo_vec_exchange_C(domain%partition, domain%parcomm, &
                                                             domain%topology, halo_width, 'full')
    end if

    call create_grid_field(div%Ju, halo_width + 1, 0, div%mesh_u)
    call create_grid_field(div%Jv, halo_width + 1, 0, div%mesh_v)
    call create_grid_field(div%Dx, 0             , 0, div%mesh_p)
    call create_grid_field(div%Dy, 0             , 0, div%mesh_p)

end function create_div_c_sbp_sat_operator

function create_div_ch_sbp_sat_operator(div_operator_name, domain, is_z) result(div)

    use div_sbp_SAT_mod,        only : div_sbp_SAT_t
    use sbp_diff_c2i_21_mod,    only : sbp_diff_c2i_21_t
    use sbp_diff_c2i_42_mod,    only : sbp_diff_c2i_42_t
    use sbp_diff_c2i_63_mod,    only : sbp_diff_c2i_63_t
    use grid_field_factory_mod, only : create_grid_field
    use exchange_factory_mod,   only : create_symmetric_halo_vec_exchange_Ch
    use halo_factory_mod,       only : create_halo_procedure

    character(len=*),   intent(in) :: div_operator_name
    type(domain_t),     intent(in) :: domain
    logical,            intent(in) :: is_z
    type(div_sbp_SAT_t)            :: div

    integer(kind=4) :: halo_width

    select case(div_operator_name)

    case("div_ch_sbp_sat_21", "div_ch_sbp_sat_21_proj")
        halo_width = 2
        div%sbp_diff = sbp_diff_c2i_21_t()
    case("div_ch_sbp_sat_42", "div_ch_sbp_sat_42_proj")
        halo_width = 3
        div%sbp_diff  = sbp_diff_c2i_42_t()
    case("div_ch_sbp_sat_63", "div_ch_sbp_sat_63_proj")
        halo_width = 6
        div%sbp_diff  = sbp_diff_c2i_63_t()
    case default
        call parcomm_global%abort("div_factory_mod, create_div_sbp_sat_operator"// &
                                  " - unknown div operator: "//div_operator_name)
    end select

    if (is_z) then
        ! This branch are not implemented yet!!!
        ! call domain%get_mesh(div%mesh_u, "yz")
        ! call domain%get_mesh(div%mesh_v, "xz")
        ! call domain%get_mesh(div%mesh_p, "xyz")
        !div%exch_halo = create_symmetric_halo_vec_exchange_Chz(domain%partition, domain%parcomm, &
        !                                                     domain%topology, halo_width, 'full')
    else
        call domain%get_mesh(div%mesh_u, "y")
        call domain%get_mesh(div%mesh_v, "x")
        call domain%get_mesh(div%mesh_p, "xy")
        div%exch_halo = create_symmetric_halo_vec_exchange_Ch(domain%partition, domain%parcomm, &
                                                              domain%topology, halo_width, 'full')
    end if

    call create_grid_field(div%Ju, halo_width + 1, 0, div%mesh_u)
    call create_grid_field(div%Jv, halo_width + 1, 0, div%mesh_v)
    call create_grid_field(div%Dx, 0             , 0, div%mesh_p)
    call create_grid_field(div%Dy, 0             , 0, div%mesh_p)

    select case (div_operator_name)
    case("div_ch_sbp_sat_63_proj", "div_ch_sbp_sat_42_proj", "div_ch_sbp_sat_21_proj")
        call create_halo_procedure(div%proj_op, domain, 1, "Ah_scalar_sync")
    end select

end function create_div_ch_sbp_sat_operator

function create_div_a2_operator(domain, div_operator_name) result(div)

    use div_a2_mod,       only : div_a2_t
    use halo_factory_mod, only : create_vector_halo_procedure

    type(domain_t),   intent(in)  :: domain
    character(len=*), intent(in)  :: div_operator_name
    type(div_a2_t) :: div

    integer(kind=4) :: ecs_halo_width, default_halo_width

    ecs_halo_width = 2
    default_halo_width = 1

    if (div_operator_name == "divergence_a2_ecs") then
        call create_vector_halo_procedure(div%halo_procedure, domain, ecs_halo_width, "ecs_A_vec")
        div%subtype = "default"
    else if (div_operator_name == "divergence_a2_cons") then
        call create_vector_halo_procedure(div%halo_procedure, domain, default_halo_width, "A_vec_default")
        div%subtype = "cons"
    else if (div_operator_name == "divergence_a2_fv") then
        call create_vector_halo_procedure(div%halo_procedure, domain, default_halo_width, "A_vec_default")
        div%subtype = "fv"
    else
        call parcomm_global%abort("div_factory_mod, create_div_a2_operator, "// &
                                  "unknown divergence_a2 operator subtype: "// div_operator_name)
    end if

end function create_div_a2_operator

end module div_factory_mod
