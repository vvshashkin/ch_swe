module curl_factory_mod

use domain_mod,         only : domain_t
use abstract_curl_mod,  only : curl_operator_t
use parcomm_mod,        only : parcomm_global

implicit none

contains

subroutine create_curl_operator(curl, curl_operator_name, domain)

    class(curl_operator_t), allocatable, intent(out) :: curl
    character(len=*),                    intent(in)  :: curl_operator_name
    type(domain_t),                      intent(in)  :: domain

    if (curl_operator_name(1:8) == "curl_div") then
        call create_curl_operator_div_based(curl, curl_operator_name(6:), domain)
    else if (curl_operator_name == "curl_c_sbp_sat_21" .or. &
             curl_operator_name == "curl_c_sbp_sat_42" .or. &
             curl_operator_name == "curl_c_sbp_sat_63" .or. &
             curl_operator_name == "curl_ch_sbp_sat_21" .or. &
             curl_operator_name == "curl_ch_sbp_sat_42" .or. &
             curl_operator_name == "curl_ch_sbp_sat_63") then
        curl = create_curl_c_sbp_sat_operator(curl_operator_name, domain)
    else
        call parcomm_global%abort("unknown curl operator name: "// curl_operator_name)
    end if
end subroutine create_curl_operator

subroutine create_curl_operator_div_based(curl, div_operator_name, domain)

    use curl_div_based_mod,     only : curl_div_based_t
    use div_factory_mod,        only : create_div_operator
    use grid_field_factory_mod, only : create_grid_field

    class(curl_operator_t), allocatable, intent(out) :: curl
    character(len=*),                    intent(in)  :: div_operator_name
    type(domain_t),                      intent(in)  :: domain

    type(curl_div_based_t), allocatable :: curl_div_based

    integer(kind=4) :: halo_width

    !This is temporary solution
    !We need a way to get information about required halo width from operator
    halo_width = 5

    allocate(curl_div_based)

    curl_div_based%div_op = create_div_operator(domain, div_operator_name)

    call create_grid_field(curl_div_based%ut, halo_width, 0, domain%mesh_u)
    call create_grid_field(curl_div_based%vt, halo_width, 0, domain%mesh_v)

    call move_alloc(curl_div_based, curl)

end subroutine create_curl_operator_div_based

function create_curl_c_sbp_sat_operator(curl_operator_name, domain) result(curl)

    use curl_sbp_sat_mod,       only : curl_sbp_SAT_t
    use sbp_diff_c2i_21_mod,    only : sbp_diff_c2i_21_t
    use sbp_diff_c2i_42_mod,    only : sbp_diff_c2i_42_t
    use sbp_diff_c2i_63_mod,    only : sbp_diff_c2i_63_t
    use sbp_diff_i2c_21_mod,    only : sbp_diff_i2c_21_t
    use sbp_diff_i2c_42_mod,    only : sbp_diff_i2c_42_t
    use sbp_diff_i2c_63_mod,    only : sbp_diff_i2c_63_t
    use grid_field_factory_mod, only : create_grid_field
    use exchange_factory_mod,   only : create_symmetric_halo_vec_exchange_C, create_symmetric_halo_vec_exchange_Ch
    use halo_factory_mod,       only : create_halo_procedure

    character(len=*),   intent(in) :: curl_operator_name
    type(domain_t),     intent(in) :: domain
    type(curl_sbp_sat_t)           :: curl

    integer(kind=4) :: halo_width

    select case(curl_operator_name)
    case("curl_c_sbp_sat_21")
        halo_width = 2
        curl%sbp_diff = sbp_diff_c2i_21_t()
    case("curl_c_sbp_sat_42")
        halo_width = 4
        curl%sbp_diff  = sbp_diff_c2i_42_t()
    case("curl_c_sbp_sat_63")
        halo_width = 6
        curl%sbp_diff = sbp_diff_c2i_63_t()
    case("curl_ch_sbp_sat_21")
        halo_width = 2
        curl%sbp_diff = sbp_diff_i2c_21_t()
    case("curl_ch_sbp_sat_42")
        halo_width = 4
        curl%sbp_diff  = sbp_diff_i2c_42_t()
    case("curl_ch_sbp_sat_63")
        halo_width = 6
        curl%sbp_diff = sbp_diff_i2c_63_t()
    case default
        call parcomm_global%abort("curl_factory_mod, create_curl_sbp_sat_operator"// &
                                  " - unknown curl operator: "//curl_operator_name)
    end select

    select case(curl_operator_name)
    case("curl_c_sbp_sat_21", "curl_c_sbp_sat_42", "curl_c_sbp_sat_63")

        call domain%get_mesh(curl%mesh_w, "xy")
        curl%exch_halo = create_symmetric_halo_vec_exchange_C(domain%partition, domain%parcomm, &
                                                              domain%topology, halo_width, 'full')

        call create_halo_procedure(curl%sync_edges, domain, 1, "Ah_scalar_sync")

    case("curl_ch_sbp_sat_21", "curl_ch_sbp_sat_42", "curl_ch_sbp_sat_63")

        call domain%get_mesh(curl%mesh_w, "o")
        curl%exch_halo = create_symmetric_halo_vec_exchange_Ch(domain%partition, domain%parcomm, &
                                                               domain%topology, halo_width, 'full')
    end select

    call create_grid_field(curl%Dx, 0, 0, curl%mesh_w)
    call create_grid_field(curl%Dy, 0, 0, curl%mesh_w)

end function create_curl_c_sbp_sat_operator

end module curl_factory_mod
