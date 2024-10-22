module laplace_factory_mod

use domain_mod,                   only : domain_t
use abstract_laplace_mod,         only : laplace_operator_t
use parcomm_mod,                  only : parcomm_global
use laplace_true_hor_factory_mod, only : create_laplace_true_hor

implicit none

contains

subroutine create_laplace_operator(laplace_operator,laplace_operator_name,domain)
    character(len=*), intent(in) :: laplace_operator_name
    type(domain_t),   intent(in) :: domain
    !output:
    class(laplace_operator_t), allocatable, intent(out) :: laplace_operator

    if(laplace_operator_name(1:15) == "divgrad_laplace") then
        call create_divgrad_laplace(laplace_operator,laplace_operator_name, domain)
    else if(laplace_operator_name == "laplace_ch_halo2" .or. &
            laplace_operator_name == "laplace_ch_halo4") then
        call create_laplace_ch_halo(laplace_operator,laplace_operator_name, domain)
    else if(laplace_operator_name(1:23) == "laplace_ah_sbp21_narrow" .or. &
            laplace_operator_name(1:23) == "laplace_ah_sbp42_narrow" .or. &
            laplace_operator_name(1:23) == "laplace_ah_sbp63_narrow") then
        call create_laplace_ah_sbp_narrow(laplace_operator, laplace_operator_name, domain)
    else if(laplace_operator_name(1:20) == "laplace_3d_true_hor_") then
        call create_laplace_true_hor(laplace_operator,laplace_operator_name,domain)
    else
        call parcomm_global%abort("create laplace_operator, unknown laplace_operator_name:"//&
                                   laplace_operator_name)
    end if
end subroutine create_laplace_operator

subroutine create_divgrad_laplace(laplace_operator,laplace_operator_name, domain)

    use divgrad_laplace_mod,    only : divgrad_laplace_t
    use grad_factory_mod,       only : create_grad_operator
    use div_factory_mod,        only : create_div_operator
    use co2contra_factory_mod,  only : create_co2contra_operator
    use grid_field_factory_mod, only : create_grid_field
    use mesh_mod,               only : mesh_t

    character(len=*), intent(in) :: laplace_operator_name
    type(domain_t),   intent(in) :: domain
    !output:
    class(laplace_operator_t), allocatable, intent(out) :: laplace_operator
    !local
    type(divgrad_laplace_t), allocatable :: laplace
    character(len=:), allocatable :: grad_name, div_name, co2contra_name
    integer(kind=4)  :: halo_width
    type(mesh_t), pointer :: mesh_x, mesh_y

    allocate(laplace)

    call domain%get_mesh(mesh_x,"u")
    call domain%get_mesh(mesh_y,"v")

    select case(laplace_operator_name)
    case("divgrad_laplace_c_sbp21")
        div_name = "div_c_sbp_sat_21"
        co2contra_name = "co2contra_c_sbp21"
        grad_name = "grad_c_sbp_sat_21"
    case("divgrad_laplace_c_sbp21_z")
        div_name = "div_c_sbp_sat_21_z"
        co2contra_name = "co2contra_c_sbp21_new_z"
        grad_name = "grad_c_sbp_sat_21_z"
        call domain%get_mesh(mesh_x,"xz")
        call domain%get_mesh(mesh_y,"yz")
    case("divgrad_laplace_c_sbp42")
        div_name = "div_c_sbp_sat_42"
        co2contra_name = "co2contra_c_sbp42"
        grad_name = "grad_c_sbp_sat_42"
    case("divgrad_laplace_c_sbp63")
        div_name = "div_c_sbp_sat_63"
        co2contra_name = "co2contra_c_sbp63_new"
        grad_name = "grad_c_sbp_sat_63"
    case("divgrad_laplace_c_sbp42_z")
        div_name = "div_c_sbp_sat_42_z"
        co2contra_name = "co2contra_c_sbp42_new_z"
        grad_name = "grad_c_sbp_sat_42_z"
        call domain%get_mesh(mesh_x,"xz")
        call domain%get_mesh(mesh_y,"yz")
    case("divgrad_laplace_c_sbp63_z")
        div_name = "div_c_sbp_sat_63_z"
        co2contra_name = "co2contra_c_sbp63_new_z"
        grad_name = "grad_c_sbp_sat_63_z"
    case("divgrad_laplace_ch_sbp21")
        div_name = "div_ch_sbp_sat_21"
        co2contra_name = "co2contra_ch_sbp21"
        grad_name = "grad_ch_sbp_sat_21"
    case("divgrad_laplace_ch_sbp42")
        div_name = "div_ch_sbp_sat_42"
        co2contra_name = "co2contra_ch_sbp42"
        grad_name = "grad_ch_sbp_sat_42"
    case("divgrad_laplace_ah_sbp21")
        div_name = "div_ah_sbp_proj_21"
        co2contra_name = "co2contra_colocated"
        grad_name = "grad_ah_sbp_proj_21"
    case("divgrad_laplace_ah_sbp42")
        div_name = "div_ah_sbp_proj_42"
        co2contra_name = "co2contra_colocated"
        grad_name = "grad_ah_sbp_proj_42"
    case default
        call parcomm_global%abort("create_divgrad_laplace, incorrect laplace_operator_name:"//&
                                  laplace_operator_name)
    end select

    laplace%div_operator       = create_div_operator(domain, div_name)
    laplace%co2contra_operator = create_co2contra_operator(domain, co2contra_name)
    laplace%grad_operator      = create_grad_operator(domain,grad_name)

    !WORKAROUND
    halo_width = 6
    call create_grid_field(laplace%gx,  halo_width, 0, mesh_x)
    call create_grid_field(laplace%gy,  halo_width, 0, mesh_y)
    call create_grid_field(laplace%gxt, halo_width, 0, mesh_x)
    call create_grid_field(laplace%gyt, halo_width, 0, mesh_y)

    call move_alloc(laplace, laplace_operator)

end subroutine create_divgrad_laplace

subroutine create_laplace_ch_halo(laplace_operator,laplace_operator_name, domain)

    use laplace_ch_halo_mod, only : laplace_ch_halo_t
    use halo_factory_mod,    only : create_halo_procedure

    character(len=*), intent(in) :: laplace_operator_name
    type(domain_t),   intent(in) :: domain
    !output:
    class(laplace_operator_t), allocatable, intent(out) :: laplace_operator
    !local
    type(laplace_ch_halo_t), allocatable :: laplace

    allocate(laplace)
    select case(laplace_operator_name)
    case("laplace_ch_halo2")
        laplace%halo_width = 1
        laplace%order      = 2
    case("laplace_ch_halo4")
        laplace%halo_width = 3
        laplace%order      = 4
    case default
        call parcomm_global%abort("create_laplace_ch_halo, incorrect laplace_operator_name:"//&
                                  laplace_operator_name)
    end select


    call create_halo_procedure(laplace%halo_f,   domain,laplace%halo_width,"ECS_xy")
    call create_halo_procedure(laplace%edge_sync,domain,1,"Ah_scalar_sync")

    call move_alloc(laplace, laplace_operator)

    contains
    function find_position_in_str(char,str) result(pos)
        character(len=1), intent(in) :: char
        character(len=*), intent(in) :: str
        integer(kind=4)  :: pos
        integer(kind=4)  :: i
        pos = -1
        do i=1, len(str)
            if(str(i:i) == char) then
                pos = i
                return
            end if
        end do
        end
end subroutine create_laplace_ch_halo
subroutine create_laplace_ah_sbp_narrow(laplace_operator, laplace_operator_name, domain)

    use laplace_Ah_sbp21_narrow_mod, only : laplace_Ah_sbp21_narrow_t
    use laplace_Ah_sbp42_narrow_mod, only : laplace_Ah_sbp42_narrow_t
    use laplace_Ah_sbp63_narrow_mod, only : laplace_Ah_sbp63_narrow_t
    use exchange_factory_mod,        only : create_xy_points_halo_exchange, &
                                            create_xyz_points_halo_exchange
    use grid_field_factory_mod,      only : create_grid_field
    use sbp_diff_21_mod,             only : sbp_diff_21_t
    use sbp_diff_42_mod,             only : sbp_diff_42_t
    use sbp_diff_63_mod,             only : sbp_diff_63_t
    use halo_factory_mod,            only : create_halo_procedure

    character(len=*), intent(in) :: laplace_operator_name
    type(domain_t),   intent(in) :: domain
    !output:
    class(laplace_operator_t), allocatable, intent(out) :: laplace_operator
    !local
    class(laplace_Ah_sbp21_narrow_t), allocatable :: laplace
    integer(kind=4)  :: halo_width
    logical          :: is_z

    select case(laplace_operator_name)
    case("laplace_ah_sbp21_narrow", "laplace_ah_sbp21_narrow_z")
        allocate(laplace_Ah_sbp21_narrow_t :: laplace)
        laplace%sbp_d1_op = sbp_diff_21_t()
        halo_width = 1
    case("laplace_ah_sbp42_narrow", "laplace_ah_sbp42_narrow_z")
        allocate(laplace_Ah_sbp42_narrow_t :: laplace)
        laplace%sbp_d1_op = sbp_diff_42_t()
        halo_width = 6
    case("laplace_ah_sbp63_narrow", "laplace_ah_sbp63_narrow_z")
        allocate(laplace_Ah_sbp63_narrow_t :: laplace)
        laplace%sbp_d1_op = sbp_diff_63_t()
        halo_width = 9
    case default
        call parcomm_global%abort("create_laplace_ah_sbp_narrow, unknown operator name "// laplace_operator_name)
    end select

    is_z = .false.
    if(len(laplace_operator_name) >= 25) is_z = (laplace_operator_name(24:25)=="_z")

    if(.not. is_z) then
        laplace%exchange = create_xy_points_halo_exchange(domain%partition, domain%parcomm, &
                                                          domain%topology,  halo_width, 'full')
        call domain%get_mesh(laplace%mesh,"xy")
        call create_halo_procedure(laplace%edge_sync,domain,1,"Ah_scalar_sync")
    else
        laplace%exchange = create_xyz_points_halo_exchange(domain%partition, domain%parcomm, &
                                                          domain%topology,  halo_width, 'full')
        call domain%get_mesh(laplace%mesh,"xyz")
        call create_halo_procedure(laplace%edge_sync,domain,1,"Ah_scalar_sync_z")
    end if

    call create_grid_field(laplace%d2f_x, 0, 0, laplace%mesh)
    call create_grid_field(laplace%d2f_y, 0, 0, laplace%mesh)

    call create_grid_field(laplace%d1f_x, halo_width, 0, laplace%mesh)
    call create_grid_field(laplace%d1f_y, halo_width, 0, laplace%mesh)

    call create_grid_field(laplace%d2f_xy, 0, 0, laplace%mesh)
    call create_grid_field(laplace%d2f_yx, 0, 0, laplace%mesh)

    call move_alloc(laplace, laplace_operator)

end subroutine create_laplace_ah_sbp_narrow

end module laplace_factory_mod
