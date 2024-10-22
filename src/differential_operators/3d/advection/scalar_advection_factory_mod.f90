module scalar_advection_factory_mod

use abstract_scalar_advection3d_mod, only : scalar_advection3d_t
use domain_mod,                      only : domain_t
use parcomm_mod,                     only : parcomm_global
use generic_config_mod,              only : generic_config_t
use v_nabla_factory_mod,             only : create_v_nabla_hor_operator
use adv_z_factory_mod,               only : create_adv_z_operator
use halo_factory_mod,                only : create_halo_procedure
use grid_field_factory_mod,          only : create_grid_field
use interpolator_uv2w_factory_mod,   only : create_uv2w_interpolator
use interpolator2d_factory_mod,      only : create_vec2vec_interpolator2d
use vertical_operator_factory_mod,   only : create_vertical_operator

implicit none

contains

subroutine create_scalar_advection3d_operator(adv_op, config, domain)
    class(scalar_advection3d_t), allocatable, intent(out)   :: adv_op
    class(generic_config_t),                  intent(inout) :: config
    type(domain_t),                           intent(in)    :: domain

    character(len=:), allocatable :: operator_name
    call config%get(operator_name,"operator_name")

    select case(operator_name)
    case("advection_p_staggered")
        call create_p_3d_advection_C(adv_op, config, domain)
    case("advection_p_Ah")
        call create_p_3d_advection_Ah(adv_op, config, domain)
    case("advection_w_staggered")
        call create_w_3d_advection_C(adv_op, config, domain)
    case("advection_w_Ah")
        call create_w_3d_advection_Ah(adv_op, config, domain)
    case default
        call parcomm_global%abort("create_scalar_advection3d_operator, unknown "//&
                                  "scalar_advection_op_name: "// operator_name)
    end select
end subroutine create_scalar_advection3d_operator

subroutine create_p_3d_advection_C(adv_op, config, domain)

    use advection_p_3d_mod,      only : advection_p_C3d_t

    class(scalar_advection3d_t), allocatable, intent(out)    :: adv_op
    class(generic_config_t),                  intent(inout)  :: config
    type(domain_t),                           intent(in)     :: domain

    type(advection_p_C3d_t), allocatable :: adv_p3d
    integer(kind=4) :: halo_width
    character(len=:), allocatable :: name

    allocate(adv_p3d)

    call config%get(name,"uv2p_operator_name")
    call create_vec2vec_interpolator2d(adv_p3d%interp_uv2p_op,name,domain)
    call config%get(name,"w2p_operator_name")
    call create_vertical_operator(adv_p3d%interp_w2p_op, name)

    call config%get(name,"hor_advection_oper_name")
    call create_v_nabla_hor_operator(adv_p3d%v_nabla_op,halo_width, name)

    adv_p3d%halo_width = halo_width
    call config%get(name,"halo_procedure")
    call create_halo_procedure(adv_p3d%halo_f,domain,max(halo_width,2),name)

    call config%get(name,"z_advection_oper_name")
    call create_adv_z_operator(adv_p3d%adv_z, name)
    
    call create_grid_field(adv_p3d%up, 0, 0, domain%mesh_p)
    call create_grid_field(adv_p3d%vp, 0, 0, domain%mesh_p)
    call create_grid_field(adv_p3d%eta_dot_p, 0, 0, domain%mesh_p)
    call create_grid_field(adv_p3d%f_tend_z, 0, 0, domain%mesh_p)

    call move_alloc(adv_p3d, adv_op)

end subroutine create_p_3d_advection_C

subroutine create_p_3d_advection_Ah(adv_op, config, domain)

    use advection_p_3d_mod,      only : advection_p_Ah3d_t
    use exchange_factory_mod,    only : create_xy_points_halo_exchange

    class(scalar_advection3d_t), allocatable, intent(out)   :: adv_op
    class(generic_config_t),                  intent(inout) :: config
    type(domain_t),                           intent(in)    :: domain

    type(advection_p_Ah3d_t), allocatable :: adv_p3d
    integer(kind=4) :: halo_width
    character(len=:), allocatable :: name

    allocate(adv_p3d)

    call config%get(name,"w2p_operator_name")
    call create_vertical_operator(adv_p3d%interp_w2p_op,name)

    call config%get(name,"hor_advection_oper_name")
    call create_v_nabla_hor_operator(adv_p3d%v_nabla_op,halo_width,name)

    adv_p3d%halo_width = halo_width
    adv_p3d%exchange_p_interior =  &
              create_xy_points_halo_exchange(domain%partition, domain%parcomm, &
                                             domain%topology,  halo_width, 'full')
    call create_halo_procedure(adv_p3d%p_edge_sync,domain,1,"Ah_scalar_sync")

    call config%get(name,"z_advection_oper_name")
    call create_adv_z_operator(adv_p3d%adv_z, name)

    call create_grid_field(adv_p3d%eta_dot_p, 0, 0, domain%mesh_p)
    call create_grid_field(adv_p3d%f_tend_z, 0, 0, domain%mesh_p)

    call move_alloc(adv_p3d, adv_op)

end subroutine create_p_3d_advection_Ah

subroutine create_w_3d_advection_C(adv_op, config, domain)

    use advection_w_3d_mod,      only : advection_w_C3d_t

    class(scalar_advection3d_t), allocatable, intent(out)    :: adv_op
    class(generic_config_t),                  intent(inout)  :: config
    type(domain_t),                           intent(in)     :: domain

    type(advection_w_C3d_t), allocatable :: adv_w3d
    integer(kind=4) :: halo_width
    character(len=:), allocatable :: name, hor_part, vert_part

    allocate(adv_w3d)

    call config%get(name, "uv2w_operator_name")
    call config%get(hor_part, "uv2w_hor_part_name")
    call config%get(vert_part,"uv2w_vert_part_name")
    call create_uv2w_interpolator(adv_w3d%interp_uv2w_op, name, &
                                  hor_part, vert_part, domain)

    call config%get(name,"hor_advection_oper_name")
    call create_v_nabla_hor_operator(adv_w3d%v_nabla_op,halo_width,name)

    adv_w3d%halo_width = halo_width
    call config%get(name,"halo_procedure")
    call create_halo_procedure(adv_w3d%halo_f,domain,max(halo_width,2),name)

    call config%get(name,"z_advection_oper_name")
    call create_adv_z_operator(adv_w3d%adv_z, name)
    call create_grid_field(adv_w3d%uw, 0, 0, domain%mesh_w)
    call create_grid_field(adv_w3d%vw, 0, 0, domain%mesh_w)
    call create_grid_field(adv_w3d%f_tend_z, 0, 0, domain%mesh_w)

    call move_alloc(adv_w3d, adv_op)

end subroutine create_w_3d_advection_C

subroutine create_w_3d_advection_Ah(adv_op, config, domain)

    use advection_w_3d_mod,      only : advection_w_Ah3d_t
    use exchange_factory_mod,    only : create_xyz_points_halo_exchange

    class(scalar_advection3d_t), allocatable, intent(out)   :: adv_op
    class(generic_config_t),                  intent(inout) :: config
    type(domain_t),                           intent(in)    :: domain

    type(advection_w_Ah3d_t), allocatable :: adv_w3d
    integer(kind=4) :: halo_width
    character(len=:), allocatable :: name, hor_part, vert_part

    allocate(adv_w3d)

    call config%get(name,"uv2w_operator_name")
    call config%get(vert_part,"uv2w_vert_part_name")
    hor_part = "None"
    call create_uv2w_interpolator(adv_w3d%interp_uv2w_op, name, &
                                  hor_part, vert_part, domain)

    call config%get(name,"hor_advection_oper_name")
    call create_v_nabla_hor_operator(adv_w3d%v_nabla_op,halo_width,name)

    adv_w3d%halo_width = halo_width
    adv_w3d%exchange_interior =  &
              create_xyz_points_halo_exchange(domain%partition, domain%parcomm, &
                                             domain%topology,  halo_width, 'full')
    call create_halo_procedure(adv_w3d%edge_sync,domain,1,"Ah_scalar_sync_z")

    call config%get(name,"z_advection_oper_name")
    call create_adv_z_operator(adv_w3d%adv_z, name)

    call create_grid_field(adv_w3d%uw, 0, 0, domain%mesh_w)
    call create_grid_field(adv_w3d%vw, 0, 0, domain%mesh_w)
    call create_grid_field(adv_w3d%f_tend_z, 0, 0, domain%mesh_w)

    call move_alloc(adv_w3d, adv_op)

end subroutine create_w_3d_advection_Ah

end module scalar_advection_factory_mod
