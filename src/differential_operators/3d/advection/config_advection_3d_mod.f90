module config_advection_3d_mod

use config_mod,  only : config_t
use parcomm_mod, only : parcomm_global

implicit none

type, extends(config_t) :: config_p_advection_t
    character(len=:), allocatable :: hor_advection_oper_name
    character(len=:), allocatable :: z_advection_oper_name
    character(len=:), allocatable :: w2p_operator_name
    character(len=:), allocatable :: uv2p_operator_name
    character(len=:), allocatable :: p_halo
    contains
    procedure :: parse => parse_p_advection_config
end type

type, extends(config_t) :: config_w_advection_t
    character(len=:), allocatable :: hor_advection_oper_name
    character(len=:), allocatable :: z_advection_oper_name
    character(len=:), allocatable :: uv2w_operator_name
    character(len=:), allocatable :: uv2w_hor_part_name
    character(len=:), allocatable :: uv2w_vert_part_name
    character(len=:), allocatable :: w_halo
    contains
    procedure :: parse => parse_w_advection_config
end type

type, extends(config_t) :: config_vector_advection_3d_t
    character(len=:), allocatable :: uv_hor_advection_oper_name
    character(len=:), allocatable :: uv_ver_advection_oper_name
    character(len=:), allocatable :: w2uv_operator_name
    character(len=:), allocatable :: w2uv_hor_part_name
    character(len=:), allocatable :: w2uv_vert_part_name
    character(len=:), allocatable :: w_advection_oper_name
    class(config_t),  allocatable :: w_advection_config
    contains
    procedure :: parse => parse_vector_advection_config
end type

contains

subroutine get_advection_3d_config(config, advection_3d_operator_name)
    character(len=*), intent(in) :: advection_3d_operator_name
    class(config_t), allocatable, intent(out) :: config

    select case(advection_3d_operator_name)
    case("advection_p_staggered", "advection_p_Ah")
        allocate(config_p_advection_t :: config)
    case("advection_w_staggered", "advection_w_Ah")
        allocate(config_w_advection_t :: config)
    case("shallow_atm_staggered_vector_advection")
        allocate(config_vector_advection_3d_t :: config)
    case default
        call parcomm_global%abort("get_advection_3d_config error: "     // &
                                  "unknown advection_3d_operator name " // &
                                   advection_3d_operator_name)
    end select
end subroutine get_advection_3d_config

subroutine parse_p_advection_config(this, config_string)

    class(config_p_advection_t),  intent(inout) :: this
    character(len=*),             intent(in)    :: config_string

    character(len=256) :: hor_advection_oper_name, z_advection_oper_name, &
                          w2p_operator_name, uv2p_operator_name, p_halo

    namelist /p_advection_conf/ hor_advection_oper_name, z_advection_oper_name, &
                                w2p_operator_name, uv2p_operator_name, p_halo

    read(config_string, p_advection_conf)

    this%hor_advection_oper_name = trim(hor_advection_oper_name)
    this%z_advection_oper_name   = trim(z_advection_oper_name)
    this%w2p_operator_name       = trim(w2p_operator_name)
    this%uv2p_operator_name      = trim(uv2p_operator_name)
    this%p_halo                  = trim(p_halo)

end subroutine parse_p_advection_config

subroutine parse_w_advection_config(this, config_string)

    class(config_w_advection_t),  intent(inout) :: this
    character(len=*),             intent(in)    :: config_string

    character(len=256) :: hor_advection_oper_name, z_advection_oper_name, &
                          uv2w_operator_name, uv2w_hor_part_name, uv2w_vert_part_name, &
                          w_halo

    namelist /w_advection_conf/ hor_advection_oper_name, z_advection_oper_name, &
                                uv2w_operator_name, uv2w_hor_part_name, uv2w_vert_part_name, &
                                w_halo

    read(config_string, w_advection_conf)

    this%hor_advection_oper_name = trim(hor_advection_oper_name)
    this%z_advection_oper_name   = trim(z_advection_oper_name)
    this%uv2w_operator_name      = trim(uv2w_operator_name)
    this%uv2w_hor_part_name      = trim(uv2w_hor_part_name)
    this%uv2w_vert_part_name     = trim(uv2w_vert_part_name)
    this%w_halo                  = trim(w_halo)

end subroutine parse_w_advection_config

subroutine parse_vector_advection_config(this, config_string)

    class(config_vector_advection_3d_t),  intent(inout) :: this
    character(len=*),                     intent(in)    :: config_string

    character(len=256) :: uv_hor_advection_oper_name, uv_ver_advection_oper_name, &
                          w2uv_operator_name, w2uv_hor_part_name, w2uv_vert_part_name, &
                          w_advection_oper_name
    character(len=512) :: w_advection_oper_config_str

    namelist /vector_advection_3d_conf/  uv_hor_advection_oper_name, uv_ver_advection_oper_name, &
                                         w2uv_operator_name, w2uv_hor_part_name, w2uv_vert_part_name, &
                                         w_advection_oper_name, w_advection_oper_config_str

    read(config_string, vector_advection_3d_conf)

    this%uv_hor_advection_oper_name = trim(uv_hor_advection_oper_name)
    this%uv_ver_advection_oper_name = trim(uv_ver_advection_oper_name)
    this%w2uv_operator_name         = trim(w2uv_operator_name)
    this%w2uv_hor_part_name         = trim(w2uv_hor_part_name)
    this%w2uv_vert_part_name        = trim(w2uv_vert_part_name)
    this%w_advection_oper_name      = trim(w_advection_oper_name)
    call get_advection_3d_config(this%w_advection_config,this%w_advection_oper_name)
    call this%w_advection_config%parse(trim(w_advection_oper_config_str))

end subroutine parse_vector_advection_config

end module config_advection_3d_mod
