module config_mixvec_transform_mod

use config_mod, only : config_t

implicit none

type, extends(config_t) :: config_mixvec_transform_t
    character(:), allocatable :: w2p_interp_name, p2w_interp_name
    character(:), allocatable :: uv2p_interp_name, p2uv_interp_name
    contains
        procedure :: parse
end type config_mixvec_transform_t

contains

subroutine parse(this, config_string)

    class(config_mixvec_transform_t),  intent(inout) :: this
    character(len=*),                  intent(in)    :: config_string

    character(len=256) :: w2p_interp_name, p2w_interp_name, &
                          uv2p_interp_name, p2uv_interp_name

    namelist /mixvec_transform_config/ w2p_interp_name, p2w_interp_name, &
                                       uv2p_interp_name, p2uv_interp_name

    read(config_string, mixvec_transform_config)

    this%w2p_interp_name  = w2p_interp_name
    this%p2w_interp_name  = p2w_interp_name
    this%uv2p_interp_name = uv2p_interp_name
    this%p2uv_interp_name = p2uv_interp_name

end subroutine parse

end module config_mixvec_transform_mod
