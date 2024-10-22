module mixvec_transform_factory_mod

use abstract_mixvec_transform_mod,      only : mixvec_transform_t
use mixvec_transform_colocated_mod,     only : mixvec_transform_colocated_t
use mixvec_transform_hor_colocated_mod, only : mixvec_transform_hor_colocated_t
use mixvec_transform_staggered_mod,     only : mixvec_transform_staggered_t
use domain_mod,                         only : domain_t
use config_mod,                         only : config_t
use vertical_operator_factory_mod,      only : create_vertical_operator
use interpolator2d_factory_mod,         only : create_vec2vec_interpolator2d
use config_mixvec_transform_mod,        only : config_mixvec_transform_t
use grid_field_factory_mod,             only : create_grid_field
use parcomm_mod,                        only : parcomm_global

implicit none

contains

subroutine create_mixvec_transform(mixvec_transform,mixvec_transform_name,config,domain)
    class(mixvec_transform_t), allocatable, intent(out) :: mixvec_transform
    character(len=*),                       intent(in)  :: mixvec_transform_name
    class(config_t),                        intent(in)  :: config
    type(domain_t),                         intent(in)  :: domain

    type(config_mixvec_transform_t) :: config_loc

    select type(config)
    type is (config_mixvec_transform_t)
        config_loc = config
    class default
        call parcomm_global%abort("wrong config type in create_mixvec_transform")
    end select

    select case(mixvec_transform_name)
    case("mixvec_colocated")
        mixvec_transform = mixvec_transform_colocated_t()
    case("mixvec_hor_colocated")
        call create_hor_colocated_mixvec_transform(mixvec_transform,config_loc,domain)
    case("mixvec_staggered")
        call create_staggered_mixvec_transform(mixvec_transform,config_loc,domain)
    case default
        call parcomm_global%abort("mixvec_transform_factory_mod, incorrect mixvec_transform_name: "//mixvec_transform_name)
    end select
end subroutine create_mixvec_transform

subroutine create_hor_colocated_mixvec_transform(mixvec_transform,config,domain)
    class(mixvec_transform_t), allocatable, intent(out) :: mixvec_transform
    type(config_mixvec_transform_t),        intent(in)  :: config
    type(domain_t),                         intent(in)  :: domain

    type(mixvec_transform_hor_colocated_t), allocatable :: mixvec

    allocate(mixvec)
    call create_vertical_operator(mixvec%w2p_interp,config%w2p_interp_name)
    call create_vertical_operator(mixvec%p2w_interp,config%p2w_interp_name)
    call create_grid_field(mixvec%wp,0,0,domain%mesh_p)

    call move_alloc(mixvec,mixvec_transform)
end subroutine create_hor_colocated_mixvec_transform

subroutine create_staggered_mixvec_transform(mixvec_transform,config,domain)
    class(mixvec_transform_t), allocatable, intent(out) :: mixvec_transform
    type(config_mixvec_transform_t),        intent(in)  :: config
    type(domain_t),                         intent(in)  :: domain

    type(mixvec_transform_staggered_t), allocatable :: mixvec
    !WORKAROUND:
    integer(kind=4), parameter :: halo_width = 6

    allocate(mixvec)
    call create_vertical_operator(mixvec%w2p_interp,config%w2p_interp_name)
    call create_vertical_operator(mixvec%p2w_interp,config%p2w_interp_name)
    call create_vec2vec_interpolator2d(mixvec%uv2p_interp, config%uv2p_interp_name, domain)
    call create_vec2vec_interpolator2d(mixvec%p2uv_interp, config%p2uv_interp_name, domain)
    call create_grid_field(mixvec%up,halo_width,0,domain%mesh_p)
    call create_grid_field(mixvec%vp,halo_width,0,domain%mesh_p)
    call create_grid_field(mixvec%wp,0,0,domain%mesh_p)

    call move_alloc(mixvec,mixvec_transform)
end subroutine create_staggered_mixvec_transform

end module mixvec_transform_factory_mod
