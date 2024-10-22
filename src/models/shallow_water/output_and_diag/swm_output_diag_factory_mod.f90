module swm_output_diag_factory_mod

use swm_output_diag_mod,        only : swm_output_diag_t
use domain_mod,                 only : domain_t
use outputer_factory_mod,       only : create_master_paneled_outputer,&
                                       create_latlon_outputer, create_latlon_vec_outputer
use generic_config_mod,         only : generic_config_t

use curl_factory_mod,           only : create_curl_operator
use co2contra_factory_mod,      only : create_co2contra_operator
use div_factory_mod,            only : create_div_operator
use interpolator2d_factory_mod, only : create_scalar2vec_interpolator2d, &
                                       create_vec2vec_interpolator2d
use grid_field_factory_mod,     only : create_grid_field
use quadrature_factory_mod,     only : create_quadrature
use test_fields_mod,            only : scalar_field_generator_t, set_scalar_test_field
use coriolis_factory_mod,       only : calc_coriolis_parameter

implicit none

!WORKAROUND
integer(kind=4), parameter :: halo_width = 7

contains

subroutine create_swm_output_diag(output_diag, config, orography_generator, &
                                  v_components_type, domain)

    type(swm_output_diag_t),         intent(out)   :: output_diag
    class(generic_config_t),         intent(inout) :: config
    class(scalar_field_generator_t), intent(in)    :: orography_generator
    character(len=*),                intent(in)    :: v_components_type
    type(domain_t),                  intent(in)    :: domain

    integer(kind=4) :: Nlon, Nlat
    character(len=:), allocatable :: u_points, v_points, h_points, curl_points
    character(len=:), allocatable :: op_name

    Nlon = 4*domain%partition%Nh
    Nlat = 2*domain%partition%Nh+1

    select case(domain%horizontal_staggering)
    case("Ah")
        h_points    = "xy"
        u_points    = "xy"
        v_points    = "xy"
        curl_points = "xy"
    case("Ch")
        h_points    = "xy"
        u_points    = "y"
        v_points    = "x"
        curl_points = "o"
    case("C")
        h_points    = "o"
        u_points    = "x"
        v_points    = "y"
        curl_points = "xy"
    case default
        call domain%parcomm%abort("This staggering is not implemented create_swm_output_diag"//&
                                   domain%horizontal_staggering)
    end select

    output_diag%v_components_type = v_components_type

    call config%get(output_diag%output_div, "output_div", default=.false.)
    call config%get(output_diag%output_curl,"output_curl",default=.false.)

    call create_latlon_outputer(output_diag%outputer_h, Nlat, Nlon, &
                                h_points, domain)
    call create_latlon_vec_outputer(output_diag%outputer_v, Nlat, Nlon, &
                                    u_points, v_points, v_components_type, &
                                    domain)
 
    if(output_diag%output_curl) &
        call create_latlon_outputer(output_diag%outputer_curl, Nlat, Nlon, &
                                    curl_points, domain)

    call create_grid_field(output_diag%curl, 1, 0, domain%mesh_q)
    call create_grid_field(output_diag%u, halo_width, 0, domain%mesh_u)
    call create_grid_field(output_diag%v, halo_width, 0, domain%mesh_v)
    call create_grid_field(output_diag%h, 0,0, domain%mesh_p)
    call create_grid_field(output_diag%hu, halo_width,0, domain%mesh_u)
    call create_grid_field(output_diag%hv, halo_width,0, domain%mesh_v)
    call create_grid_field(output_diag%hqu, 0,0, domain%mesh_q)
    call create_grid_field(output_diag%hqv, 0,0, domain%mesh_q)
    call create_grid_field(output_diag%hq, 0,0, domain%mesh_q)
    call create_grid_field(output_diag%hq_tend, 0,0, domain%mesh_q)
    call create_grid_field(output_diag%curl_tend, 1,0, domain%mesh_q)
    call create_grid_field(output_diag%fcori, 0, 0, domain%mesh_q)
    call create_grid_field(output_diag%h_surf, 0, 0, domain%mesh_p)
    call create_grid_field(output_diag%div, 1, 0, domain%mesh_p)

    call set_scalar_test_field(output_diag%h_surf, orography_generator, domain%mesh_p, 0)

    call config%get(op_name,"curl_operator_name",default="None")
    if(op_name /= "None") &
        call create_curl_operator(output_diag%curl_op, op_name, domain)

    call config%get(op_name,"div_operator_name",default="None")
    if(op_name /= "None" .and. output_diag%output_div) &
        output_diag%div_op = create_div_operator(domain, op_name)

    call config%get(op_name,"co2contra_operator_name",default="None")
    if(op_name /= "None") &
         output_diag%co2contra_op = create_co2contra_operator(domain, op_name)

    call config%get(op_name,"h2uv_operator_name",default="None")
    if(op_name /= "None") &
        call create_scalar2vec_interpolator2d(output_diag%interp_h2uv, op_name, domain)

    call config%get(op_name,"uv2q_operator_name",default="None")
    if(op_name /= "None") &
        call create_vec2vec_interpolator2d(output_diag%interp_uv2q, op_name, domain)

    call calc_coriolis_parameter(output_diag%fcori, domain%mesh_q, domain%metric)

    call config%get(op_name,"quadrature_name")
    call create_quadrature(output_diag%quadrature_h, op_name, domain%mesh_p)
    call create_quadrature(output_diag%quadrature_u, op_name, domain%mesh_u)
    call create_quadrature(output_diag%quadrature_v, op_name, domain%mesh_v)
    call create_quadrature(output_diag%quadrature_q, op_name, domain%mesh_q)

end subroutine

end module swm_output_diag_factory_mod