module operator_swm_factory_mod

use domain_mod,         only : domain_t
use operator_mod,       only : operator_t
use generic_config_mod, only : generic_config_t
use parcomm_mod,        only : parcomm_global
use test_fields_mod,    only : scalar_field_generator_t, zero_scalar_field_generator, &
                               set_scalar_test_field

implicit none

contains

subroutine create_swm_operator(operator, grav, config, domain, orography_generator)

    type(domain_t),                  intent(in)     :: domain
    class(generic_config_t),         intent(inout)  :: config
    real(kind=8),                    intent(in)     :: grav

    class(scalar_field_generator_t), intent(in), optional :: orography_generator

    class(operator_t), allocatable,  intent(out) :: operator

    class(scalar_field_generator_t), allocatable :: orography_generator_loc
    character(len=:), allocatable :: swm_operator_type

    if (present(orography_generator)) then
        orography_generator_loc = orography_generator
    else
        orography_generator_loc = zero_scalar_field_generator
    end if

    call config%get(swm_operator_type, "swm_operator_type")

    select case (swm_operator_type)
    case ("vector_invariant_swm_operator")
        call create_vector_invariant_swm_operator(operator, grav, config, domain, &
                                                  orography_generator_loc)
    case ("vector_invariant_imex_swm_operator")
          call create_vector_invariant_imex_swm_operator(operator, grav, config, domain, &
                                                    orography_generator_loc)
    case ("advective_swm_operator")
        call create_advective_swm_operator(operator, grav, config, domain, &
                                           orography_generator_loc)
    case ("SISL_swm_operator")
        call create_SISL_swm_operator(operator,grav,config,domain,&
                                      orography_generator_loc)
    case ("skew_swm_operator")
        call create_skew_swm_operator(operator,grav,config,domain,&
                                      orography_generator_loc)
    case ("linear_swm_operator")
        call create_linear_swm_operator(operator, grav, config, domain)
    case default
        call parcomm_global%abort("unknown swm_op_type: "//swm_operator_type)
    end select

end subroutine create_swm_operator

subroutine create_vector_invariant_swm_operator(operator, grav, config, domain, &
                                                orography_generator)

    use operator_swm_mod,       only : operator_swm_t
    use div_factory_mod,        only : create_div_operator
    use grad_factory_mod,       only : create_grad_operator
    use curl_factory_mod,       only : create_curl_operator
    use coriolis_factory_mod,   only : create_coriolis
    use KE_factory_mod,         only : create_KE_operator
    use massflux_factory_mod,   only : create_massflux_operator
    use co2contra_factory_mod,  only : create_co2contra_operator
    use quadrature_factory_mod, only : create_quadrature
    use hordiff_factory_mod,    only : create_hordiff_operator

    use grid_field_factory_mod, only : create_grid_field

    type(domain_t),                  intent(in)     :: domain
    class(generic_config_t),         intent(inout)  :: config
    real(kind=8),                    intent(in)     :: grav
    class(scalar_field_generator_t), intent(in)     :: orography_generator

    class(operator_t), allocatable,  intent(out)    :: operator

    type(operator_swm_t), allocatable :: swm_op
    integer(kind=4) :: halo_width_xy
    character(len = :), allocatable :: name

    !WORKAROUND
    halo_width_xy = 8

    allocate(swm_op)

    call config%get(name,"div_op_name")
    swm_op%div_op = create_div_operator(domain, name)

    call config%get(name,"grad_op_name")
    swm_op%grad_op =  create_grad_operator(domain, name)

    call config%get(name,"coriolis_op_name")
    call create_coriolis(swm_op%coriolis_op, name, domain)

    call config%get(name,"curl_op_name")
    call create_curl_operator(swm_op%curl_op, name, domain)

    call config%get(name,"KE_op_name")
    call create_KE_operator(swm_op%KE_op, name, domain)

    call config%get(name,"massflux_op_name")
    swm_op%massflux_op = create_massflux_operator(domain, name)

    call config%get(name,"co2contra_op_name")
    swm_op%co2contra_op = create_co2contra_operator(domain, name)

    call config%get(name,"quadrature_name")
    call create_quadrature(swm_op%quadrature_h, name, domain%mesh_p)
    call create_quadrature(swm_op%quadrature_u, name, domain%mesh_u)
    call create_quadrature(swm_op%quadrature_v, name, domain%mesh_v)
    call create_quadrature(swm_op%quadrature_w, name, domain%mesh_q)

    call create_grid_field(swm_op%KE,  halo_width_xy, 0, domain%mesh_p)
    call create_grid_field(swm_op%div, halo_width_xy, 0, domain%mesh_p)


    call create_grid_field(swm_op%curl, halo_width_xy, 0, domain%mesh_q)

    call create_grid_field(swm_op%hu, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%hv, halo_width_xy, 0, domain%mesh_v)

    call create_grid_field(swm_op%hu_diag, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%hv_diag, halo_width_xy, 0, domain%mesh_v)

    call create_grid_field(swm_op%cor_u, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%cor_v, halo_width_xy, 0, domain%mesh_v)

    call create_grid_field(swm_op%ut, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%vt, halo_width_xy, 0, domain%mesh_v)

    call create_grid_field(swm_op%grad_x, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%grad_y, halo_width_xy, 0, domain%mesh_v)

    call create_grid_field(swm_op%h_surf,    halo_width_xy, 0, domain%mesh_p)
    call set_scalar_test_field(swm_op%h_surf, orography_generator, domain%mesh_p, 0)

    !WORKAROUND
    swm_op%grav = grav

    call create_grid_field(swm_op%KE_diag_u,  0, 0, domain%mesh_u)
    call create_grid_field(swm_op%KE_diag_v,  0, 0, domain%mesh_v)
    call create_grid_field(swm_op%PE_diag,    0, 0, domain%mesh_p)

    call move_alloc(swm_op, operator)

end subroutine create_vector_invariant_swm_operator

subroutine create_vector_invariant_imex_swm_operator(operator, grav, config, domain, &
                                                orography_generator)

    use operator_swm_vec_inv_imex_mod, only : operator_swm_vec_inv_imex_t
    use div_factory_mod,               only : create_div_operator
    use grad_factory_mod,              only : create_grad_operator
    use curl_factory_mod,              only : create_curl_operator
    use coriolis_factory_mod,          only : create_coriolis
    use KE_factory_mod,                only : create_KE_operator
    use massflux_factory_mod,          only : create_massflux_operator
    use co2contra_factory_mod,         only : create_co2contra_operator
    use quadrature_factory_mod,        only : create_quadrature
    use hordiff_factory_mod,           only : create_hordiff_operator

    use stvec_swm_factory_mod, only : create_stvec_swm

    use cg_solver_mod,          only : cg_solver_t
    use bicgstab_solver_mod,    only : bicgstab_solver_t
    use outputer_factory_mod,   only : create_master_paneled_outputer
    use swm_helm_operator_mod,  only : swm_helm_operator_t

    use grid_field_based_vector_mod, only : grid_field_based_vector_t

    use grid_field_factory_mod, only : create_grid_field

    use iterative_solver_factory_mod, only : create_iterative_solver

    type(domain_t),                  intent(in)     :: domain
    class(generic_config_t),         intent(inout)  :: config
    real(kind=8),                    intent(in)     :: grav
    class(scalar_field_generator_t), intent(in)     :: orography_generator

    class(operator_t), allocatable,  intent(out)    :: operator

    type(operator_swm_vec_inv_imex_t), allocatable :: swm_op
    integer(kind=4) :: halo_width_xy
    character(len = :), allocatable :: name
    class(generic_config_t), allocatable :: solver_config

    logical :: is_found

    !WORKAROUND
    halo_width_xy = 8

    allocate(swm_op)

    call config%get(name, "div_op_name")
    swm_op%div_op = create_div_operator(domain, name)

    call config%get(name, "grad_op_name")
    swm_op%grad_op =  create_grad_operator(domain, name)

    call config%get(name, "coriolis_op_name")
    call create_coriolis(swm_op%coriolis_op, name, domain)

    call config%get(name, "curl_op_name")
    call create_curl_operator(swm_op%curl_op, name, domain)

    call config%get(name, "KE_op_name")
    call create_KE_operator(swm_op%KE_op, name, domain)

    call config%get(name, "massflux_op_name")
    swm_op%massflux_op = create_massflux_operator(domain, name)

    call config%get(name, "co2contra_op_name")
    swm_op%co2contra_op = create_co2contra_operator(domain, name)

    call config%get(name, "quadrature_name")
    call create_quadrature(swm_op%quadrature_h, name, domain%mesh_p)
    call create_quadrature(swm_op%quadrature_u, name, domain%mesh_u)
    call create_quadrature(swm_op%quadrature_v, name, domain%mesh_v)
    call create_quadrature(swm_op%quadrature_w, name, domain%mesh_q)

    call create_grid_field(swm_op%KE,  halo_width_xy, 0, domain%mesh_p)
    call create_grid_field(swm_op%div, halo_width_xy, 0, domain%mesh_p)


    call create_grid_field(swm_op%curl, halo_width_xy, 0, domain%mesh_q)

    call create_grid_field(swm_op%hu, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%hv, halo_width_xy, 0, domain%mesh_v)

    call create_grid_field(swm_op%hu_diag, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%hv_diag, halo_width_xy, 0, domain%mesh_v)

    call create_grid_field(swm_op%cor_u, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%cor_v, halo_width_xy, 0, domain%mesh_v)

    call create_grid_field(swm_op%ut, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%vt, halo_width_xy, 0, domain%mesh_v)

    call create_grid_field(swm_op%grad_x, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%grad_y, halo_width_xy, 0, domain%mesh_v)

    call create_grid_field(swm_op%h_surf,    halo_width_xy, 0, domain%mesh_p)
    call set_scalar_test_field(swm_op%h_surf, orography_generator, domain%mesh_p, 0)

    !WORKAROUND
    swm_op%grav = grav

    call create_grid_field(swm_op%KE_diag_u,  0, 0, domain%mesh_u)
    call create_grid_field(swm_op%KE_diag_v,  0, 0, domain%mesh_v)
    call create_grid_field(swm_op%PE_diag,    0, 0, domain%mesh_p)

    call config%get(solver_config, "newton_solver_cfg", found = is_found)
    if (is_found) then
        call solver_config%get(swm_op%Newton_iterations_num,   "iter_num")
        call solver_config%get(swm_op%is_Newton_solver_verbose, "verbose")
    else
        swm_op%Newton_iterations_num = 3
        swm_op%is_Newton_solver_verbose = .false.
    end if

    call config%get(solver_config, "helm_solver_cfg")
    call create_iterative_solver(swm_op%helm_solver, solver_config, domain)

    swm_op%helm_oper%div_op       = swm_op%div_op
    swm_op%helm_oper%grad_op      = swm_op%grad_op
    swm_op%helm_oper%co2contra_op = swm_op%co2contra_op
    swm_op%helm_oper%massflux_op  = swm_op%massflux_op

    call create_grid_field(swm_op%helm_oper%grad_x,      halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%helm_oper%grad_y,      halo_width_xy, 0, domain%mesh_v)
    call create_grid_field(swm_op%helm_oper%grad_x_t,    halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%helm_oper%grad_y_t,    halo_width_xy, 0, domain%mesh_v)
    call create_grid_field(swm_op%helm_oper%gamma_h_ref, halo_width_xy, 0, domain%mesh_p)

    call create_grid_field(swm_op%helm_oper%hu,   halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%helm_oper%hv,   halo_width_xy, 0, domain%mesh_v)

    swm_op%helm_sol%mesh => domain%mesh_p
    swm_op%helm_sol%quadrature = swm_op%quadrature_h
    call create_grid_field(swm_op%helm_sol%grid_field, halo_width_xy, 0, domain%mesh_p)

    swm_op%helm_rhs%mesh => domain%mesh_p
    swm_op%helm_rhs%quadrature = swm_op%quadrature_h
    call create_grid_field(swm_op%helm_rhs%grid_field, halo_width_xy, 0, domain%mesh_p)

    call create_stvec_swm(swm_op%residual, domain, halo_width_xy, 0)
    call create_stvec_swm(swm_op%delta_v,  domain, halo_width_xy, 0)

    call move_alloc(swm_op, operator)

end subroutine create_vector_invariant_imex_swm_operator

subroutine create_advective_swm_operator(operator, grav, config, domain, &
                                         orography_generator)

    use operator_adv_swm_mod,   only : operator_adv_swm_t
    use div_factory_mod,        only : create_div_operator
    use grad_factory_mod,       only : create_grad_operator
    use coriolis_factory_mod,   only : create_coriolis
    use massflux_factory_mod,   only : create_massflux_operator
    use co2contra_factory_mod,  only : create_co2contra_operator
    use quadrature_factory_mod, only : create_quadrature
    use hordiff_factory_mod,    only : create_hordiff_operator

    use grid_field_factory_mod, only : create_grid_field

    use vector_advection_factory_mod, only : create_vector_advection_operator

    type(domain_t),                  intent(in)     :: domain
    class(generic_config_t),         intent(inout)  :: config
    real(kind=8),                    intent(in)     :: grav
    class(operator_t), allocatable,  intent(out)    :: operator
    class(scalar_field_generator_t), intent(in)     :: orography_generator

    type(operator_adv_swm_t), allocatable :: swm_op
    integer(kind=4) :: halo_width_xy
    character(len=:), allocatable :: name

    !WORKAROUND
    halo_width_xy = 8

    allocate(swm_op)

    call config%get(swm_op%v_components_type,"v_components_type")

    call config%get(name,"div_op_name")
    swm_op%div_op = create_div_operator(domain, name)

    call config%get(name,"grad_op_name")
    swm_op%grad_op =  create_grad_operator(domain, name)

    call config%get(name,"coriolis_op_name")
    call create_coriolis(swm_op%coriolis_op, name, domain)

    call config%get(name,"massflux_op_name")
    swm_op%massflux_op = create_massflux_operator(domain, name)

    call config%get(name,"co2contra_op_name")
    swm_op%co2contra_op = create_co2contra_operator(domain, name)

    call config%get(name,"vector_advection_op_name")
    call create_vector_advection_operator(swm_op%adv_uv, name, domain)

    call config%get(name,"quadrature_name")
    call create_quadrature(swm_op%quadrature_h, name, domain%mesh_p)
    call create_quadrature(swm_op%quadrature_u, name, domain%mesh_u)
    call create_quadrature(swm_op%quadrature_v, name, domain%mesh_v)

    call create_grid_field(swm_op%div, halo_width_xy, 0, domain%mesh_p)
    call create_grid_field(swm_op%h_total, halo_width_xy, 0, domain%mesh_p)
    call create_grid_field(swm_op%hu, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%hv, halo_width_xy, 0, domain%mesh_v)
    call create_grid_field(swm_op%cor_u, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%cor_v, halo_width_xy, 0, domain%mesh_v)
    call create_grid_field(swm_op%ut, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%vt, halo_width_xy, 0, domain%mesh_v)
    call create_grid_field(swm_op%grad_x, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%grad_y, halo_width_xy, 0, domain%mesh_v)

    call create_grid_field(swm_op%h_surf,    halo_width_xy, 0, domain%mesh_p)
    call set_scalar_test_field(swm_op%h_surf, orography_generator, domain%mesh_p, 0)

    !WORKAROUND !Why workarond stays here?
    swm_op%grav = grav

    call create_grid_field(swm_op%KE_diag_u,  0, 0, domain%mesh_u)
    call create_grid_field(swm_op%KE_diag_v,  0, 0, domain%mesh_v)
    call create_grid_field(swm_op%PE_diag,    0, 0, domain%mesh_p)

    call move_alloc(swm_op, operator)

end subroutine create_advective_swm_operator

subroutine create_SISL_swm_operator(operator, grav, config, domain, &
    orography_generator)

use operator_swm_SISL_mod,  only : operator_swm_SISL_t
use div_factory_mod,        only : create_div_operator
use grad_factory_mod,       only : create_grad_operator
use curl_factory_mod,       only : create_curl_operator
use coriolis_factory_mod,   only : create_coriolis
use KE_factory_mod,         only : create_KE_operator
use massflux_factory_mod,   only : create_massflux_operator
use co2contra_factory_mod,  only : create_co2contra_operator
use quadrature_factory_mod, only : create_quadrature
use hordiff_factory_mod,    only : create_hordiff_operator
use grid_field_factory_mod, only : create_grid_field

use swm_helm_operator_mod,         only : swm_helm_lin_operator_t
use cg_solver_mod,                 only : cg_solver_t
use bicgstab_solver_mod,           only : bicgstab_solver_t
use iterative_solver_factory_mod,  only : create_iterative_solver
use grid_field_based_vector_mod,   only : grid_field_based_vector_t

use dep_points_interp_driver_factory_mod, only : create_dep_points_interp_driver
use trajectory_solver_factory_mod,        only : create_trajectory_solver
use string_mod,                           only : strings
use stvec_flexible_factory_mod,           only : create_stvec_flexible_allocated

type(domain_t),                  intent(in)     :: domain
class(generic_config_t),         intent(inout)  :: config
real(kind=8),                    intent(in)     :: grav
class(scalar_field_generator_t), intent(in)     :: orography_generator

class(operator_t), allocatable,  intent(out)    :: operator

type(operator_swm_SISL_t), allocatable :: swm_op
integer(kind=4) :: halo_width_xy
character(len = :), allocatable :: name
class(generic_config_t), allocatable :: dp_driver_config, trajectory_solver_config, &
                                        lin_solver_config

!WORKAROUND
halo_width_xy = 8

allocate(swm_op)

call config%get(name,"div_op_name")
swm_op%div_op = create_div_operator(domain, name)

call config%get(name,"grad_op_name")
swm_op%grad_op =  create_grad_operator(domain, name)

call config%get(name,"coriolis_op_name")
call create_coriolis(swm_op%coriolis_op, name, domain)

call config%get(name,"co2contra_op_name")
swm_op%co2contra_op = create_co2contra_operator(domain, name)

call config%get(name,"quadrature_name")
call create_quadrature(swm_op%quadrature_h, name, domain%mesh_p)
call create_quadrature(swm_op%quadrature_u, name, domain%mesh_u)
call create_quadrature(swm_op%quadrature_v, name, domain%mesh_v)
call create_quadrature(swm_op%quadrature_w, name, domain%mesh_q)

call create_grid_field(swm_op%h_tot,  halo_width_xy, 0, domain%mesh_p)
call create_grid_field(swm_op%div, halo_width_xy, 0, domain%mesh_p)

call create_grid_field(swm_op%cor_u, halo_width_xy, 0, domain%mesh_u)
call create_grid_field(swm_op%cor_v, halo_width_xy, 0, domain%mesh_v)

call create_grid_field(swm_op%grad_x, halo_width_xy, 0, domain%mesh_u)
call create_grid_field(swm_op%grad_y, halo_width_xy, 0, domain%mesh_v)

call create_grid_field(swm_op%h_surf,    halo_width_xy, 0, domain%mesh_p)
call set_scalar_test_field(swm_op%h_surf, orography_generator, domain%mesh_p, 0)

!WORKAROUND
swm_op%grav = grav

call create_grid_field(swm_op%KE_diag_u,  0, 0, domain%mesh_u)
call create_grid_field(swm_op%KE_diag_v,  0, 0, domain%mesh_v)
call create_grid_field(swm_op%PE_diag,    0, 0, domain%mesh_p)
call create_grid_field(swm_op%hu_diag, 0, 0, domain%mesh_u)
call create_grid_field(swm_op%hv_diag, 0, 0, domain%mesh_v)

call config%get(lin_solver_config, "helm_solver_cfg")
call create_iterative_solver(swm_op%helm_solver, lin_solver_config, domain)

!Helm solver

call config%get(swm_op%H0,"H0")

swm_op%helm_oper%div_op       = swm_op%div_op
swm_op%helm_oper%grad_op      = swm_op%grad_op
swm_op%helm_oper%co2contra_op = swm_op%co2contra_op

call create_grid_field(swm_op%helm_oper%grad_x,      halo_width_xy, 0, domain%mesh_u)
call create_grid_field(swm_op%helm_oper%grad_y,      halo_width_xy, 0, domain%mesh_v)
call create_grid_field(swm_op%helm_oper%grad_x_t,    halo_width_xy, 0, domain%mesh_u)
call create_grid_field(swm_op%helm_oper%grad_y_t,    halo_width_xy, 0, domain%mesh_v)
call create_grid_field(swm_op%helm_oper%gamma_h_ref,             0, 0, domain%mesh_p)

swm_op%helm_sol%mesh => domain%mesh_p
swm_op%helm_sol%quadrature = swm_op%quadrature_h
call create_grid_field(swm_op%helm_sol%grid_field, halo_width_xy, 0, domain%mesh_p)

swm_op%helm_rhs%mesh => domain%mesh_p
swm_op%helm_rhs%quadrature = swm_op%quadrature_h
call create_grid_field(swm_op%helm_rhs%grid_field, halo_width_xy, 0, domain%mesh_p)

!semi-Lagrangian
call config%get(dp_driver_config, "dep_points_interp_driver")
call create_dep_points_interp_driver(swm_op%dp_interp_driver, domain, &
                                     dp_driver_config)

call config%get(trajectory_solver_config, "trajectory_solver")
call create_trajectory_solver(swm_op%trajectory_solver, trajectory_solver_config, domain)

call create_grid_field(swm_op%vx, 0, 0, domain%mesh_p)
call create_grid_field(swm_op%vy, 0, 0, domain%mesh_p)
call create_grid_field(swm_op%vz, 0, 0, domain%mesh_p)

call create_stvec_flexible_allocated(swm_op%dp_coords,strings("x","y","z","alpha","beta","panel_ind"),strings("p","p","p","p","p","p"),0,0,domain)

call move_alloc(swm_op, operator)

end subroutine create_SISL_swm_operator

subroutine create_skew_swm_operator(operator, grav, config, domain, &
                                    orography_generator)

    use operator_skew_swm_mod,      only : operator_skew_swm_t
    use div_factory_mod,            only : create_div_operator
    use grad_factory_mod,           only : create_grad_operator
    use coriolis_factory_mod,       only : create_coriolis
    use massflux_factory_mod,       only : create_massflux_operator
    use co2contra_factory_mod,      only : create_co2contra_operator
    use quadrature_factory_mod,     only : create_quadrature
    use flux_div_factory_mod,       only : create_flux_div_operator
    use interpolator2d_factory_mod, only : create_vec2vec_interpolator2d
    use grid_field_factory_mod,     only : create_grid_field

    type(domain_t),                  intent(in)     :: domain
    class(generic_config_t),         intent(inout)  :: config
    real(kind=8),                    intent(in)     :: grav
    class(operator_t), allocatable,  intent(out)    :: operator
    class(scalar_field_generator_t), intent(in)     :: orography_generator

    type(operator_skew_swm_t), allocatable :: swm_op
    integer(kind=4) :: halo_width_xy
    character(len=:), allocatable :: name

    !WORKAROUND
    halo_width_xy = 8

    allocate(swm_op)

    call config%get(swm_op%v_components_type,"v_components_type")

    call config%get(name,"h_flux_div_name")
    call create_flux_div_operator(swm_op%h_flux_div, name, domain)

    call config%get(name,"div_op_name")
    swm_op%div_op = create_div_operator(domain, name)

    call config%get(name,"grad_op_name")
    swm_op%grad_op =  create_grad_operator(domain, name)

    call config%get(name,"coriolis_op_name")
    call create_coriolis(swm_op%coriolis_op, name, domain)

    ! call config%get(name,"massflux_op_name")
    ! swm_op%massflux_op = create_massflux_operator(domain, name)

    call config%get(name,"co2contra_op_name")
    swm_op%co2contra_op = create_co2contra_operator(domain, name)

    call create_grid_field(swm_op%hx, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%hy, halo_width_xy, 0, domain%mesh_v)

    call create_grid_field(swm_op%up, halo_width_xy, 0, domain%mesh_p)
    call create_grid_field(swm_op%vp, halo_width_xy, 0, domain%mesh_p)

    call create_grid_field(swm_op%vx, halo_width_xy, 0, domain%mesh_p)
    call create_grid_field(swm_op%vy, halo_width_xy, 0, domain%mesh_p)
    call create_grid_field(swm_op%vz, halo_width_xy, 0, domain%mesh_p)
    call create_grid_field(swm_op%vx_tend_adv, halo_width_xy, 0, domain%mesh_p)
    call create_grid_field(swm_op%vy_tend_adv, halo_width_xy, 0, domain%mesh_p)
    call create_grid_field(swm_op%vz_tend_adv, halo_width_xy, 0, domain%mesh_p)
    call create_grid_field(swm_op%vx_tend_metric, halo_width_xy, 0, domain%mesh_p)
    call create_grid_field(swm_op%vy_tend_metric, halo_width_xy, 0, domain%mesh_p)
    call create_grid_field(swm_op%vz_tend_metric, halo_width_xy, 0, domain%mesh_p)
    call create_grid_field(swm_op%div, halo_width_xy, 0, domain%mesh_p)
    call create_grid_field(swm_op%p_buff, halo_width_xy, 0, domain%mesh_p)

    call config%get(name, "v2p_interp_name")
    call create_vec2vec_interpolator2d(swm_op%v2p_interp, name, domain)
    call config%get(name, "p2v_interp_name")
    call create_vec2vec_interpolator2d(swm_op%p2v_interp, name, domain)

    ! call create_grid_field(swm_op%h_surf,    halo_width_xy, 0, domain%mesh_p)
    ! call set_scalar_test_field(swm_op%h_surf, orography_generator, domain%mesh_p, 0)

    swm_op%grav = grav

    call move_alloc(swm_op, operator)

end subroutine create_skew_swm_operator

subroutine create_linear_swm_operator(operator, grav, config, domain)

    use operator_swm_lin_mod,   only : operator_swm_lin_t
    use div_factory_mod,        only : create_div_operator
    use grad_factory_mod,       only : create_grad_operator
    use curl_factory_mod,       only : create_curl_operator
    use coriolis_factory_mod,   only : create_coriolis
    use KE_factory_mod,         only : create_KE_operator
    use massflux_factory_mod,   only : create_massflux_operator
    use co2contra_factory_mod,  only : create_co2contra_operator
    use quadrature_factory_mod, only : create_quadrature
    use hordiff_factory_mod,    only : create_hordiff_operator
    use cg_solver_mod,          only : cg_solver_t
    use bicgstab_solver_mod,    only : bicgstab_solver_t
    use outputer_factory_mod,   only : create_master_paneled_outputer
    use swm_helm_operator_mod,  only : swm_helm_lin_operator_t

    use iterative_solver_factory_mod, only : create_iterative_solver

    use grid_field_based_vector_mod, only : grid_field_based_vector_t

    use grid_field_factory_mod, only : create_grid_field

    type(domain_t),           intent(in)    :: domain
    class(generic_config_t),  intent(inout) :: config
    real(kind=8),             intent(in)    :: grav

    class(operator_t), allocatable, intent(out) :: operator

    type(operator_swm_lin_t),        allocatable :: swm_op
    type(swm_helm_lin_operator_t),   allocatable :: helm_oper

    class(generic_config_t), allocatable :: solver_config
    character(len=:), allocatable :: name
    integer(kind=4) :: halo_width_xy, t

    !WORKAROUND
    halo_width_xy = 8

    allocate(swm_op)

    call config%get(swm_op%h_ref, "h_mean")

    call config%get(name, "div_op_name")
    swm_op%div_op = create_div_operator(domain, name)

    call config%get(name, "grad_op_name")
    swm_op%grad_op = create_grad_operator(domain, name)

    call config%get(name,"co2contra_op_name")
    swm_op%co2contra_op = create_co2contra_operator(domain, name)

    call config%get(name,"coriolis_op_name")
    call create_coriolis(swm_op%coriolis_op, name, domain)

    call config%get(name,"quadrature_name")
    call create_quadrature(swm_op%quadrature_h, name, domain%mesh_p)
    call create_quadrature(swm_op%quadrature_u, name, domain%mesh_u)
    call create_quadrature(swm_op%quadrature_v, name, domain%mesh_v)
    call create_quadrature(swm_op%quadrature_w, name, domain%mesh_q)


    call create_grid_field(swm_op%div, halo_width_xy, 0, domain%mesh_p)

    call create_grid_field(swm_op%ut, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%vt, halo_width_xy, 0, domain%mesh_v)

    call create_grid_field(swm_op%grad_x, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%grad_y, halo_width_xy, 0, domain%mesh_v)

    call config%get(solver_config, "helm_solver_cfg")
    call create_iterative_solver(swm_op%helm_solver, solver_config, domain)


    swm_op%helm_oper%div_op       = swm_op%div_op
    swm_op%helm_oper%grad_op      = swm_op%grad_op
    swm_op%helm_oper%co2contra_op = swm_op%co2contra_op

    call create_grid_field(swm_op%helm_oper%grad_x,      halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%helm_oper%grad_y,      halo_width_xy, 0, domain%mesh_v)
    call create_grid_field(swm_op%helm_oper%grad_x_t,    halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%helm_oper%grad_y_t,    halo_width_xy, 0, domain%mesh_v)
    call create_grid_field(swm_op%helm_oper%gamma_h_ref,             0, 0, domain%mesh_p)

    swm_op%helm_sol%mesh => domain%mesh_p
    swm_op%helm_sol%quadrature = swm_op%quadrature_h
    call create_grid_field(swm_op%helm_sol%grid_field, halo_width_xy, 0, domain%mesh_p)

    swm_op%helm_rhs%mesh => domain%mesh_p
    swm_op%helm_rhs%quadrature = swm_op%quadrature_h
    call create_grid_field(swm_op%helm_rhs%grid_field, halo_width_xy, 0, domain%mesh_p)

    call create_master_paneled_outputer(swm_op%outputer, "p", domain)

    !WORKAROUND
    swm_op%grav = grav

    call move_alloc(swm_op, operator)

end subroutine create_linear_swm_operator

end module operator_swm_factory_mod
