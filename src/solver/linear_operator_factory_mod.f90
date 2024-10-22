module linear_operator_factory_mod

use domain_mod,          only : domain_t
use linear_operator_mod, only : linear_operator_t
use generic_config_mod,  only : generic_config_t

implicit none

contains

subroutine create_linear_operator(linear_operator, domain, config)

    class(linear_operator_t), allocatable, intent(out)   :: linear_operator
    type(domain_t),                        intent(in)    :: domain
    class(generic_config_t),               intent(inout) :: config

    character(len=:), allocatable :: linear_operator_name


    call config%get(linear_operator_name, "linear_operator_name")

    select case (linear_operator_name)
    case ("swm_helm_lin_operator")
        call create_swm_helm_lin_operator(linear_operator, domain, config)
    case default
        call domain%parcomm%abort("unknown linear_operator_name: "//linear_operator_name)
    end select

end subroutine create_linear_operator

subroutine create_swm_helm_lin_operator(linear_operator, domain, config)

    use swm_helm_operator_mod,  only : swm_helm_lin_operator_t
    use grid_field_factory_mod, only : create_grid_field
    use div_factory_mod,        only : create_div_operator
    use grad_factory_mod,       only : create_grad_operator
    use co2contra_factory_mod,  only : create_co2contra_operator

    class(linear_operator_t), allocatable, intent(out)   :: linear_operator
    type(domain_t),                        intent(in)    :: domain
    class(generic_config_t),               intent(inout) :: config

    type(swm_helm_lin_operator_t), allocatable :: helm_oper

    integer(kind=4) :: hw_xy
    character(len=:), allocatable :: name
    real(kind=8) :: helm_const

    allocate(helm_oper)

    call config%get(name, "div_op_name")
    helm_oper%div_op = create_div_operator(domain, name)

    call config%get(name, "grad_op_name")
    helm_oper%grad_op = create_grad_operator(domain, name)

    call config%get(name,"co2contra_op_name")
    helm_oper%co2contra_op = create_co2contra_operator(domain, name)

    call config%get(hw_xy,"halo_width_xy")

    call create_grid_field(helm_oper%grad_x,      hw_xy, 0, domain%mesh_u)
    call create_grid_field(helm_oper%grad_y,      hw_xy, 0, domain%mesh_v)
    call create_grid_field(helm_oper%grad_x_t,    hw_xy, 0, domain%mesh_u)
    call create_grid_field(helm_oper%grad_y_t,    hw_xy, 0, domain%mesh_v)
    call create_grid_field(helm_oper%gamma_h_ref, hw_xy, 0, domain%mesh_p)

    call config%get(helm_const, "helm_const")

    call helm_oper%gamma_h_ref%assign(helm_const, domain%mesh_p)

    call move_alloc(helm_oper, linear_operator)

end subroutine create_swm_helm_lin_operator

end module linear_operator_factory_mod
