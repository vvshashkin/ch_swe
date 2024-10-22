module swm_helm_operator_mod

use linear_operator_mod,         only : linear_operator_t
use abstract_grad_mod,           only : grad_operator_t
use abstract_div_mod,            only : div_operator_t
use abstract_co2contra_mod,      only : co2contra_operator_t
use domain_mod,                  only : domain_t
use grid_field_mod,              only : grid_field_t
use vector_mod,                  only : abstract_vector_t
use parcomm_mod,                 only : parcomm_global
use grid_field_based_vector_mod, only : grid_field_based_vector_t
use abstract_massflux_mod,       only : massflux_operator_t

implicit none

type, public, extends(linear_operator_t) :: swm_helm_lin_operator_t

    ! (I - gamma_h_ref * div(grad)) operator
    ! gamma_h_ref - spatially varying field usually
    ! scaled by timescheme specific const gamma

    class(div_operator_t),       allocatable :: div_op
    class(grad_operator_t),      allocatable :: grad_op
    class(co2contra_operator_t), allocatable :: co2contra_op

    type(grid_field_t) :: grad_x, grad_y, grad_x_t, grad_y_t
    type(grid_field_t) :: gamma_h_ref

contains
    procedure, public :: apply => apply_swm_helm_lin_operator
end type swm_helm_lin_operator_t


type, public, extends(linear_operator_t) :: swm_helm_operator_t

    ! (I - div(gamma_h_ref * grad)) operator
    ! gamma_h_ref - spatially varying field usually
    ! scaled by timescheme specific const gamma

    class(div_operator_t),       allocatable :: div_op
    class(grad_operator_t),      allocatable :: grad_op
    class(co2contra_operator_t), allocatable :: co2contra_op
    class(massflux_operator_t),  allocatable :: massflux_op

    type(grid_field_t) :: grad_x, grad_y, grad_x_t, grad_y_t
    type(grid_field_t) :: hu, hv !mass fluxes
    type(grid_field_t) :: gamma_h_ref


contains
    procedure, public :: apply => apply_swm_helm_operator
end type swm_helm_operator_t

contains

subroutine apply_swm_helm_lin_operator(this, vout, vin, domain)
    class(swm_helm_lin_operator_t), intent(inout) :: this
    class(abstract_vector_t),       intent(inout) :: vout, vin
    type(domain_t),                 intent(in)    :: domain

    select type (vin)
    class is (grid_field_based_vector_t)
        select type (vout)
        class is (grid_field_based_vector_t)
            call this%grad_op%calc_grad(this%grad_x, this%grad_y, vin%grid_field, domain)

            call this%co2contra_op%transform(this%grad_x_t, this%grad_y_t, &
                                     this%grad_x, this%grad_y, domain)

            call this%div_op%calc_div(vout%grid_field, this%grad_x_t, this%grad_y_t, domain)

            call vout%grid_field%assign_prod(-1.0_8, this%gamma_h_ref, vout%grid_field, vout%mesh)

            call vout%grid_field%update(1.0_8, vin%grid_field, vout%mesh)
        class default
            call parcomm_global%abort("swm_helm_operator_t error: vout of wrong type")
        end select
    class default
        call parcomm_global%abort("swm_helm_operator_t error: vin of wrong type")
    end select
end subroutine apply_swm_helm_lin_operator

subroutine apply_swm_helm_operator(this, vout, vin, domain)
    class(swm_helm_operator_t), intent(inout) :: this
    class(abstract_vector_t),   intent(inout) :: vout, vin
    type(domain_t),             intent(in)    :: domain

    select type (vin)
    class is (grid_field_based_vector_t)
        select type (vout)
        class is (grid_field_based_vector_t)
            call this%grad_op%calc_grad(this%grad_x, this%grad_y, vin%grid_field, domain)

            call this%co2contra_op%transform(this%grad_x_t, this%grad_y_t, &
                                     this%grad_x, this%grad_y, domain)

            call this%massflux_op%calc_massflux(this%hu, this%hv, this%gamma_h_ref, &
                                     this%grad_x_t, this%grad_y_t, domain)

            call this%div_op%calc_div(vout%grid_field, this%hu, this%hv, domain)

            call vout%grid_field%assign(1.0_8, vin%grid_field, -1.0_8,  vout%grid_field, vout%mesh)
        class default
            call parcomm_global%abort("swm_helm_operator_t error: vout of wrong type")
        end select
    class default
        call parcomm_global%abort("swm_helm_operator_t error: vin of wrong type")
    end select
end subroutine apply_swm_helm_operator

end module swm_helm_operator_mod
