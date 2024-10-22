module abstract_mixvec_transform_mod

use grid_field_mod, only : grid_field_t
use domain_mod,     only : domain_t
use parcomm_mod,    only : parcomm_global

implicit none

type, abstract, public :: mixvec_transform_t
    contains
        procedure(transform_co2mix),        deferred :: transform_co2mix
        procedure(transform_mix2contra),    deferred :: transform_mix2contra
end type mixvec_transform_t

abstract interface
    subroutine transform_co2mix(this, u_contra, v_contra, w, u_cov, v_cov, w_cov, domain)
        import mixvec_transform_t, grid_field_t, domain_t
        class(mixvec_transform_t), intent(inout) :: this
        type(domain_t),            intent(in)    :: domain
        type(grid_field_t),        intent(inout) :: u_cov, v_cov, w_cov
        !output:
        type(grid_field_t),        intent(inout) :: u_contra, v_contra, w
    end subroutine transform_co2mix
    subroutine transform_mix2contra(this, w_contra, u_contra, v_contra, w, domain)
        import mixvec_transform_t, grid_field_t, domain_t
        class(mixvec_transform_t), intent(inout) :: this
        type(grid_field_t),        intent(inout) :: u_contra, v_contra, w
        type(domain_t),            intent(in)    :: domain
        !output:
        type(grid_field_t),        intent(inout) :: w_contra
    end subroutine transform_mix2contra
end interface

contains

end module abstract_mixvec_transform_mod
