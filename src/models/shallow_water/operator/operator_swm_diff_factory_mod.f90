module operator_swm_diff_factory_mod

use domain_mod,         only : domain_t
use operator_mod,       only : operator_t
use generic_config_mod, only : generic_config_t

implicit none

contains

subroutine create_swm_diff_operator(operator, config, domain)

    use operator_swm_diff_mod,  only : operator_swm_diff_t
    use hordiff_factory_mod,    only : create_hordiff_operator

    type(domain_t),                 intent(in)     :: domain
    class(generic_config_t),        intent(inout)  :: config
    type(operator_swm_diff_t),      allocatable, intent(out) :: operator

    character(len=:), allocatable :: name
    real(kind=8) :: coef

    allocate(operator)

    call config%get(name,"hordiff_uv_name")
    call config%get(coef,"uv_diff_coeff")
    call create_hordiff_operator(operator%hordiff_uv, name, coef, domain)

    call config%get(name,"hordiff_h_name")
    call config%get(coef,"h_diff_coeff")
    call create_hordiff_operator(operator%hordiff, name, coef, domain)

end subroutine create_swm_diff_operator

end module operator_swm_diff_factory_mod
