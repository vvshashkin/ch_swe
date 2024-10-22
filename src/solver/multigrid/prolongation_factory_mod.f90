module prolongation_factory_mod

use abstract_prolongation_mod, only : prolongation_operator_t
use domain_mod,                only : domain_t
use parcomm_mod,               only : parcomm_global

implicit none

contains

subroutine create_prolongation_operator(prolongation_op, prolongation_name, domain_coarse)

    class(prolongation_operator_t), allocatable, intent(out)   :: prolongation_op
    type(domain_t),                              intent(inout) :: domain_coarse
    character(len=*),                            intent(in)    :: prolongation_name

    select case (prolongation_name)
    case ("bilinear_const")
        call create_bilinear_const_prolongation(prolongation_op, domain_coarse)
    case default
        call parcomm_global%abort("Unknown prolongation operator name "//prolongation_name)
    end select

end subroutine create_prolongation_operator

subroutine create_bilinear_const_prolongation(prolongation_op, domain_coarse)

    use prolongation_bilinear_const_mod, only : prolongation_bilinear_const_t
    use exchange_factory_mod,            only : create_o_points_halo_exchange

    class(prolongation_operator_t), allocatable, intent(out)   :: prolongation_op
    type(domain_t),                              intent(inout) :: domain_coarse

    type(prolongation_bilinear_const_t), allocatable :: prolongation_bilinear_const

    integer(kind=4) :: halo_width

    halo_width = 1

    allocate(prolongation_bilinear_const)

    prolongation_bilinear_const%halo_exchange = create_o_points_halo_exchange( &
                                                domain_coarse%partition,  &
                                                domain_coarse%parcomm,    &
                                                domain_coarse%topology,   &
                                                halo_width, "full")

    call move_alloc(prolongation_bilinear_const, prolongation_op)

end subroutine create_bilinear_const_prolongation


end module prolongation_factory_mod
