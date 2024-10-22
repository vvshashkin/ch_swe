module restriction_factory_mod

use abstract_restriction_mod, only : restriction_operator_t
use domain_mod,               only : domain_t
use parcomm_mod,              only : parcomm_global

implicit none

contains

subroutine create_restriction_operator(restriction_op, restriction_name, domain_fine)

    class(restriction_operator_t), allocatable, intent(out)   :: restriction_op
    type(domain_t),                             intent(inout) :: domain_fine
    character(len=*),                           intent(in)    :: restriction_name

    select case (restriction_name)
    case ("4point_cell_average")
        call create_4point_cell_average_restriction(restriction_op, domain_fine)
    case default
        call parcomm_global%abort("Unknown restriction operator name "//restriction_name)
    end select

end subroutine create_restriction_operator


subroutine create_4point_cell_average_restriction(restriction_op, domain_fine)

    use restriction_4point_cell_average_mod, only : restriction_4point_cell_average_t
    use exchange_factory_mod,                only : create_o_points_halo_exchange

    class(restriction_operator_t), allocatable, intent(out)   :: restriction_op
    type(domain_t),                             intent(inout) :: domain_fine

    type(restriction_4point_cell_average_t), allocatable :: restriction_4point

    integer(kind=4) :: halo_width

    halo_width = 1

    allocate(restriction_4point)

    restriction_4point%halo_exchange = create_o_points_halo_exchange( &
                                              domain_fine%partition,  &
                                              domain_fine%parcomm,    &
                                              domain_fine%topology,   &
                                              halo_width, "full")

    call move_alloc(restriction_4point, restriction_op)

end subroutine create_4point_cell_average_restriction

end module restriction_factory_mod
