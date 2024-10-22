module massflux_3d_factory_mod

use domain_mod,                    only : domain_t
use abstract_massflux_3d_mod,      only : massflux_3d_operator_t
use massflux_factory_mod,          only : create_massflux_operator
use massflux_vertical_factory_mod, only : create_vertical_massflux
use massflux_3d_mod,               only : massflux_horvert_t
use generic_config_mod,            only : generic_config_t
use string_mod,                    only : string_t

implicit none

contains

subroutine create_massflux_3d_operator(massflux,config,domain)

    !input:
    class(generic_config_t), intent(inout) :: config
    type(domain_t),          intent(in)    :: domain

    !output:
    class(massflux_3d_operator_t), allocatable, intent(out) :: massflux

    !local
    type(massflux_horvert_t), allocatable :: massflux_horvert
    type(string_t) :: massflux_hor_name, massflux_ver_name

    call config%get(massflux_hor_name,"massflux_hor_name")
    call config%get(massflux_ver_name,"massflux_ver_name")

    allocate(massflux_horvert)

    massflux_horvert%massflux_hor = create_massflux_operator(domain,massflux_hor_name%str)
    call create_vertical_massflux(massflux_horvert%massflux_vert,massflux_ver_name%str)

    call move_alloc(massflux_horvert, massflux)

end subroutine

end module
