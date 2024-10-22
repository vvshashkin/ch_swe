module stvec_swm_factory_mod


use domain_mod,                 only : domain_t
use stvec_mod,                  only : stvec_t
use stvec_swm_mod,              only : stvec_swm_t
use stvec_flexible_factory_mod, only : create_stvec_flexible_allocated
use grid_field_factory_mod,     only : create_grid_field
use string_mod,                 only : string_t, strings

implicit none

contains

subroutine create_stvec_swm(stvec, domain, halo_width_xy, halo_width_z)

    class(stvec_t), allocatable, intent(out) :: stvec
    type(domain_t),              intent(in)  :: domain
    integer(kind=4),             intent(in)  :: halo_width_xy, halo_width_z

    type(stvec_swm_t), allocatable :: stvec_swm
    type(string_t) :: field_names(3), mesh_names(3)

    allocate(stvec_swm)

    field_names = strings("u","v","h")
    mesh_names  = strings("u","v","p")

    call create_stvec_flexible_allocated(stvec_swm, field_names, mesh_names, &
                                         halo_width_xy, halo_width_z, domain)

    call stvec_swm%associate_pointers()

    call move_alloc(stvec_swm, stvec)

end subroutine create_stvec_swm

end module stvec_swm_factory_mod
