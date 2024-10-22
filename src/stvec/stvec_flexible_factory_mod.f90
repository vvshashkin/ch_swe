module stvec_flexible_factory_mod

implicit none

contains

subroutine create_stvec_flexible(stvec, field_names, mesh_names, &
                                 halo_width_xy, halo_width_z, domain)

    use stvec_flexible_mod,     only : stvec_flexible_t
    use domain_mod,             only : domain_t
    use stvec_mod,              only : stvec_t
    use string_mod,             only : string_t
    use grid_field_factory_mod, only : create_grid_field
    use grid_field_mod,         only : grid_field_t
    use mesh_mod,               only : mesh_t

    class(stvec_t), allocatable, intent(out) :: stvec
    type(string_t),              intent(in)  :: mesh_names(:), field_names(:)
    integer(kind=4),             intent(in)  :: halo_width_xy, halo_width_z
    type(domain_t),              intent(in)  :: domain


    type(stvec_flexible_t), allocatable :: stvec_loc

    allocate(stvec_loc)

    call create_stvec_flexible_allocated(stvec_loc, field_names, mesh_names, &
                                         halo_width_xy, halo_width_z, domain)

    call move_alloc(stvec_loc, stvec)

end subroutine create_stvec_flexible

subroutine create_stvec_flexible_allocated(stvec, field_names, mesh_names, &
                                 halo_width_xy, halo_width_z, domain)

    use stvec_flexible_mod,     only : stvec_flexible_t
    use domain_mod,             only : domain_t
    use stvec_mod,              only : stvec_t
    use string_mod,             only : string_t
    use grid_field_factory_mod, only : create_grid_field
    use grid_field_mod,         only : grid_field_t
    use mesh_mod,               only : mesh_t

    class(stvec_flexible_t),     intent(inout) :: stvec
    type(string_t),              intent(in)  :: mesh_names(:), field_names(:)
    integer(kind=4),             intent(in)  :: halo_width_xy, halo_width_z
    type(domain_t),              intent(in)  :: domain

    type(grid_field_t),     allocatable :: tmp_field
    type(mesh_t),           pointer     :: tmp_mesh

    integer(kind=4) :: field_idx

    do field_idx = 1, ubound(field_names, 1)
        call domain%get_mesh(tmp_mesh, mesh_names(field_idx)%str)
        allocate(tmp_field)
        call create_grid_field(tmp_field, halo_width_xy, halo_width_z, tmp_mesh)
        call stvec%fields%add_grid_field(tmp_field, field_names(field_idx)%str)
        deallocate(tmp_field)
    end do

    allocate(stvec%field_names, source = field_names)
    allocate(stvec%mesh_names,  source = mesh_names )

end subroutine create_stvec_flexible_allocated

end module stvec_flexible_factory_mod
