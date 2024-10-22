module stvec_flexible_mod

use grid_field_collection_mod, only : grid_field_collection_t
use stvec_mod,                 only : stvec_t
use grid_field_mod,            only : grid_field_t
use parcomm_mod,               only : parcomm_global
use domain_mod,                only : domain_t
use string_mod,                only : string_t
use mesh_mod,                  only : mesh_t
use basic_collection_mod,      only : basic_collection_t

implicit none

type, extends(stvec_t) :: stvec_flexible_t
    type(grid_field_collection_t) :: fields
    type(string_t), allocatable   :: field_names(:)
    type(string_t), allocatable   :: mesh_names(:)
    type(basic_collection_t)      :: meta_information
contains
    procedure, public :: get_field
    procedure, public :: get_field_and_mesh
    procedure, public :: get_field_names
    procedure, public :: create_similar
    procedure, public :: create_similar_allocated
    procedure, public :: assign_s1
    procedure, public :: assign_v1
    procedure, public :: assign_s1v1
    procedure, public :: assign_s1v1s2v2
    procedure, public :: update_s1v1s2v2
    procedure, public :: update_s1v1
    procedure, public :: is_similar
end type stvec_flexible_t

contains

subroutine get_field(this, field, field_name, abort_on_fail)

    class(stvec_flexible_t),          intent(in)    :: this
    type(grid_field_t),      pointer, intent(inout) :: field
    character(len=*),                 intent(in)    :: field_name
    logical,                optional, intent(in)    :: abort_on_fail

    integer(kind=4) :: field_idx
    logical :: fail

    field => null()

    do field_idx = 1, ubound(this%field_names, 1)
        if (field_name == this%field_names(field_idx)%str) then
            call this%fields%get_grid_field(field, field_name)
            exit
        end if
    end do

    fail = .false.
    if(present(abort_on_fail)) fail = abort_on_fail .and. .not. associated(field)

    if(fail) &
        call parcomm_global%abort("stvec flexible get_field, cannot find field "//field_name)

end subroutine get_field

subroutine get_field_and_mesh(this, field, mesh, field_name, domain)

    class(stvec_flexible_t),          intent(in)    :: this
    type(grid_field_t),      pointer, intent(inout) :: field
    type(mesh_t),            pointer, intent(inout) :: mesh
    character(len=*),                 intent(in)    :: field_name
    type(domain_t),                   intent(in)    :: domain

    integer(kind=4) :: field_idx

    field => null()
    mesh  => null()

    do field_idx = 1, ubound(this%field_names, 1)
        if (field_name == this%field_names(field_idx)%str) then
            call this%fields%get_grid_field(field, field_name)
            call domain%get_mesh(mesh, this%mesh_names(field_idx)%str)
            exit
        end if
    end do

end subroutine get_field_and_mesh

function get_field_names(this) result(field_names)
    class(stvec_flexible_t), intent(in) :: this
    type(string_t), allocatable :: field_names(:)

    field_names = this%fields%get_grid_field_names()
end function

subroutine create_similar(this, destination)
    class(stvec_flexible_t),     intent(in)    :: this
    class(stvec_t), allocatable, intent(inout) :: destination

    type(stvec_flexible_t), allocatable :: dest_loc
    type(grid_field_t), pointer :: source_field
    type(grid_field_t), allocatable :: target_field
    integer(kind=4) :: field_idx

    allocate(dest_loc)

    call this%create_similar_allocated(dest_loc)

    call move_alloc(dest_loc, destination)

end subroutine create_similar

subroutine create_similar_allocated(this, destination)
    class(stvec_flexible_t),     intent(in)    :: this
    class(stvec_flexible_t),     intent(inout) :: destination

    type(grid_field_t), pointer :: source_field
    type(grid_field_t), allocatable :: target_field
    integer(kind=4) :: field_idx

    destination%field_names = this%field_names
    destination%mesh_names  = this%mesh_names

    do field_idx = 1, ubound(this%field_names, 1)
        call this%fields%get_grid_field(source_field, this%field_names(field_idx)%str)
        target_field = source_field%create_similar()
        call destination%fields%add_grid_field(target_field, this%field_names(field_idx)%str)
    end do

    call this%meta_information%copy_to(destination%meta_information)

end subroutine create_similar_allocated

subroutine assign_v1(this, v1, domain)

    class(stvec_flexible_t), intent(inout) :: this
    class(stvec_t),          intent(in)    :: v1
    class(domain_t),         intent(in)    :: domain

    type(mesh_t), pointer :: mesh
    type(grid_field_t), pointer :: gf1, out
    integer(kind=4) :: field_idx
    logical :: is_stvecs_compatible

    select type (v1)
    class is (stvec_flexible_t)
        is_stvecs_compatible = this%is_similar(v1)
        if ( .not. is_stvecs_compatible ) then
            call parcomm_global%abort("stvec_flexible_t%assign_s1v1: uncompatible inputs")
        end if

        do field_idx = 1, ubound(this%field_names, 1)
            call domain%get_mesh(mesh, this%mesh_names(field_idx)%str)

            call   v1%fields%get_grid_field(gf1, v1%field_names(field_idx)%str)
            call this%fields%get_grid_field(out, this%field_names(field_idx)%str)

            call out%assign(gf1, mesh)
        end do
    class default
        call parcomm_global%abort("stvec_flexible_t%assign_s1v1: types error")
    end select

end subroutine assign_v1

subroutine assign_s1v1(this, scalar1, v1, domain)

    class(stvec_flexible_t), intent(inout) :: this
    real(kind=8),            intent(in)    :: scalar1
    class(stvec_t),          intent(in)    :: v1
    class(domain_t),         intent(in)    :: domain

    type(mesh_t), pointer :: mesh
    type(grid_field_t), pointer :: gf1, out
    integer(kind=4) :: field_idx
    logical :: is_stvecs_compatible

    select type (v1)
    class is (stvec_flexible_t)
        is_stvecs_compatible = this%is_similar(v1)
        if ( .not. is_stvecs_compatible ) then
            call parcomm_global%abort("stvec_flexible_t%assign_s1v1: uncompatible inputs")
        end if

        do field_idx = 1, ubound(this%field_names, 1)
            call domain%get_mesh(mesh, this%mesh_names(field_idx)%str)

            call   v1%fields%get_grid_field(gf1, v1%field_names(field_idx)%str)
            call this%fields%get_grid_field(out, this%field_names(field_idx)%str)

            call out%assign(scalar1, gf1, mesh)
        end do
    class default
        call parcomm_global%abort("stvec_flexible_t%assign_s1v1: types error")
    end select

end subroutine assign_s1v1

subroutine assign_s1v1s2v2(this, scalar1, v1, scalar2, v2, domain)

    class(stvec_flexible_t), intent(inout) :: this
    real(kind=8),            intent(in)    :: scalar1, scalar2
    class(stvec_t),          intent(in)    :: v1, v2
    class(domain_t),         intent(in)    :: domain

    type(mesh_t), pointer :: mesh
    type(grid_field_t), pointer :: gf1, gf2, out
    integer(kind=4) :: field_idx
    logical :: is_stvecs_compatible

    select type (v1)
    class is (stvec_flexible_t)
        select type(v2)
        class is (stvec_flexible_t)

            is_stvecs_compatible = this%is_similar(v1) .and. this%is_similar(v2)
            if ( .not. is_stvecs_compatible ) then
                call parcomm_global%abort("stvec_flexible_t%assign_s1v1s2v2: uncompatible inputs")
            end if

            do field_idx = 1, ubound(this%field_names, 1)
                call domain%get_mesh(mesh, this%mesh_names(field_idx)%str)

                call   v1%fields%get_grid_field(gf1, v1%field_names(field_idx)%str)
                call   v2%fields%get_grid_field(gf2, v2%field_names(field_idx)%str)
                call this%fields%get_grid_field(out, this%field_names(field_idx)%str)

                call out%assign(scalar1, gf1, scalar2, gf2, mesh)
            end do
    class default
        call parcomm_global%abort("stvec_flexible_t%assign_s1v1s2v2: types error")
    end select
    class default
        call parcomm_global%abort("stvec_flexible_t%assign_s1v1s2v2: types error")
    end select

end subroutine assign_s1v1s2v2
subroutine update_s1v1s2v2(this, scalar1, v1, scalar2, v2, domain)

    class(stvec_flexible_t), intent(inout) :: this
    real(kind=8),            intent(in)    :: scalar1, scalar2
    class(stvec_t),          intent(in)    :: v1, v2
    class(domain_t),         intent(in)    :: domain

    type(mesh_t), pointer :: mesh
    type(grid_field_t), pointer :: gf1, gf2, out
    integer(kind=4) :: field_idx
    logical :: is_stvecs_compatible

    select type (v1)
    class is (stvec_flexible_t)
        select type(v2)
        class is (stvec_flexible_t)

            is_stvecs_compatible = this%is_similar(v1) .and. this%is_similar(v2)
            if ( .not. is_stvecs_compatible ) then
                call parcomm_global%abort("stvec_flexible_t%update_s1v1s2v2: uncompatible inputs")
            end if

            do field_idx = 1, ubound(this%field_names, 1)
                call domain%get_mesh(mesh, this%mesh_names(field_idx)%str)

                call   v1%fields%get_grid_field(gf1,   v1%field_names(field_idx)%str)
                call   v2%fields%get_grid_field(gf2,   v2%field_names(field_idx)%str)
                call this%fields%get_grid_field(out, this%field_names(field_idx)%str)

                call out%update(scalar1, gf1, scalar2, gf2, mesh)
            end do
    class default
        call parcomm_global%abort("stvec_flexible_t%update_s1v1s2v2: types error")
    end select
    class default
        call parcomm_global%abort("stvec_flexible_t%update_s1v1s2v2: types error")
    end select

end subroutine update_s1v1s2v2
subroutine update_s1v1(this, scalar1, v1, domain)

    class(stvec_flexible_t), intent(inout) :: this
    real(kind=8),            intent(in)    :: scalar1
    class(stvec_t),          intent(in)    :: v1
    class(domain_t),         intent(in)    :: domain

    type(mesh_t), pointer :: mesh
    type(grid_field_t), pointer :: gf1, gf2, out
    integer(kind=4) :: field_idx
    logical :: is_stvecs_compatible

    select type (v1)
    class is (stvec_flexible_t)

        is_stvecs_compatible = this%is_similar(v1)
        if ( .not. is_stvecs_compatible ) then
            call parcomm_global%abort("stvec_flexible_t%update_s1v1: uncompatible inputs")
        end if

        do field_idx = 1, ubound(this%field_names, 1)
            call domain%get_mesh(mesh, this%mesh_names(field_idx)%str)

            call   v1%fields%get_grid_field(gf1,   v1%field_names(field_idx)%str)
            call this%fields%get_grid_field(out, this%field_names(field_idx)%str)

            call out%update(scalar1, gf1, mesh)
        end do
    class default
        call parcomm_global%abort("stvec_flexible_t%update_s1v1: types error")
    end select

end subroutine update_s1v1
subroutine assign_s1(this, scalar1, domain)

    class(stvec_flexible_t), intent(inout) :: this
    real(kind=8),       intent(in)       :: scalar1
    class(domain_t),    intent(in)       :: domain

    type(mesh_t),       pointer :: mesh
    type(grid_field_t), pointer :: out
    integer(kind=4) :: field_idx

    do field_idx = 1, ubound(this%field_names, 1)
        call domain%get_mesh(mesh, this%mesh_names(field_idx)%str)
        call this%fields%get_grid_field(out, this%field_names(field_idx)%str)
        call out%assign(scalar1, mesh)
    end do

end subroutine assign_s1

function is_similar(this, other) result(out)

    class(stvec_flexible_t), intent(in) :: this
    type (stvec_flexible_t), intent(in) :: other
    logical :: out

    integer(kind=4) :: field_idx

    out = .false.

    if ( ubound(this%field_names, 1) /= ubound(other%field_names, 1) ) return

    do field_idx = 1, ubound(this%field_names, 1)
        if (this%field_names(field_idx)%str /= other%field_names(field_idx)%str) return
        if (this%mesh_names(field_idx)%str /= other%mesh_names(field_idx)%str)   return
    end do

    out = .true.

end function is_similar

end module stvec_flexible_mod
