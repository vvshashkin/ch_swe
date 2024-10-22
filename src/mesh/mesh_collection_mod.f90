module mesh_collection_mod

use linked_list_mod, only : linked_list_t
use mesh_mod,        only : mesh_t

implicit none

type, public :: mesh_collection_t
    private
    type(linked_list_t) :: mesh_list
contains
    procedure, public :: add_mesh
    procedure, public :: add_mesh_reference
    procedure, public :: get_mesh
end type mesh_collection_t

type, public :: mesh_reference_t
    type(mesh_t), pointer :: mesh_pointer
contains
end type mesh_reference_t

contains

subroutine add_mesh(this, new_mesh, mesh_name)

    class(mesh_collection_t), intent(inout) :: this
    type(mesh_t),             intent(in)    :: new_mesh
    character(len=*),         intent(in)    :: mesh_name

    call this%mesh_list%append(new_mesh, mesh_name)

end subroutine add_mesh

subroutine add_mesh_reference(this, mesh, mesh_name)

    class(mesh_collection_t), intent(inout) :: this
    type(mesh_t),    pointer, intent(in)    :: mesh
    character(len=*),         intent(in)    :: mesh_name
    type(mesh_reference_t), pointer :: mesh_reference

    allocate(mesh_reference)
    mesh_reference%mesh_pointer => mesh
    call this%mesh_list%append(mesh_reference, mesh_name)

end subroutine add_mesh_reference

subroutine get_mesh(this, mesh, mesh_name)

    class(mesh_collection_t),          intent(in)    :: this
    type(mesh_t),             pointer, intent(inout) :: mesh
    character(len=*),                  intent(in)    :: mesh_name

    class(*), pointer :: data

    data => null()
    mesh => null()

    call this%mesh_list%get_item(data, mesh_name)

    select type (data)
    type is (mesh_t)
        mesh => data
    type is (mesh_reference_t)
        mesh => data%mesh_pointer
    end select

end subroutine get_mesh

end module mesh_collection_mod
