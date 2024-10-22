module domain_mod

use topology_mod,        only : topology_t
use metric_mod,          only : metric_t
use partition_mod,       only : partition_t
use mesh_mod,            only : mesh_t
use parcomm_mod,         only : parcomm_t
use orography_mod,       only : orography_collection_t
use mesh_collection_mod, only : mesh_collection_t
use generic_config_mod,  only : generic_config_t

implicit none

type, public :: domain_t

    class(generic_config_t), allocatable :: config

    class(topology_t), allocatable  :: topology
    class(metric_t),   allocatable  :: metric
    character(len=:),  allocatable  :: horizontal_staggering
    type(domain_t),    allocatable  :: domain_2d
    type(parcomm_t)                 :: parcomm
    type(partition_t)               :: partition
    type(orography_collection_t)    :: orography

    type(mesh_collection_t) :: mesh_depository ! physical storage of domain meshes
    type(mesh_collection_t) :: mesh_collection ! references to meshes from mesh_depository

    ! pointers to meshes from mesh_depository. To be removed
    type(mesh_t), pointer :: mesh_o, mesh_x, mesh_y, mesh_xy, mesh_z, mesh_xyz
    type(mesh_t), pointer :: mesh_u, mesh_v, mesh_p, mesh_q, mesh_w

contains
    procedure, public :: get_mesh
end type domain_t

contains

subroutine get_mesh(this, mesh, name)

    class(domain_t),  intent(in)             :: this
    type(mesh_t),     intent(inout), pointer :: mesh
    character(len=*), intent(in)             :: name

    ! mesh => null()
    call this%mesh_collection%get_mesh(mesh, name)

end subroutine get_mesh

end module domain_mod
