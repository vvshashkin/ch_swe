module domain_factory_mod

use domain_mod,                only : domain_t
use topology_factory_mod,      only : create_topology
use cubed_sphere_topology_mod, only : cubed_sphere_topology_t
use metric_mod,                only : metric_t
use metric_factory_mod,        only : create_metric
use orography_mod,             only : orography_collection_t, orography_t
use orography_factory_mod,     only : create_orography
use generic_config_mod,        only : generic_config_t
use geosci_config_mod,         only : geosci_config_t
use mesh_factory_mod,          only : create_mesh
use parcomm_factory_mod,       only : create_parcomm
use parcomm_mod,               only : parcomm_global, parcomm_t
use string_mod,                only : string_t
use partition_mod,             only : partition_t
use mpi

implicit none

character(len=*), private, parameter :: VERTICAL_STAGGERING_DEFAULT = "None"
character(len=*), private, parameter :: METRIC_DEFAULT = "shallow_atmosphere_metric"
real(kind=8),     private, parameter :: HTOP_DEFAULT = 1.0_8
integer(kind=4),  private, parameter :: HALO_WIDTH_DEFAULT = 8

interface create_domain
    module procedure :: create_domain_by_config
    module procedure :: create_domain_by_arguments
end interface

contains

subroutine create_domain_by_arguments(domain, topology_type, staggering_type, nh, nz, &
                                      parcomm)

    type(domain_t),   intent(out) :: domain
    character(len=*), intent(in)  :: topology_type
    character(len=*), intent(in)  :: staggering_type
    integer(kind=4),  intent(in)  :: nh, nz
    type(parcomm_t),  optional, intent(in) :: parcomm

    type(geosci_config_t) :: config_domain

    call config_domain%set("N",nh)
    call config_domain%set("Nz",nz)
    call config_domain%set("horizontal_staggering",staggering_type)
    call config_domain%set("vertical_staggering", VERTICAL_STAGGERING_DEFAULT)
    call config_domain%set("metric%metric_type",METRIC_DEFAULT)
    call config_domain%set("topology_type",topology_type)

    call create_domain_by_config(domain,config_domain,parcomm)

end subroutine create_domain_by_arguments

subroutine create_domain_by_config(domain, config, parcomm, partition)

    use partition_factory_mod, only : extract_2d_partition

    type(domain_t),             intent(out)   :: domain
    class(generic_config_t),    intent(inout) :: config
    type(parcomm_t),  optional, intent(in)    :: parcomm
    type(partition_t), optional, intent(in)   :: partition

    type(partition_t) :: partition_2d

    integer(kind=4) :: halo_width
    integer(kind=4) :: mpi_comm_local
    integer(kind=4) :: Nz_inter

    type(domain_t)  :: domain_2d
    class(generic_config_t), allocatable :: config_metric, config_orography
    character(len=:), allocatable :: vertical_staggering, topology_type
    integer(kind=4) :: N, Nz
    logical :: found

    type(string_t), allocatable :: unexpected_vars(:), expected_vars(:)
    logical         :: check
    integer(kind=4) :: i

    allocate(domain%config, source = config)

    expected_vars = [string_t("metric"),&
                     string_t("orography"), &
                     string_t("N"), string_t("Nz"), &
                     string_t("horizontal_staggering"), &
                     string_t("vertical_staggering"), &
                     string_t("topology_type")]

    check = config%check_no_unexpected(expected_vars,&
                        unexpected_varnames=unexpected_vars)
    if(.not. check) then
        do i = 1, size(unexpected_vars,1)
            call parcomm_global%print("metric factory WARNING: config variable "//unexpected_vars(i)%str//" presumably is not expected")
        end do
    end if

    call config%get(N,"N")
    call config%get(Nz,"Nz")
    call config%get(halo_width,"halo_width",default=HALO_WIDTH_DEFAULT)

    call config%get(topology_type,"topology_type")
    domain%topology = create_topology(topology_type)

    call config%get(config_metric,"metric")
    call create_metric(domain%metric, domain%topology, config_metric)

    call config%get(domain%horizontal_staggering, "horizontal_staggering")

    mpi_comm_local = parcomm_global%comm_w
    if(present(parcomm)) then
        domain%parcomm = parcomm
    else
        !call create_parcomm(parcomm_glocal%comm_w, domain%parcomm)
        domain%parcomm = parcomm_global
    end if

    call config%get(vertical_staggering, "vertical_staggering", default=VERTICAL_STAGGERING_DEFAULT)
    if(vertical_staggering == "CharneyPhilips") then
        Nz_inter = Nz+1
    else
        Nz_inter = max(Nz-1,1)
    end if

    if (present(partition)) then
        domain%partition = partition
    else
        call domain%partition%init(N, Nz, Nz_inter, max(1,domain%parcomm%np/6), &
                                   domain%parcomm%myid, domain%parcomm%Np,             &
                                   domain%horizontal_staggering,  &
                                   strategy = 'default')
    end if

    !initialize orographic curvilinear domain if orography config is present
    call config%get(config_orography,"orography",found=found)
    if(found) then
        allocate(domain%domain_2d)
        call extract_2d_partition(partition_2d, domain%partition)
        call create_2d_subdomain(domain%domain_2d, config, partition_2d, domain%parcomm)
        call create_orography(domain%orography,config_orography, &
                              domain%domain_2d,halo_width)
    end if

    call create_domain_meshes(domain,domain%horizontal_staggering, &
                              vertical_staggering,halo_width, &
                              domain%orography)

end subroutine create_domain_by_config

subroutine create_2d_subdomain(domain, config, partition_2d, parcomm)

    use mesh_factory_mod,    only : create_mesh
    use parcomm_factory_mod, only : create_parcomm
    use parcomm_mod,         only : parcomm_global, parcomm_t

    type(domain_t),             intent(out)    :: domain
    class(generic_config_t),    intent(inout)  :: config
    type(partition_t),          intent(in)     :: partition_2d
    type(parcomm_t),  optional, intent(in)     :: parcomm

    integer(kind=4) :: N, Nz, halo_width
    character(len=:), allocatable :: topology_type, horizontal_staggering
    character(len=*), parameter :: vertical_staggering="None"
    integer(kind=4) :: mpi_comm_local
    class(generic_config_t), allocatable :: config_metric

    type(domain_t)  :: domain_2d

    !have to be passed as an argument in future
    call config%get(N,"N")
    Nz = 1
    call config%get(halo_width,"halo_width",default=HALO_WIDTH_DEFAULT)
    call config%get(topology_type,"topology_type")
    call config%get(horizontal_staggering,"horizontal_staggering")

    domain%topology = create_topology(topology_type)

    call config%get(config_metric,"metric")
    call create_metric(domain%metric, domain%topology, config_metric)

    call config%get(domain%horizontal_staggering, "horizontal_staggering")

    if(present(parcomm)) then
        domain%parcomm = parcomm
    else
        domain%parcomm = parcomm_global
    end if

    domain%partition = partition_2d

    !passing empty orography is also ok
    call create_domain_meshes(domain,domain%horizontal_staggering, &
                              vertical_staggering,halo_width, &
                              domain%orography)

end subroutine create_2d_subdomain

subroutine create_domain_meshes(domain,horizontal_staggering, &
                                vertical_staggering,halo_width,orography)

    use mesh_mod,            only : mesh_t
    use mesh_collection_mod, only : mesh_reference_t
    use tiles_mod,           only : tiles_t

    type(domain_t),               intent(inout) :: domain
    character(len=*),             intent(in)    :: horizontal_staggering,&
                                                   vertical_staggering
    integer(kind=4),              intent(in)    :: halo_width
    type(orography_collection_t), intent(in)    :: orography

    type(mesh_t),  allocatable :: tmp_mesh
    type(orography_t), pointer :: tmp_orog
    type(tiles_t) :: tmp_tiles

    type(mesh_t), pointer :: mesh_ptr
    real(kind=8)    :: shift_zc, shift_zi!, shift_xyz(3)
    integer(kind=4) :: iname
    character(len=3), allocatable :: mesh_names(:), orography_names(:)
    real(kind=8) :: mesh_shifts(3,8)

    if (vertical_staggering == "CharneyPhilips") then
        shift_zc = 0.5_8
        shift_zi = 0.0_8
        mesh_names = ["o  ", "x  ", "y  ", "xy ", "z  ", "xz ", "yz ", "xyz"]
    else
        shift_zc = 0.0_8
        shift_zi = 0.5_8
        mesh_names = ["o  ", "x  ", "y  ", "xy ", "z  ", "xz ", "yz ", "xyz"]
    end if
    orography_names = ["o  ", "x  ", "y  ", "xy ", "o  ", "x  ","y  ", "xy "]
    mesh_shifts = reshape([0.5_8, 0.5_8, shift_zc, & !mesh_o
                           0.0_8, 0.5_8, shift_zc, & !mesh_x
                           0.5_8, 0.0_8, shift_zc, & !mesh_y
                           0.0_8, 0.0_8, shift_zc, & !mesh_xy
                           0.5_8, 0.5_8, shift_zi, & !mesh_z
                           0.0_8, 0.5_8, shift_zi, & !mesh_xz
                           0.5_8, 0.0_8, shift_zi, & !mesh_yz
                           0.0_8, 0.0_8, shift_zi],& !mesh_xyz
                           [3,8])

    do iname = 1, size(mesh_names)
        allocate(tmp_mesh)
        call orography%get_orography(tmp_orog, trim(orography_names(iname)))
        call domain%partition%get_tiles(trim(mesh_names(iname)), tmp_tiles)
        call create_mesh(tmp_mesh,  domain%partition, domain%metric, halo_width, &
                         tmp_tiles, mesh_shifts(1:3,iname), tmp_orog)
        call domain%mesh_depository%add_mesh(tmp_mesh, trim(mesh_names(iname)))
        deallocate(tmp_mesh)
        call domain%mesh_depository%get_mesh(mesh_ptr, trim(mesh_names(iname)))
        call domain%mesh_collection%add_mesh_reference(mesh_ptr, trim(mesh_names(iname)))
    end do

    call domain%mesh_depository%get_mesh(domain%mesh_o,   'o')
    call domain%mesh_depository%get_mesh(domain%mesh_x,   'x')
    call domain%mesh_depository%get_mesh(domain%mesh_y,   'y')
    call domain%mesh_depository%get_mesh(domain%mesh_xy,  'xy')
    call domain%mesh_depository%get_mesh(domain%mesh_z,   'z')
    call domain%mesh_depository%get_mesh(domain%mesh_xyz, 'xyz')

    select case(horizontal_staggering)
    case ('A')
        call domain%mesh_depository%get_mesh(mesh_ptr, 'o')
        call domain%mesh_collection%add_mesh_reference(mesh_ptr, 'p')
        call domain%mesh_collection%add_mesh_reference(mesh_ptr, 'u')
        call domain%mesh_collection%add_mesh_reference(mesh_ptr, 'v')
        call domain%mesh_collection%add_mesh_reference(mesh_ptr, 'q')
        domain%mesh_p => domain%mesh_o
        domain%mesh_u => domain%mesh_o
        domain%mesh_v => domain%mesh_o
        domain%mesh_q => domain%mesh_o

    case ('Ah') !all degrees of freedom at corner points
        call domain%mesh_depository%get_mesh(mesh_ptr, 'xy')
        call domain%mesh_collection%add_mesh_reference(mesh_ptr, 'p')
        call domain%mesh_collection%add_mesh_reference(mesh_ptr, 'u')
        call domain%mesh_collection%add_mesh_reference(mesh_ptr, 'v')
        call domain%mesh_collection%add_mesh_reference(mesh_ptr, 'q')
        domain%mesh_p => domain%mesh_xy
        domain%mesh_u => domain%mesh_xy
        domain%mesh_v => domain%mesh_xy
        domain%mesh_q => domain%mesh_xy
    case ('C')
        call domain%mesh_depository%get_mesh(mesh_ptr, 'o')
        call domain%mesh_collection%add_mesh_reference(mesh_ptr, 'p')

        call domain%mesh_depository%get_mesh(mesh_ptr, 'x')
        call domain%mesh_collection%add_mesh_reference(mesh_ptr, 'u')

        call domain%mesh_depository%get_mesh(mesh_ptr, 'y')
        call domain%mesh_collection%add_mesh_reference(mesh_ptr, 'v')

        call domain%mesh_depository%get_mesh(mesh_ptr, 'xy')
        call domain%mesh_collection%add_mesh_reference(mesh_ptr, 'q')
        domain%mesh_p => domain%mesh_o
        domain%mesh_u => domain%mesh_x
        domain%mesh_v => domain%mesh_y
        domain%mesh_q => domain%mesh_xy
    case ('Ch')
        call domain%mesh_depository%get_mesh(mesh_ptr, 'xy')
        call domain%mesh_collection%add_mesh_reference(mesh_ptr, 'p')

        call domain%mesh_depository%get_mesh(mesh_ptr, 'y')
        call domain%mesh_collection%add_mesh_reference(mesh_ptr, 'u')

        call domain%mesh_depository%get_mesh(mesh_ptr, 'x')
        call domain%mesh_collection%add_mesh_reference(mesh_ptr, 'v')

        call domain%mesh_depository%get_mesh(mesh_ptr, 'o')
        call domain%mesh_collection%add_mesh_reference(mesh_ptr, 'q')

        domain%mesh_p => domain%mesh_xy
        domain%mesh_u => domain%mesh_y
        domain%mesh_v => domain%mesh_x
        domain%mesh_q => domain%mesh_o
    case default
        call parcomm_global%abort("domain_factory_mod, unknown horizontal staggering: "//horizontal_staggering)
    end select
    select case(vertical_staggering)
    case("None")
        domain%mesh_w => domain%mesh_p
    case("CharneyPhilips")
            select case(horizontal_staggering)
            case ('A','C')
                domain%mesh_w => domain%mesh_z
            case ('Ah','Ch')
                domain%mesh_w => domain%mesh_xyz
            case default
                call parcomm_global%abort("domain_factory_mod, unknown horizontal staggering: "//horizontal_staggering)
            end select
    case default
        call parcomm_global%abort("domain_factory_mod, unknown vertical staggering type: "//vertical_staggering)
    end select
    mesh_ptr => domain%mesh_w
    call domain%mesh_collection%add_mesh_reference(mesh_ptr, 'w')

end subroutine create_domain_meshes

end module domain_factory_mod
