module orography_factory_mod

use orography_mod,              only : orography_collection_t, orography_t
use domain_mod,                 only : domain_t
use mesh_mod,                   only : mesh_t
use generic_config_mod,         only : generic_config_t
use grid_field_factory_mod,     only : create_grid_field
use parcomm_mod,                only : parcomm_global
use orography_test_field_mod,   only : orography_test_field_t, orography_test_grad_t
use Schar_orography_field_mod,  only : Schar_orography_grad_t, Schar_orography_field_t
use MIRW3d_orography_field_mod, only : MIRW3d_orography_field_t, MIRW3d_orography_grad_t
use baroclinic_instability_test_orography_mod, &
                                only : baroclinic_instability_test_orography_t, &
                                       baroclinic_instability_test_orography_grad_t
use test_fields_3d_mod,         only : scalar_field3d_t, vector_field3d_t
use parcomm_mod,                only : parcomm_global

implicit none

character(len=2), parameter :: orography_mesh_names(4) = ["o ","x ","y ","xy"]

contains

subroutine create_orography(orography, config, domain, halo_width)
    type(orography_collection_t), intent(out)   :: orography
    class(generic_config_t),      intent(inout) :: config
    type(domain_t),               intent(in)    :: domain
    integer(kind=4),              intent(in)    :: halo_width

    type(orography_t), pointer :: tmp_orog
    type(mesh_t),      pointer :: tmp_mesh
    real(kind=8) :: max_height
    class(scalar_field3d_t), allocatable :: orog_gen
    class(vector_field3d_t), allocatable :: orog_grad_gen
    integer(kind=4) :: iname
    character(len=:), allocatable :: orography_name

    do iname = 1, size(orography_mesh_names, 1)
        call domain%get_mesh(tmp_mesh, trim(orography_mesh_names(iname)))
        allocate(tmp_orog)
        call create_grid_field(tmp_orog%h,        halo_width, 0, tmp_mesh)
        call create_grid_field(tmp_orog%dh_alpha, halo_width, 0, tmp_mesh)
        call create_grid_field(tmp_orog%dh_beta,  halo_width, 0, tmp_mesh)
        call orography%add_orography(tmp_orog, trim(orography_mesh_names(iname)))
        deallocate(tmp_orog)
    end do

    call config%get(orography_name, "orography_name", default="zero_orography")
    call config%get(max_height,"max_height",default = 0.0_8)

    if(orography_name == "zero_orography") then
        do iname = 1, size(orography_mesh_names,1)
            call domain%get_mesh(tmp_mesh,trim(orography_mesh_names(iname)))
            call orography%get_orography(tmp_orog,trim(orography_mesh_names(iname)))
            call tmp_orog%h       %assign(0.0_8,tmp_mesh,halo_width)
            call tmp_orog%dh_alpha%assign(0.0_8,tmp_mesh,halo_width)
            call tmp_orog%dh_beta %assign(0.0_8,tmp_mesh,halo_width)
        end do

    else
        select case(orography_name)
        case("test_orography")
            orog_gen      = orography_test_field_t(max_height)
            orog_grad_gen = orography_test_grad_t(max_height, domain%mesh_o%scale)
        case("Schar_orography")
            orog_gen      = Schar_orography_field_t(max_height, domain%mesh_o%scale)
            orog_grad_gen = Schar_orography_grad_t(max_height, domain%mesh_o%scale)
        case("MIRW3d_orography")
            orog_gen      = MIRW3d_orography_field_t(max_height, domain%mesh_o%scale)
            orog_grad_gen = MIRW3d_orography_grad_t(max_height, domain%mesh_o%scale)
        case("Baroclinic_instability_test_orography")
            orog_gen      = baroclinic_instability_test_orography_t()
            orog_grad_gen = baroclinic_instability_test_orography_grad_t()
        case default
            call parcomm_global%abort("orography factory mod error - unknown orography name:"//&
                                      orography_name)
        end select

        do iname = 1, size(orography_mesh_names, 1)
            call domain%get_mesh(tmp_mesh, trim(orography_mesh_names(iname)))
            call orography%get_orography(tmp_orog, trim(orography_mesh_names(iname)))
            call orog_gen%get_scalar_field(tmp_orog%h, tmp_mesh, halo_width)
            call orog_grad_gen%get_x_component(tmp_orog%dh_alpha, tmp_mesh, halo_width,'covariant')
            call orog_grad_gen%get_y_component(tmp_orog%dh_beta , tmp_mesh, halo_width,'covariant')
        end do
    end if

end subroutine create_orography

end module orography_factory_mod
