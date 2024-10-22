module dp_interpolator_factory_mod

use abstract_dp_interpolator_mod, only : dp_interpolator_t
use dp_interpolators_mod,         only : nearest_neighbour_interpolator_t, &
                                         bilinear_interpolator_t,          &
                                         bicubic_Ah_interpolator_t,        &
                                         trilinear_Ah_interpolator_t,         &
                                         tricubic_Ah_interpolator_t
use exchange_factory_mod,         only : create_xy_points_halo_exchange,   &
                                         create_xyz_points_halo_exchange
use domain_mod,                   only : domain_t
use parcomm_mod,                  only : parcomm_global

implicit none

contains

subroutine create_dp_interpolator(interpolator,domain,interpolator_name)
    class(dp_interpolator_t), allocatable, intent(out) :: interpolator
    type(domain_t),   intent(in) :: domain
    character(len=*), intent(in) :: interpolator_name

    select case(interpolator_name)
    case("nearest_neighbour")
        interpolator = nearest_neighbour_interpolator_t()
        call interpolator%init(mesh_name="xy",domain=domain)
    case("bilinear")
        allocate(bilinear_interpolator_t :: interpolator)
        select type (interpolator)
        type is (bilinear_interpolator_t)
            interpolator%exchange_halo = create_xy_points_halo_exchange(domain%partition, domain%parcomm, &
                                                                        domain%topology,  1, 'full')
        end select
        call interpolator%init(mesh_name="xy", domain=domain)
    case ("bicubic_Ah")
        call create_bicubic_Ah_interp(interpolator,domain)
    case ("trilinear_Ah", "trilinear_Ah_z")
        call create_trilinear_Ah_interp(interpolator,interpolator_name,domain)
    case ("tricubic_Ah", "tricubic_Ah_z")
        call create_tricubic_Ah_interp(interpolator,interpolator_name,domain)
    case default
        call parcomm_global%abort(__FILE__//" unknown interpolator name "//interpolator_name)
    end select

end subroutine

subroutine create_bicubic_Ah_interp(interpolator,domain)
    class(dp_interpolator_t), allocatable, intent(out) :: interpolator
    type(domain_t), intent(in) :: domain

    type(bicubic_Ah_interpolator_t), allocatable :: bicubic

    allocate(bicubic)

    bicubic%exchange_halo = create_xy_points_halo_exchange(domain%partition, domain%parcomm, &
                                                domain%topology,  3, 'full')

    call bicubic%init(mesh_name="xy", domain=domain)

    call move_alloc(bicubic,interpolator)
end subroutine

subroutine create_trilinear_Ah_interp(interpolator,interpolator_name,domain)
    class(dp_interpolator_t), allocatable, intent(out) :: interpolator
    character(len=*), intent(in) :: interpolator_name
    type(domain_t),   intent(in) :: domain

    type(trilinear_Ah_interpolator_t), allocatable :: trilinear

    allocate(trilinear)

    select case (interpolator_name)
    case ("trilinear_Ah")
        trilinear%exchange_halo = create_xy_points_halo_exchange(domain%partition, domain%parcomm, &
                                                                domain%topology,  3, 'full')
        call trilinear%init(mesh_name="xy", domain=domain)
    case ("trilinear_Ah_z")
        trilinear%exchange_halo = create_xyz_points_halo_exchange(domain%partition, domain%parcomm, &
                                                                 domain%topology,  3, 'full')
        call trilinear%init(mesh_name="xyz", domain=domain)
    case default
        call parcomm_global%abort("create_trilinear_Ah_interp, unknown interp name: "//interpolator_name)
    end select

    call move_alloc(trilinear,interpolator)
end subroutine

subroutine create_tricubic_Ah_interp(interpolator,interpolator_name,domain)
    class(dp_interpolator_t), allocatable, intent(out) :: interpolator
    character(len=*), intent(in) :: interpolator_name
    type(domain_t),   intent(in) :: domain

    type(tricubic_Ah_interpolator_t), allocatable :: tricubic

    allocate(tricubic)

    select case (interpolator_name)
    case ("tricubic_Ah")
        tricubic%exchange_halo = create_xy_points_halo_exchange(domain%partition, domain%parcomm, &
                                                                domain%topology,  3, 'full')
        call tricubic%init(mesh_name="xy", domain=domain)
    case ("tricubic_Ah_z")
        tricubic%exchange_halo = create_xyz_points_halo_exchange(domain%partition, domain%parcomm, &
                                                                 domain%topology,  3, 'full')
        call tricubic%init(mesh_name="xyz", domain=domain)
    case default
        call parcomm_global%abort("create_tricubic_Ah_interp, unknown interp name: "//interpolator_name)
    end select

    call move_alloc(tricubic,interpolator)
end subroutine

end module