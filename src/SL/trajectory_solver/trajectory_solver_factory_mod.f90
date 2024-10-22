module trajectory_solver_factory_mod

use domain_mod,                           only : domain_t
use abstract_trajectory_solver_mod,       only : trajectory_solver_t
use generic_config_mod,                   only : generic_config_t
use geosci_config_mod,                    only : geosci_config_t
use grid_field_factory_mod,               only : create_grid_field
use mesh_mod,                             only : mesh_t
use dep_points_interp_driver_factory_mod, only : create_dep_points_interp_driver
use stvec_flexible_factory_mod,           only : create_stvec_flexible_allocated
use string_mod,                           only : strings, string_t
use parcomm_mod,                          only : parcomm_global
use shallow_spherical_trajectory_solver_mod, only : shallow_spherical_trajectory_solver_t

implicit none

contains

subroutine create_trajectory_solver(trajectory_solver, config, domain)
    class(trajectory_solver_t), allocatable, intent(out) :: trajectory_solver
    class(generic_config_t), intent(inout) :: config
    type(domain_t),          intent(in)    :: domain

    character(len=:), allocatable :: trajectory_solver_name

    call config%get(trajectory_solver_name,"trajectory_solver_name")

    select case(trajectory_solver_name)
    case("shallow_atmosphere_trajectory_solver")
        call create_shallow_atm_trajectory_solver(trajectory_solver, config, domain)
    case default
        call parcomm_global%abort("create_trajectory_solver, unknown solver name "//trajectory_solver_name)
    end select

end subroutine

subroutine create_shallow_atm_trajectory_solver(trajectory_solver, config, domain)
    class(trajectory_solver_t), allocatable, intent(out) :: trajectory_solver
    class(generic_config_t), intent(inout) :: config
    type(domain_t),          intent(in)    :: domain

    type(shallow_spherical_trajectory_solver_t), allocatable :: shallow_spherical
    type(mesh_t), pointer :: mesh
    type(geosci_config_t) :: interp_driver_config, interp_config_u, interp_config_v, interp_config_eta_dot
    character(len=:), allocatable :: arrival_mesh, hor_wind_interp_name, &
                                     ver_wind_interp_name, dp_interp_driver_name
    integer(kind=4) :: max_cfl
    logical :: is_conf_correct
    type(string_t), allocatable :: unexpected_config_vars(:)
    type(string_t), allocatable :: wind_fields(:), wind_meshes(:)

    allocate(shallow_spherical)

    is_conf_correct = config%check_no_unexpected(&
                              strings("trajectory_solver_name",            &
                                      "dep_points_interp_driver_name",     &
                                      "hor_wind_interp_name",              &
                                      "ver_wind_interp_name",              &
                                      "max_cfl", "num_iter", "mode_2d",    &
                                      "arrival_mesh", "use_first_guess" ), &
                                       unexpected_config_vars)

    if(.not.is_conf_correct) call parcomm_global%abort("create trajectory solver, unexpected varnames in config: " // unexpected_config_vars(1)%str)

    call config%get(shallow_spherical%num_iter,"num_iter",default=3)
    call config%get(shallow_spherical%mode_2d,"mode_2d",default=.false.)
    call config%get(shallow_spherical%use_first_guess,"use_first_guess",default=.false.)
    call config%get(dp_interp_driver_name,"dep_points_interp_driver_name",default="WhiteDongarra")
    call config%get(hor_wind_interp_name,"hor_wind_interp_name",default="bicubic_Ah")
    call config%get(ver_wind_interp_name,"ver_wind_interp_name",default="bicubic_Ah")
    call config%get(arrival_mesh,"arrival_mesh")
    call config%get(max_cfl, "max_cfl", default = 9)

    call domain%get_mesh(shallow_spherical%arrival_mesh,arrival_mesh)
    call domain%get_mesh(mesh,arrival_mesh)

    call create_grid_field(shallow_spherical%vx_arr,0,0,mesh)
    call create_grid_field(shallow_spherical%vy_arr,0,0,mesh)
    call create_grid_field(shallow_spherical%vz_arr,0,0,mesh)
    call create_grid_field(shallow_spherical%vx    ,0,0,mesh)
    call create_grid_field(shallow_spherical%vy    ,0,0,mesh)
    call create_grid_field(shallow_spherical%vz    ,0,0,mesh)

    if(shallow_spherical%mode_2d) then
        wind_fields = strings("u","v")
        wind_meshes = strings(arrival_mesh,arrival_mesh)
    else
        wind_fields = strings("u","v","eta_dot")
        wind_meshes = strings(arrival_mesh,arrival_mesh,arrival_mesh)
    end if

    call create_stvec_flexible_allocated(shallow_spherical%wind_dp, &
                                         wind_fields, wind_meshes,0, 0, domain)

    call interp_driver_config%set("dep_points_interp_driver_name",dp_interp_driver_name)
    call interp_driver_config%set("max_cfl",max_cfl)
    call interp_driver_config%set("arrival_mesh",arrival_mesh)

    call interp_config_u%set("field_name","u")
    call interp_config_u%set("interpolation_name",hor_wind_interp_name)
    call interp_config_v%set("field_name","v")
    call interp_config_v%set("interpolation_name",hor_wind_interp_name)
    if( .not. shallow_spherical%mode_2d) then
        call interp_config_eta_dot%set("field_name","eta_dot")
        call interp_config_eta_dot%set("interpolation_name",ver_wind_interp_name)
        call interp_driver_config%set("interp_configs",[interp_config_u,interp_config_v,interp_config_eta_dot])
    else
        call interp_driver_config%set("interp_configs",[interp_config_u,interp_config_v])
    end if
 
    call create_dep_points_interp_driver(shallow_spherical%dp_interp_driver, &
                                         domain,interp_driver_config)

    call move_alloc(shallow_spherical, trajectory_solver)
    
end subroutine create_shallow_atm_trajectory_solver

end module