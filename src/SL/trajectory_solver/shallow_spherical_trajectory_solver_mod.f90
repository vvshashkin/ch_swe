module shallow_spherical_trajectory_solver_mod

use domain_mod,                            only : domain_t
use stvec_flexible_mod,                    only : stvec_flexible_t
use abstract_trajectory_solver_mod,        only : trajectory_solver_t
use mesh_mod,                              only : mesh_t
use grid_field_mod,                        only : grid_field_t
use shallow_atm_trajectory_tools_mod,      only : calc_shallow_atm_xyz_wind,    &
                                                  displace_by_wind_shallow_atm, &
                                                  rotate_shallow_atm_hor_wind,  &
                                                  displace_by_vertical_wind_eta
use stvec_flexible_mod,                    only : stvec_flexible_t
use abstract_dep_points_interp_driver_mod, only : dep_points_interp_driver_t

use timer_mod, only : start_timer, stop_timer

implicit none

type, extends(trajectory_solver_t) :: shallow_spherical_trajectory_solver_t
    integer(kind=4)    :: num_iter
    logical            :: mode_2d, use_first_guess
    type(mesh_t), pointer :: arrival_mesh => null()
    type(grid_field_t) :: vx_arr, vy_arr, vz_arr ! cartesian wind at arrival point
    type(grid_field_t) :: vx, vy, vz ! total wind (arrival+departure)
    class(dep_points_interp_driver_t), allocatable :: dp_interp_driver
    type(stvec_flexible_t) :: wind_dp !buffer for wind values interpolated to departure point guess
    contains
        procedure :: find_departure_points
end type

contains
subroutine find_departure_points(this, dep_points_coords, dt, wind_arr, wind_dep, domain)

    class(shallow_spherical_trajectory_solver_t), intent(inout) :: this

    class(stvec_flexible_t), intent(inout) :: dep_points_coords
    real(kind=8),            intent(in)    :: dt
    class(stvec_flexible_t), intent(inout) :: wind_arr, wind_dep
    type(domain_t),          intent(in)    :: domain

    type(mesh_t), pointer :: mesh
    type(grid_field_t), pointer :: x, y, z, alpha, beta, panel_ind, eta, u, v, eta_dot, eta_dot_dp
    real(kind=8) :: scale
    integer(kind=4) :: iter

    mesh => this%arrival_mesh

    call dep_points_coords%get_field(x,"x")
    call dep_points_coords%get_field(y,"y")
    call dep_points_coords%get_field(z,"z")
    call dep_points_coords%get_field(alpha,"alpha")
    call dep_points_coords%get_field(beta,"beta")
    call dep_points_coords%get_field(panel_ind,"panel_ind")
    if(.not. this%mode_2d) call dep_points_coords%get_field(eta,"eta")

    call wind_arr%get_field(u,"u")
    call wind_arr%get_field(v,"v")
    if(.not. this%mode_2d) call wind_arr%get_field(eta_dot,"eta_dot")

    scale = domain%metric%scale

    call calc_shallow_atm_xyz_wind(this%vx_arr, this%vy_arr, this%vz_arr, u, v, mesh)

    !if: use zero displacement as trajectory first-guess
    !else: coordinates in dep_point_coords are used as first guess
    if(.not. this%use_first_guess) then

        call displace_by_wind_shallow_atm(x, y, z, this%vx_arr, this%vy_arr, &
                                          this%vz_arr, dt, scale, mesh)

        if(.not. this%mode_2d) &
            call displace_by_vertical_wind_eta(eta,eta_dot,dt,mesh)

        call domain%metric%transform_xyz_to_native(alpha,beta,panel_ind, x, y, z, mesh)

    end if

    call this%wind_dp%get_field(u,"u")
    call this%wind_dp%get_field(v,"v")
    if(.not. this%mode_2d) call this%wind_dp%get_field(eta_dot_dp,"eta_dot")

    do iter = 1, this%num_iter

        call this%dp_interp_driver%do(this%wind_dp, wind_dep, dep_points_coords, domain)

        call calc_shallow_atm_xyz_wind(this%vx, this%vy, this%vz, u, v, &
                                       panel_ind, alpha, beta, mesh, domain%metric)
        call rotate_shallow_atm_hor_wind(this%vx, this%vy, this%vz, x, y, z, mesh)

        call this%vx%assign(0.5_8, this%vx, 0.5_8, this%vx_arr, mesh)
        call this%vy%assign(0.5_8, this%vy, 0.5_8, this%vy_arr, mesh)
        call this%vz%assign(0.5_8, this%vz, 0.5_8, this%vz_arr, mesh)

        call displace_by_wind_shallow_atm(x, y, z, this%vx, this%vy, &
                                          this%vz, dt, scale, mesh)

        call domain%metric%transform_xyz_to_native(alpha,beta,panel_ind, x, y, z, mesh)

        if(.not. this%mode_2d) then
            call eta_dot_dp%assign(0.5_8,eta_dot_dp,0.5_8,eta_dot,mesh)
            call displace_by_vertical_wind_eta(eta,eta_dot_dp,dt,mesh)
        end if
    end do

end subroutine find_departure_points

end module