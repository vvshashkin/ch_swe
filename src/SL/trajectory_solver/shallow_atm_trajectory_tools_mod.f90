module shallow_atm_trajectory_tools_mod

use grid_field_mod,         only : grid_field_t
use mesh_mod,               only : mesh_t
use metric_mod,             only : metric_t
use shallow_atm_metric_mod, only : shallow_atm_metric_t
use parcomm_mod,            only : parcomm_global

implicit none

interface calc_shallow_atm_xyz_wind
    module procedure :: calc_shallow_atm_xyz_wind_at_mesh_points
    module procedure :: calc_shallow_atm_xyz_wind_at_departure_points
end interface

contains

subroutine calc_shallow_atm_xyz_wind_at_mesh_points(vx, vy, vz, u, v, mesh)
    type(grid_field_t),   intent(inout) :: vx, vy, vz
    type(grid_field_t),   intent(in)    :: u, v
    type(mesh_t),         intent(in)    :: mesh

    integer(kind=4) :: i, j, k, t
    real(kind=8)    :: a1(3), a2(3)

    do t = mesh%ts, mesh%te
        do k = mesh%tile(t)%ks, mesh%tile(t)%ke
            do j = mesh%tile(t)%js, mesh%tile(t)%je
                do i = mesh%tile(t)%is, mesh%tile(t)%ie

                    a1(1:3) = mesh%tile(t)%a1(1:3,i,j,k)
                    a2(1:3) = mesh%tile(t)%a2(1:3,i,j,k)

                    vx%tile(t)%p(i,j,k) = u%tile(t)%p(i,j,k)*a1(1)+v%tile(t)%p(i,j,k)*a2(1)
                    vy%tile(t)%p(i,j,k) = u%tile(t)%p(i,j,k)*a1(2)+v%tile(t)%p(i,j,k)*a2(2)
                    vz%tile(t)%p(i,j,k) = u%tile(t)%p(i,j,k)*a1(3)+v%tile(t)%p(i,j,k)*a2(3)

                end do
            end do
        end do
    end do

end subroutine

subroutine calc_shallow_atm_xyz_wind_at_departure_points(vx, vy, vz, u, v, &
                                                         panel_ind, alpha, beta, mesh, metric)
    type(grid_field_t),   intent(inout) :: vx, vy, vz
    type(grid_field_t),   intent(in)    :: u, v, alpha, beta, panel_ind
    type(mesh_t),         intent(in)    :: mesh
    class(metric_t),      intent(in)    :: metric

    integer(kind=4) :: i, j, k, t
    real(kind=8)    :: a1(4), a2(4), r(3)

    select type (metric)
    class is (shallow_atm_metric_t)
        call metric%calc_cartesian_hor_wind(vx, vy, vz, u, v, &
                                            panel_ind, alpha, beta, &
                                            mesh, "contravariant")
    class default
        call parcomm_global%abort("calc_shallow_atm_xyz_wind_at_departure_points: metric is not shallow atmosphere metric")
    end select

end subroutine

subroutine calc_hor_native_wind(u, v, vx, vy, vz, mesh)
    type(grid_field_t),   intent(inout) :: u, v
    type(grid_field_t),   intent(in)    :: vx, vy, vz
    type(mesh_t),         intent(in)    :: mesh
    
    integer(kind=4) :: i, j, k, t
    real(kind=8)    :: b1(3), b2(3)

    do t = mesh%ts, mesh%te
        do k = mesh%tile(t)%ks, mesh%tile(t)%ke
            do j = mesh%tile(t)%js, mesh%tile(t)%je
                do i = mesh%tile(t)%is, mesh%tile(t)%ie

                    b1(1:3) = mesh%tile(t)%b1(1:3,i,j,k)
                    b2(1:3) = mesh%tile(t)%b2(1:3,i,j,k)

                    u%tile(t)%p(i,j,k) = vx%tile(t)%p(i,j,k)*b1(1)+vy%tile(t)%p(i,j,k)*b1(2)+vz%tile(t)%p(i,j,k)*b1(3)
                    v%tile(t)%p(i,j,k) = vx%tile(t)%p(i,j,k)*b2(1)+vy%tile(t)%p(i,j,k)*b2(2)+vz%tile(t)%p(i,j,k)*b2(3)

                end do
            end do
        end do
    end do
    
end subroutine

subroutine displace_by_wind_shallow_atm(x, y, z, vx, vy, vz, dt, scale, mesh)
    type(grid_field_t), intent(inout) :: x, y, z
    type(grid_field_t), intent(in)    :: vx, vy, vz
    real(kind=8),       intent(in)    :: dt, scale
    type(mesh_t),       intent(in)    :: mesh

    integer(kind=4) :: i, j, k, t
    real(kind=8)    :: d(3),psi, d_abs

    do t = mesh%ts, mesh%te
        do k = mesh%tile(t)%ks, mesh%tile(t)%ke
            do j = mesh%tile(t)%js, mesh%tile(t)%je
                do i = mesh%tile(t)%is, mesh%tile(t)%ie

                    d(1:3) = [vx%tile(t)%p(i,j,k), vy%tile(t)%p(i,j,k), &
                              vz%tile(t)%p(i,j,k)] * (dt/scale)
                    psi = max(sqrt(d(1)**2 + d(2)**2 + d(3)**2),1e-16_8)
                    d(1:3) = d(1:3)/psi*tan(psi)

                    d(1) = mesh%tile(t)%rx(i,j,k) - d(1)
                    d(2) = mesh%tile(t)%ry(i,j,k) - d(2)
                    d(3) = mesh%tile(t)%rz(i,j,k) - d(3)

                    d_abs = sqrt(d(1)**2 + d(2)**2 + d(3)**2)

                    x%tile(t)%p(i,j,k) = d(1) / d_abs
                    y%tile(t)%p(i,j,k) = d(2) / d_abs
                    z%tile(t)%p(i,j,k) = d(3) / d_abs

                end do
            end do
        end do
    end do

end subroutine

subroutine rotate_shallow_atm_hor_wind(vx, vy, vz, x, y, z, mesh)

    type(grid_field_t), intent(inout) :: vx, vy, vz
    type(grid_field_t), intent(in)    :: x, y, z
    type(mesh_t),       intent(in)    :: mesh

    integer(kind=4) :: i, j, k, t
    real(kind=8) :: xa, ya, za, xd, yd, zd
    real(kind=8) :: tdx, tdy, tdz, tax, tay, taz
    real(kind=8) :: dot_prod, c, tau2, vt

    do t = mesh%ts, mesh%te
        do k = mesh%tile(t)%ks, mesh%tile(t)%ke
            do j = mesh%tile(t)%js, mesh%tile(t)%je
                do i = mesh%tile(t)%is, mesh%tile(t)%ie

                    xa = mesh%tile(t)%rx(i,j,k)
                    ya = mesh%tile(t)%ry(i,j,k)
                    za = mesh%tile(t)%rz(i,j,k)

                    xd = x%tile(t)%p(i,j,k)
                    yd = y%tile(t)%p(i,j,k)
                    zd = z%tile(t)%p(i,j,k)

                    dot_prod = xa*xd+ya*yd+za*zd
                    c = 1.0_8 - dot_prod

                    tdx = xa - (1.0_8-c)*xd
                    tdy = ya - (1.0_8-c)*yd
                    tdz = za - (1.0_8-c)*zd

                    tax = (1.0_8-c)*xa - xd
                    tay = (1.0_8-c)*ya - yd
                    taz = (1.0_8-c)*za - zd

                    tau2 = max((1.0_8-c)**2-2.0_8*(1.0_8-c)*dot_prod+1.0_8,1e-16_8)

                    vt = (vx%tile(t)%p(i,j,k)*tdx+vy%tile(t)%p(i,j,k)*tdy+&
                          vz%tile(t)%p(i,j,k)*tdz)/tau2

                    vx%tile(t)%p(i,j,k) = vx%tile(t)%p(i,j,k) + vt*(tax-tdx)
                    vy%tile(t)%p(i,j,k) = vy%tile(t)%p(i,j,k) + vt*(tay-tdy)
                    vz%tile(t)%p(i,j,k) = vz%tile(t)%p(i,j,k) + vt*(taz-tdz)
                end do
            end do
        end do
    end do

end subroutine

subroutine displace_by_vertical_wind_eta(eta,eta_dot,dt,mesh)
    type(grid_field_t), intent(inout) :: eta
    type(grid_field_t), intent(in)    :: eta_dot
    real(kind=8),       intent(in)    :: dt
    type(mesh_t),       intent(in)    :: mesh

    integer(kind=4) :: i, j, k, t
    real(kind=8) :: eta_arr, eta_dp

    do t = mesh%ts, mesh%te
        do k = mesh%tile(t)%ks, mesh%tile(t)%ke
            eta_arr = mesh%tile(t)%get_eta(k)
            do j = mesh%tile(t)%js, mesh%tile(t)%je
                do i = mesh%tile(t)%is, mesh%tile(t)%ie
                    eta_dp = eta_arr - dt*eta_dot%tile(t)%p(i,j,k) / mesh%vertical_scale
                    eta%tile(t)%p(i,j,k) = max(0.0_8,min(1.0_8,eta_dp))
                end do
            end do
        end do
    end do

end subroutine

end module