program get_latlon_main

use parcomm_mod,        only : init_global_parallel_enviroment, deinit_global_parallel_enviroment
use domain_mod,         only : domain_t
use domain_factory_mod, only : create_domain
use mesh_mod,           only : mesh_t
use sph_coords_mod,     only : cart2sph

integer(kind=4),  parameter :: N = 192
character(len=*), parameter :: mesh_name = "xy"

type(domain_t) :: domain
type(mesh_t), pointer :: mesh
integer(kind=4) :: i, j, t, ip
real(kind=8)    :: lon, lat

call init_global_parallel_enviroment()

call create_domain(domain, "cube", "Ah", N, 1)

call domain%get_mesh(mesh,mesh_name)

open(117,file="lats_192.dat",access="direct",recl=2)
open(118,file="lons_192.dat",access="direct",recl=2)

print *, mesh%tile(1)%nx, mesh%tile(1)%ny

ip = 1
do t = mesh%ts, mesh%te
    do j = mesh%tile(t)%js, mesh%tile(t)%je
        do i = mesh%tile(t)%is, mesh%tile(t)%ie

            call cart2sph(mesh%tile(t)%rx(i,j,1),mesh%tile(t)%ry(i,j,1),&
                          mesh%tile(t)%rz(i,j,1),lon,lat)

            write(117,rec=ip) lat
            write(118,rec=ip) lon
            ip = ip+1

        end do
    end do
end do

call deinit_global_parallel_enviroment()

end program