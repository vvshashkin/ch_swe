module ecs_halo_vec_Ch_factory_mod

use domain_mod,                  only : domain_t
use halo_mod,                    only : halo_vec_t
use ecs_halo_vec_Ch_mod,         only : ecs_halo_vec_Ch_t
use exchange_abstract_mod,       only : exchange_t
use exchange_factory_mod,        only : create_symmetric_halo_vec_exchange_Ch
use const_mod,                   only : pi
use particles_mod,               only : particles_t
use distribute_particles_mod,    only : distribute_particles_t
use geosci_config_mod,           only : geosci_config_t
use string_mod,                  only : strings, string_t
use particle_interp_factory_mod, only : create_particle_interp

implicit none

integer(kind=4), parameter :: margin = 9

contains

subroutine create_ecs_halo_vec_Ch(halo_vec, domain, halo_width)

    class(halo_vec_t), allocatable, intent(out) :: halo_vec
    type(domain_t),                 intent(in)  :: domain
    integer(kind=4),                intent(in)  :: halo_width

    type(ecs_halo_vec_Ch_t), allocatable :: ecs_halo

    integer(kind=4)    :: t, i, j, iface, nx, ny
    real(kind=8)       :: r(3), coords(5)
    type(particles_t) :: particles_tangent(domain%partition%ts:domain%partition%te)
    type(particles_t), allocatable :: particles_tangent_dist(:)
    type(distribute_particles_t) :: distribute
    type(geosci_config_t)        :: distribute_conf
    type(string_t)               :: str_levs(domain%partition%Nz)

    allocate(ecs_halo)

    ecs_halo%exchange = create_symmetric_halo_vec_exchange_Ch(domain%partition, domain%parcomm, domain%topology, halo_width, "full")

    do t = domain%partition%ts, domain%partition%te
        call particles_tangent(t)%init(strings("alpha","beta","ax", "ay", "az"),0)

        ny = domain%mesh_y%tile(t)%ny
        nx = domain%mesh_x%tile(t)%nx

        !'south' panel edge
        if(domain%mesh_y%tile(t)%js == 1) then
            do j = 1, halo_width
                do i = domain%mesh_y%tile(t)%is, domain%mesh_y%tile(t)%ie

                    r(1) = domain%mesh_y%tile(t)%rx(i,1-j,1)
                    r(2) = domain%mesh_y%tile(t)%ry(i,1-j,1)
                    r(3) = domain%mesh_y%tile(t)%rz(i,1-j,1)

                    call domain%metric%transform_cartesian_to_native(iface, coords(1), coords(2), r)
                    
                    coords(3) = domain%mesh_y%tile(t)%a1(1,i,1-j,1)
                    coords(4) = domain%mesh_y%tile(t)%a1(2,i,1-j,1)
                    coords(5) = domain%mesh_y%tile(t)%a1(3,i,1-j,1)

                    call particles_tangent(t)%add_particle(coords, iface)

                end do
            end do
        end if

        !'north' panel edge
        if(domain%mesh_y%tile(t)%je == ny) then
            do j = 1, halo_width
                do i = domain%mesh_y%tile(t)%is, domain%mesh_y%tile(t)%ie

                    r(1) = domain%mesh_y%tile(t)%rx(i,ny+j,1)
                    r(2) = domain%mesh_y%tile(t)%ry(i,ny+j,1)
                    r(3) = domain%mesh_y%tile(t)%rz(i,ny+j,1)

                    call domain%metric%transform_cartesian_to_native(iface, coords(1), coords(2), r)
                    coords(3) = domain%mesh_y%tile(t)%a1(1,i,ny+j,1)
                    coords(4) = domain%mesh_y%tile(t)%a1(2,i,ny+j,1)
                    coords(5) = domain%mesh_y%tile(t)%a1(3,i,ny+j,1)

                    call particles_tangent(t)%add_particle(coords, iface)

                end do
            end do
        end if

        !'west' panel edge
        if(domain%mesh_x%tile(t)%is == 1) then
            do i = 1, halo_width
                do j = domain%mesh_x%tile(t)%js, domain%mesh_x%tile(t)%je

                    r(1) = domain%mesh_x%tile(t)%rx(1-i,j,1)
                    r(2) = domain%mesh_x%tile(t)%ry(1-i,j,1)
                    r(3) = domain%mesh_x%tile(t)%rz(1-i,j,1)

                    call domain%metric%transform_cartesian_to_native(iface, coords(1), coords(2), r)
                    coords(3) = domain%mesh_x%tile(t)%a2(1,1-i,j,1)
                    coords(4) = domain%mesh_x%tile(t)%a2(2,1-i,j,1)
                    coords(5) = domain%mesh_x%tile(t)%a2(3,1-i,j,1)

                    call particles_tangent(t)%add_particle(coords, iface)

                end do
            end do
        end if

        !'east' panel edge
        if(domain%mesh_x%tile(t)%ie == nx) then
            do i = 1, halo_width
                do j = domain%mesh_x%tile(t)%js, domain%mesh_x%tile(t)%je

                    r(1) = domain%mesh_x%tile(t)%rx(nx+i,j,1)
                    r(2) = domain%mesh_x%tile(t)%ry(nx+i,j,1)
                    r(3) = domain%mesh_x%tile(t)%rz(nx+i,j,1)

                    call domain%metric%transform_cartesian_to_native(iface, coords(1), coords(2), r)
                    coords(3) = domain%mesh_x%tile(t)%a2(1,nx+i,j,1)
                    coords(4) = domain%mesh_x%tile(t)%a2(2,nx+i,j,1)
                    coords(5) = domain%mesh_x%tile(t)%a2(3,nx+i,j,1)

                    call particles_tangent(t)%add_particle(coords, iface)

                end do
            end do
        end if

        ny = domain%mesh_x%tile(t)%ny
        nx = domain%mesh_y%tile(t)%nx

        !'south' panel edge
        if(domain%mesh_x%tile(t)%js == 1) then
            do j = 1, halo_width
                do i = domain%mesh_x%tile(t)%is, domain%mesh_x%tile(t)%ie

                    r(1) = domain%mesh_x%tile(t)%rx(i,1-j,1)
                    r(2) = domain%mesh_x%tile(t)%ry(i,1-j,1)
                    r(3) = domain%mesh_x%tile(t)%rz(i,1-j,1)

                    call domain%metric%transform_cartesian_to_native(iface, coords(1), coords(2), r)
                
                    coords(3) = domain%mesh_x%tile(t)%a2(1,i,1-j,1)
                    coords(4) = domain%mesh_x%tile(t)%a2(2,i,1-j,1)
                    coords(5) = domain%mesh_x%tile(t)%a2(3,i,1-j,1)

                    call particles_tangent(t)%add_particle(coords, iface)

                end do
            end do
        end if

        !'north' panel edge
        if(domain%mesh_x%tile(t)%je == ny) then
            do j = 1, halo_width
                do i = domain%mesh_x%tile(t)%is, domain%mesh_x%tile(t)%ie

                    r(1) = domain%mesh_x%tile(t)%rx(i,ny+j,1)
                    r(2) = domain%mesh_x%tile(t)%ry(i,ny+j,1)
                    r(3) = domain%mesh_x%tile(t)%rz(i,ny+j,1)

                    call domain%metric%transform_cartesian_to_native(iface, coords(1), coords(2), r)
                    coords(3) = domain%mesh_x%tile(t)%a2(1,i,ny+j,1)
                    coords(4) = domain%mesh_x%tile(t)%a2(2,i,ny+j,1)
                    coords(5) = domain%mesh_x%tile(t)%a2(3,i,ny+j,1)

                    call particles_tangent(t)%add_particle(coords, iface)

                end do
            end do
        end if

        !'west' panel edge
        if(domain%mesh_y%tile(t)%is == 1) then
            do i = 1, halo_width
                do j = domain%mesh_y%tile(t)%js, domain%mesh_y%tile(t)%je

                    r(1) = domain%mesh_y%tile(t)%rx(1-i,j,1)
                    r(2) = domain%mesh_y%tile(t)%ry(1-i,j,1)
                    r(3) = domain%mesh_y%tile(t)%rz(1-i,j,1)

                    call domain%metric%transform_cartesian_to_native(iface, coords(1), coords(2), r)
                    coords(3) = domain%mesh_y%tile(t)%a1(1,1-i,j,1)
                    coords(4) = domain%mesh_y%tile(t)%a1(2,1-i,j,1)
                    coords(5) = domain%mesh_y%tile(t)%a1(3,1-i,j,1)

                    call particles_tangent(t)%add_particle(coords, iface)

                end do
            end do
        end if

        !'east' panel edge
        if(domain%mesh_y%tile(t)%ie == nx) then
            do i = 1, halo_width
                do j = domain%mesh_y%tile(t)%js, domain%mesh_y%tile(t)%je

                    r(1) = domain%mesh_y%tile(t)%rx(nx+i,j,1)
                    r(2) = domain%mesh_y%tile(t)%ry(nx+i,j,1)
                    r(3) = domain%mesh_y%tile(t)%rz(nx+i,j,1)

                    call domain%metric%transform_cartesian_to_native(iface, coords(1), coords(2), r)
                    coords(3) = domain%mesh_y%tile(t)%a1(1,nx+i,j,1)
                    coords(4) = domain%mesh_y%tile(t)%a1(2,nx+i,j,1)
                    coords(5) = domain%mesh_y%tile(t)%a1(3,nx+i,j,1)

                    call particles_tangent(t)%add_particle(coords, iface)

                end do
            end do
        end if
    end do

    call distribute_conf%set("communication_style","local")
    call distribute_conf%set("max_dist", halo_width)
    call distribute%init(domain, distribute_conf)
    call distribute%distribute_particles(particles_tangent_dist,particles_tangent, domain)

    allocate(ecs_halo%interp_u(domain%partition%ts:domain%partition%te))
    allocate(ecs_halo%interp_v(domain%partition%ts:domain%partition%te))
    allocate(ecs_halo%values_u(domain%partition%ts:domain%partition%te))
    allocate(ecs_halo%values_v(domain%partition%ts:domain%partition%te))
    allocate(ecs_halo%values_assembled(domain%partition%ts:domain%partition%te))

    allocate(ecs_halo%u_proj(domain%partition%ts:domain%partition%te))
    allocate(ecs_halo%v_proj(domain%partition%ts:domain%partition%te))

    do t = domain%partition%ts, domain%partition%te

        call create_particle_interp(ecs_halo%interp_u(t)%interp, "bicubic", particles_tangent_dist(t)%N, "y", domain)
        call create_particle_interp(ecs_halo%interp_v(t)%interp, "bicubic", particles_tangent_dist(t)%N, "x", domain)
        
        call ecs_halo%interp_u(t)%interp%calc_w(particles_tangent_dist(t), t, domain)
        call ecs_halo%interp_v(t)%interp%calc_w(particles_tangent_dist(t), t, domain)

        call ecs_halo%values_u(t)%init(particles_tangent_dist(t)%N,str_levs)
        call ecs_halo%values_v(t)%init(particles_tangent_dist(t)%N,str_levs)
        call ecs_halo%values_assembled(t)%init(particles_tangent(t)%N,str_levs)

        allocate(ecs_halo%u_proj(t)%a(particles_tangent_dist(t)%N))
        allocate(ecs_halo%v_proj(t)%a(particles_tangent_dist(t)%N))

        call init_proj_w(ecs_halo%u_proj(t)%a, ecs_halo%v_proj(t)%a, &
                         particles_tangent_dist(t), domain%metric,   &
                         domain%partition%panel_map(t))

    end do

    call ecs_halo%assembly%prepare_distribute_inversion(distribute, domain)

    call move_alloc(ecs_halo, halo_vec)

end subroutine

subroutine init_proj_w(wu, wv, particles, metric, face_ind)

    use metric_mod, only : metric_t

    real(kind=8), intent(inout) :: wu(:), wv(:)
    type(particles_t), intent(in) :: particles
    class(metric_t),   intent(in) :: metric
    integer(kind=4),   intent(in) :: face_ind

    integer(kind=4) :: i
    real(kind=8) :: a(3), b1(3), b2(3)
    real(kind=8) :: alpha, beta

    do i =1, particles%N

        alpha = particles%coords%a(1,i)
        beta  = particles%coords%a(2,i)
        a(1)  = particles%coords%a(3,i)
        a(2)  = particles%coords%a(4,i)
        a(3)  = particles%coords%a(5,i)

        b1(1:3) = metric%calculate_b1(face_ind, alpha, beta)
        b2(1:3) = metric%calculate_b2(face_ind, alpha, beta)

        wu(i) = sum(a(1:3)*b1(1:3))
        wv(i) = sum(a(1:3)*b2(1:3))

    end do

end subroutine

end module