module distribute_particles_mod

use particles_mod,      only : particles_t
use smart_array_mod,    only : smart_array_i4_t
use domain_mod,         only : domain_t
use mesh_mod,           only : mesh_t
use generic_config_mod, only : generic_config_t
use tile_mod,           only : tile_t
use tiles_mod,          only : tiles_t

implicit none

integer(kind=4), parameter :: BUFFER_SIZE=1000

type distribute_particles_t

    integer(kind=4) :: ts, te, N_tiles_send, N_proc_recv
    integer(kind=4) :: recv_buffer_size = 0
    integer(kind=4), allocatable :: local_to_real_tile(:), &
                                    real_to_local_tile(:)

    type(smart_array_i4_t), allocatable :: go_to_tile(:)
    type(smart_array_i4_t), allocatable :: come_from_proc(:), N_from_proc(:)

    contains
        procedure, public  :: init
        procedure, public  :: distribute_particles
        procedure, private :: init_sends
        procedure, private :: send_recv_particles
end type

contains

subroutine init(this, domain, config)
    class(distribute_particles_t), intent(inout) :: this
    type(domain_t),                intent(in)    :: domain
    class(generic_config_t),       intent(inout) :: config

    character(len=:), allocatable :: communication_style
    integer(kind=4) :: t, remote_tile, max_dist, proc
    logical :: if_intersect
    type(smart_array_i4_t) :: recv_from_proc, send_to_tile

    this%ts = domain%partition%ts
    this%te = domain%partition%te

    allocate(this%go_to_tile(this%ts:this%te))
    allocate(this%come_from_proc(this%ts:this%te))
    allocate(this%N_from_proc(this%ts:this%te))
        
    do t = this%ts, this%te
        call this%go_to_tile(t)%init(0)
        call this%come_from_proc(t)%init(0)
        call this%N_from_proc(t)%init(0)
    end do

    call config%get(this%recv_buffer_size, "recv_buffer_size", default=1)

    call config%get(communication_style, "communication_style", default="all-to-all")

    select case(communication_style)
    case("all-to-all")
        this%N_tiles_send = domain%partition%Nt
        this%N_proc_recv = domain%parcomm%np - 1 ! do not receive from myself

        allocate(this%local_to_real_tile(domain%partition%Nt))
        allocate(this%real_to_local_tile(domain%partition%Nt))

        do t = 1, domain%partition%Nt
            this%local_to_real_tile(t) = t
            this%real_to_local_tile(t) = t
        end do

    case("local")

        call config%get(max_dist,"max_dist")

        allocate(this%real_to_local_tile(1:domain%partition%Nt))
        this%real_to_local_tile(:) = 0 !just fill value
        this%N_tiles_send = 0

        call send_to_tile%init(0)

        !find tiles where we could possibly send some particles
        do remote_tile = 1, domain%partition%Nt
            do t = this%ts, this%te
                if(check_intersection(t, max_dist, remote_tile, domain)) then
                    call send_to_tile%append(remote_tile)
                    this%N_tiles_send = this%N_tiles_send+1
                    this%real_to_local_tile(remote_tile) = this%N_tiles_send
                    exit
                end if
            end do
        end do

        !find processes from where we could receive some particles
        this%N_proc_recv = 0
        call recv_from_proc%init(0)
        do remote_tile = 1, domain%partition%Nt
            do t = this%ts, this%te
                if(check_intersection(remote_tile, max_dist, t, domain)) then

                    proc = domain%partition%proc_map(remote_tile)
                    if(recv_from_proc%find_val(proc) == 0) then
                        call recv_from_proc%append(proc)
                        this%N_proc_recv = this%N_proc_recv+1
                    end if

                    exit
                end if
            end do
        end do

        this%N_proc_recv = this%N_proc_recv-1 !do not receive from myself
        this%local_to_real_tile = send_to_tile%a(1:send_to_tile%size_used)

    case default
        call domain%parcomm%abort("distribute particles init error, unknown communication style: "//communication_style)
    end select

    contains

    logical function check_intersection(tile_n1, halo_width, tile_n2, domain) result(intersection)

        !checks if tile_n1 extended by halo_width intersects with tile_n2

        integer(kind=4), intent(in) :: tile_n1, tile_n2, halo_width
        type(domain_t),  intent(in) :: domain

        integer(kind=4) :: panel1, panel2
        type(tiles_t) :: tiles
        type(tile_t) :: tile1, tile2

        call domain%partition%get_tiles("xy",tiles)

        panel1 = domain%partition%panel_map(tile_n1)
        panel2 = domain%partition%panel_map(tile_n2)

        call domain%topology%transform_tile_coords(panel1, tiles%tile(tile_n1), &
                                                   panel2, tile1,            &
                                                   tiles%Nx, tiles%Ny)

        tile2 = tiles%tile(tile_n2)

        tile1%is = max(tile1%is-halo_width, tile2%is)
        tile1%ie = min(tile1%ie+halo_width, tile2%ie)
        tile1%js = max(tile1%js-halo_width, tile2%js)
        tile1%je = min(tile1%je+halo_width, tile2%je)

        intersection = (tile1%is <= tile1%ie) .and. (tile1%js <= tile1%je)

    end function

end subroutine

subroutine distribute_particles(this, particles_out, particles_in, domain)

    !input
    class(distribute_particles_t), intent(inout) :: this
    type(domain_t),                intent(in)    :: domain
    type(particles_t),             intent(in)    :: particles_in(this%ts:this%te)

    !output:
    type(particles_t),      allocatable, intent(inout) :: particles_out(:)

    !local
    type(particles_t), allocatable :: particles_sends(:)
    integer(kind=4) :: t

    if( .not. allocated(particles_out) ) then

        allocate(particles_out(this%ts:this%te))
        
        do t = this%ts, this%te
            call particles_out(t)%init(particles_in(t)%coord_names, particles_in(t)%N)
        end do

    end if

    !cleanup work arrays from previous uses, i.e.
    !set current position in smart arrays to start:
    do t = this%ts, this%te
        call particles_out(t)%reset()
        call this%go_to_tile(t)%reset()
        call this%come_from_proc(t)%reset()
    end do

    call this%init_sends(particles_sends, particles_in, domain)

    call this%send_recv_particles(particles_out, particles_sends, domain)

end subroutine

subroutine init_sends(this, particles_sends, particles_in, domain)

    !input
    class(distribute_particles_t), intent(inout) :: this
    type(domain_t),    intent(in) :: domain
    type(particles_t), intent(in) :: particles_in(this%ts:this%te)

    !output:
    type(particles_t),      allocatable, intent(inout) :: particles_sends(:)

    integer(kind=4) :: i, t, panel_ind, indx, indy, to_tile, ic
    integer(kind=4) :: first_guess_size
    real(kind=8) :: alpha, beta
    type(mesh_t), pointer :: mesh

    first_guess_size = 0
    do t = domain%partition%ts, domain%partition%te
        first_guess_size = first_guess_size+particles_in(t)%N
    end do

    first_guess_size = 2*first_guess_size / this%N_tiles_send

    allocate(particles_sends(0:this%N_tiles_send))
    call particles_sends(0)%init(particles_in(this%ts)%coord_names, 0)
    do t = 1, this%N_tiles_send
        call particles_sends(t)%init(particles_in(this%ts)%coord_names, first_guess_size)
    end do

    call domain%get_mesh(mesh, "xy")

    do t = this%ts, this%te
        do i = 1, particles_in(t)%N

            alpha     = particles_in(t)%coords%a(1,i)
            beta      = particles_in(t)%coords%a(2,i)
            panel_ind = particles_in(t)%face_ind%a(i)

            indx = int((alpha - mesh%tile(t)%alpha_0) / mesh%tile(t)%hx+1 - mesh%tile(t)%shift_i)
            indy = int((beta  - mesh%tile(t)%beta_0)  / mesh%tile(t)%hy+1 - mesh%tile(t)%shift_j)
            
            to_tile = domain%partition%tile_map_xy(panel_ind)%a(indx,indy)

            call this%go_to_tile(t)%append(to_tile)

            to_tile = this%real_to_local_tile(to_tile)
            call particles_sends(to_tile)%add_particle(particles_in(t)%coords%a(:,i),panel_ind)

        end do
    end do

    if(particles_sends(0)%N > 0) then
        call domain%parcomm%abort("distribute particles, some particles fall out of local set of tiles. Consider increasing max_dist")
    end if

end subroutine

subroutine send_recv_particles(this, particles_out, &
                               particles_sends, domain)

    use mpi

    !input
    class(distribute_particles_t), intent(inout) :: this
    type(domain_t),    intent(in) :: domain
    type(particles_t), intent(in) :: particles_sends(0:this%N_tiles_send)

    !output:
    type(particles_t),      intent(inout) :: particles_out(domain%partition%ts:domain%partition%te)

    integer(kind=4) :: t, t2, dest_tile, i, n_coords
    integer(kind=4) :: send_count, ierr, dest, src, msg_size, n_pts
    integer(kind=4) :: send_req(domain%partition%Nt)
    integer(kind=4) :: status_send(MPI_STATUS_SIZE, domain%partition%Nt)
    integer(kind=4) :: status(MPI_STATUS_SIZE)

    real(kind=8),    allocatable :: buffer(:,:)
    integer(kind=4), allocatable :: panel_ind(:)

    n_coords = particles_sends(0)%n_coords

    allocate(buffer(n_coords, this%recv_buffer_size))
    allocate(panel_ind(this%recv_buffer_size))

    do t = domain%partition%ts, domain%partition%te

        t2 = this%real_to_local_tile(t)

        call particles_out(t)%add(particles_sends(t2))

        if(particles_sends(t2)%N /= 0) then
            call this%come_from_proc(t)%append(domain%partition%proc_map(t))
            call this%N_from_proc(t)%append(particles_sends(t2)%N)
        end if

    end do

    send_count = 0
    do t = 1, this%N_tiles_send

        t2 = this%local_to_real_tile(t)

        if(t2 >= domain%partition%ts .and. t2 <= domain%partition%te) cycle

        dest = domain%partition%proc_map(t2)
        send_count = send_count + 1

        call mpi_isend(particles_sends(t)%coords%a, &
                       n_coords*particles_sends(t)%N, &
                       MPI_DOUBLE, dest, t2, domain%parcomm%comm_w, &
                       send_req(send_count), ierr)

    end do

    do t = 1, (domain%partition%te-domain%partition%ts+1)*this%N_proc_recv

        call mpi_probe(MPI_ANY_SOURCE, MPI_ANY_TAG, domain%parcomm%comm_w, status, ierr)

        src = status(MPI_SOURCE)
        dest_tile = status(MPI_TAG)

        call mpi_get_count(status, MPI_DOUBLE, msg_size, ierr)
        n_pts = msg_size / n_coords

        if(n_pts > size(buffer,2)) then
            deallocate(buffer)
            allocate(buffer(n_coords,n_pts))
            deallocate(panel_ind)
            allocate(panel_ind(n_pts))
        end if

        call mpi_recv(buffer,msg_size,MPI_DOUBLE,src,dest_tile,&
                      domain%parcomm%comm_w,status,ierr)

        panel_ind(1:n_pts) = domain%partition%panel_map(dest_tile)

        call particles_out(dest_tile)%add_particles(buffer(1:n_coords,1:n_pts),panel_ind(1:n_pts))        

        if(n_pts /= 0) then
            call this%come_from_proc(dest_tile)%append(src)
            call this%N_from_proc(dest_tile)%append(n_pts)
        end if

    end do

    call mpi_waitall(send_count, send_req, status_send, ierr)

    call mpi_barrier(domain%parcomm%comm_w, ierr)

end subroutine

end module distribute_particles_mod   