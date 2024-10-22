module assembly_particles_mod

use smart_array_mod,          only : smart_array_i4_t
use particle_values_mod,      only : particle_values_t
use distribute_particles_mod, only : distribute_particles_t
use domain_mod,               only : domain_t
use parcomm_mod,              only : parcomm_t
use mpi

implicit none

integer(kind=4), parameter :: INITIAL_SIZE = 0

type, private :: recv_map_t

    integer(kind=4) :: N_points = 0
    type(smart_array_i4_t) :: to_tile, to_ind

end type

type, private :: send_map_t

    type(smart_array_i4_t) :: to_proc, N_points

end type

type, public :: assembly_particles_t

    integer(kind=4) :: ts, te
    integer(kind=4) :: N_recv, N_send
    integer(kind=4) :: max_recv_size, max_send_size

    integer(kind=4), allocatable :: local_to_real_tile(:), &
                                    real_to_local_tile(:)

    type(send_map_t), allocatable :: send_map(:)
    type(recv_map_t), allocatable :: map_recv_from_tile(:)

    contains
        procedure, public :: prepare_distribute_inversion
        procedure, public :: do => do_assembly
end type

contains

subroutine prepare_distribute_inversion(this, distribute, domain)

    class(assembly_particles_t),  intent(inout) :: this
    type(distribute_particles_t), intent(in)    :: distribute
    type(domain_t),               intent(in)    :: domain


    integer(kind=4) :: t, t_real_num, t_loc_num, i, n_points
    logical :: is_remote

    this%ts = domain%partition%ts
    this%te = domain%partition%te


    !Init send part (inverse recieve part of distribute)

    if(.not. allocated(this%send_map)) then

        allocate(this%send_map(this%ts:this%te))
  
    end if

    this%max_send_size = 0

    do t = this%ts, this%te
        this%send_map(t)%to_proc = distribute%come_from_proc(t)
        this%send_map(t)%N_points = distribute%N_from_proc(t)

        ! currently, we can not have 0-s in distribute%N_from_proc
        ! (but there is no guarantee, that this will not ever be violated) 
        ! so the following line can be uncommented:

        ! this%N_send = this%N_send + distribute%come_from_proc%size_used !(if line is uncommented, comment the line marked 'unnesesary check' below)

        do i = 1, distribute%come_from_proc(t)%size_used

            this%max_send_size = max(this%max_send_size, distribute%N_from_proc(t)%a(i))

            !unnesesary check:
            if(distribute%N_from_proc(t)%a(i) /= 0) this%N_send = this%N_send + 1

        end do

    end do

    !Init recieve part (inverse send part of distribute)

    if( .not. allocated (this%map_recv_from_tile)) then
        
        allocate(this%map_recv_from_tile(distribute%N_tiles_send))

        do t = 1, distribute%N_tiles_send
            call this%map_recv_from_tile(t)%to_tile%init(INITIAL_SIZE)
            call this%map_recv_from_tile(t)%to_ind%init(INITIAL_SIZE)
        end do

    else

        !clean-up from previous uses
        do t = 1, distribute%N_tiles_send
            call this%map_recv_from_tile(t)%to_tile%reset()
            call this%map_recv_from_tile(t)%to_ind%reset()
        end do

    end if

    if(allocated(this%real_to_local_tile)) deallocate(this%real_to_local_tile)
    allocate(this%real_to_local_tile, source=distribute%real_to_local_tile)

    do t = this%ts, this%te
        do i = 1, distribute%go_to_tile(t)%size_used

            t_real_num = distribute%go_to_tile(t)%a(i)
            t_loc_num  = this%real_to_local_tile(t_real_num)

            call this%map_recv_from_tile(t_loc_num)%to_tile%append(t)
            call this%map_recv_from_tile(t_loc_num)%to_ind%append(i)

        end do
    end do

    this%max_recv_size = 0

    do t = 1, distribute%N_tiles_send

        !check if considered tile is located at this MPI-proc
        t_real_num = distribute%local_to_real_tile(t)
        is_remote = (t_real_num < this%ts .or. t_real_num > this%te)

        n_points = this%map_recv_from_tile(t)%to_ind%size_used

        this%map_recv_from_tile(t)%N_points = n_points

        if(n_points /= 0 .and. is_remote) then
            this%N_recv = this%N_recv+1
            this%max_recv_size = max(this%max_recv_size, n_points)
        end if

    end do

end subroutine

subroutine do_assembly(this, values_out, values_in, domain)
 
    !input
    class(assembly_particles_t), intent(in) :: this
    type(domain_t),              intent(in) :: domain
    class(particle_values_t),    intent(in) :: values_in(domain%partition%ts:domain%partition%te)
    !output
    class(particle_values_t), intent(inout) :: values_out(domain%partition%ts:domain%partition%te)

    integer(kind=4) :: send_req(this%N_send)
    real(kind=8)    :: send_buffer(values_in(domain%partition%ts)%Nfld*this%max_send_size)
    real(kind=8)    :: recv_buffer(values_in(domain%partition%ts)%Nfld*this%max_recv_size)
    integer(kind=4) :: t, t_loc_num, t2, i, j, ind, ifld, send_count, val_count, &
                       start, dest_proc, ierr, Nfld, status(MPI_STATUS_SIZE), src

    Nfld = values_in(domain%partition%ts)%Nfld

    send_count = 0
    do t = domain%partition%ts, domain%partition%te

        t_loc_num = this%real_to_local_tile(t) !use translated number for recv maps

        start = 0
        do i = 1, this%send_map(t)%to_proc%size_used

            dest_proc = this%send_map(t)%to_proc%a(i)

            !avoid self sends
            if(dest_proc == domain%parcomm%myid) then

                do ifld = 1, Nfld
                    do j = 1, this%send_map(t)%N_points%a(i)
                        t2 = this%map_recv_from_tile(t_loc_num)%to_tile%a(j)
                        ind = this%map_recv_from_tile(t_loc_num)%to_ind%a(j)
                        values_out(t2)%p(ind,ifld) = values_in(t)%p(start+j,ifld)
                    end do
                end do

            else

                val_count = 0
                do ifld = 1, Nfld
                    do j = 1, this%send_map(t)%N_points%a(i)
                        val_count = val_count+1
                        send_buffer(val_count) = values_in(t)%p(start+j,ifld)
                    end do

                end do

                send_count = send_count+1
                call mpi_isend(send_buffer, val_count, MPI_DOUBLE, dest_proc, t, &
                               domain%parcomm%comm_w, send_req(send_count), ierr)
            end if

            start = start + this%send_map(t)%N_points%a(i)
        end do

    end do

    do i = 1, this%N_recv

        call mpi_probe(MPI_ANY_SOURCE, MPI_ANY_TAG, domain%parcomm%comm_w, status, ierr)

        t   = status(MPI_TAG)
        t_loc_num = this%real_to_local_tile(t)
        src = status(MPI_SOURCE)

        call mpi_recv(recv_buffer,this%map_recv_from_tile(t_loc_num)%N_points*Nfld, &
                      MPI_DOUBLE, src, t, domain%parcomm%comm_w, status, ierr)

        val_count = 1
        do ifld = 1, Nfld
            do j = 1, this%map_recv_from_tile(t_loc_num)%N_points
                t2 = this%map_recv_from_tile(t_loc_num)%to_tile%a(j)
                ind = this%map_recv_from_tile(t_loc_num)%to_ind%a(j)
                values_out(t2)%p(ind,ifld) = recv_buffer(val_count)
                val_count = val_count+1
            end do
        end do

    end do

    call mpi_waitall(send_count, send_req, status, ierr)
    
    call mpi_barrier(domain%parcomm%comm_w, ierr)

end subroutine

end module