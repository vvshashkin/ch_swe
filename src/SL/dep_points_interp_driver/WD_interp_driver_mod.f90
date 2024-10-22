! this module implements White III -Dongarra algorithm 
! for interpolation of fields to departure points
! http://dx.doi.org/10.1016/j.jcp.2011.05.008

module WD_interp_driver_mod

use abstract_dep_points_interp_driver_mod,  only : dep_points_interp_driver_t
use abstract_dp_interpolator_mod,           only : dp_interpolator_t, dp_interpolator_container_t
use stvec_flexible_mod,                     only : stvec_flexible_t
use domain_mod,                             only : domain_t
use string_mod,                             only : string_t
use grid_field_mod,                         only : grid_field_t
use mesh_mod,                               only : mesh_t, tile_mesh_t
use parcomm_mod,                            only : parcomm_global
use mpi

implicit none


type, private :: int_2d_array_t
    integer(kind=4), allocatable :: ind(:,:)
    ! contains
    !     procedure :: init => init_int_2d_array
end type

type, private :: buffer_to_grid_map_t !maps 1d recvs buffer to 3d gridfields indices
    integer(kind=4), dimension(:), allocatable :: i, j, k, t
    integer(kind=4), dimension(:), allocatable :: start_of_tile !starting index of i, j, k, t for recv from tile
    ! contains
    !     procedure :: init => init_buffer_to_grid_map
end type

type, extends(dep_points_interp_driver_t) :: WD_interp_driver_t

    !separate communicators for departure points and interpolated values exchanges:
    integer(kind=mpi_integer_kind) :: dp_comm, val_comm 

    type(int_2d_array_t), allocatable :: tile_map(:) !establish correspondence between horizontal grid indices and tile number

    type(mesh_t), pointer :: mesh

    integer(kind=4) :: num_sends_dp, num_recvs_dp !number of send and recv calls to exchange departure points
    integer(kind=4), allocatable :: send_dp_to_tiles(:)

    real(kind=8), allocatable :: dp_send_buffer(:)
    integer(kind=4), allocatable :: dp_send_request(:)
    integer(kind=4) :: n_my_points
    type(buffer_to_grid_map_t) :: gather_map

    real(kind=8), allocatable :: val_buffer(:), dp_recv_buffer(:)
    integer(kind=4), allocatable :: val_send_request(:)

    type(dp_interpolator_container_t), allocatable :: interpolator_array(:)
    type(string_t), allocatable :: interp_fields(:)
    integer(kind=4)             :: n_fld 
    integer(kind=4), allocatable :: interpolator_index(:)

    
    contains
    procedure :: do
    procedure, private :: interp_halo
    procedure, private :: send_dp_and_prepare_gather
    procedure, private :: recv_dp_interpolate_send_vals
    procedure, private :: recv_interpolated_values
    procedure, private :: cleanup_mpi_requests
end type WD_interp_driver_t

contains

subroutine do(this, q_dep, q, dep_points_coords, domain)
    class(WD_interp_driver_t),  intent(inout) :: this
    class(stvec_flexible_t),    intent(inout) :: q_dep, q, dep_points_coords
    type(domain_t),             intent(in)    :: domain

    type(mesh_t), pointer :: mesh

    call this%interp_halo(q, domain)

    call this%send_dp_and_prepare_gather(dep_points_coords, this%mesh, domain)

    call this%recv_dp_interpolate_send_vals(q, domain)

    call this%recv_interpolated_values(q_dep, domain)

    call this%cleanup_mpi_requests()

end subroutine  

subroutine interp_halo(this, q, domain)
    class(WD_interp_driver_t),  intent(inout) :: this
    class(stvec_flexible_t),    intent(inout) :: q
    type(domain_t),             intent(in)    :: domain

    type(grid_field_t), pointer :: f
    integer(kind=4) :: i, ind

    do i = 1, size(this%interp_fields,1)
        call q%get_field(f,this%interp_fields(i)%str)
        ind = this%interpolator_index(i)
        call this%interpolator_array(ind)%interpolator%ext_halo(f,domain)
    end do

end subroutine

subroutine send_dp_and_prepare_gather(this,dep_points_coords,mesh,domain)

    class(WD_interp_driver_t),  intent(inout) :: this
    class(stvec_flexible_t),    intent(inout) :: dep_points_coords
    type(mesh_t),               intent(in)    :: mesh
    type(domain_t),             intent(in)    :: domain

    real(kind=8), dimension(this%n_my_points) :: alpha, beta, eta
    integer(kind=4) :: send_map(this%n_my_points)
    integer(kind=4), dimension(0:this%num_sends_dp) :: num_dp_on_tile, num_dp_on_tile_tmp, msg_start, msg_end
    integer(kind=4) :: i,j,k,t, panel_ind, indx, indy, to_tile, count, ierr, end_of_tile
    real(kind=8) :: r(3)
    type(grid_field_t), pointer :: alpha_dp, beta_dp, panel_ind_dp, eta_dp

    call dep_points_coords%get_field(alpha_dp,"alpha")
    call dep_points_coords%get_field(beta_dp,"beta")
    call dep_points_coords%get_field(panel_ind_dp,"panel_ind")
    call dep_points_coords%get_field(eta_dp,"eta")

    ! 1.1 Look up tile where each departure point is located and count dp for each tile
    !     (not exactly departure point itself, but the base point of interpolation stencil)
    num_dp_on_tile(0:this%num_sends_dp) = 0
    count = 0
    do t = mesh%ts, mesh%te
        do k = mesh%tile(t)%ks, mesh%tile(t)%ke
            do j = mesh%tile(t)%js, mesh%tile(t)%je
                do i = mesh%tile(t)%is, mesh%tile(t)%ie

                    count = count+1

                    ! r(1:3) = [x%tile(t)%p(i,j,k),y%tile(t)%p(i,j,k),z%tile(t)%p(i,j,k)]
                    ! call domain%metric%transform_cartesian_to_native(panel_ind, &
                        !  alpha(count),beta(count),r)
                    ! eta(count) = 0.0_8

                    alpha(count) = alpha_dp%tile(t)%p(i,j,k)
                    beta(count)  = beta_dp%tile(t)%p(i,j,k)
                    panel_ind = int(panel_ind_dp%tile(t)%p(i,j,k))

                    ! find interpolation stencil base point (indx, indy) and tile number (to_tile) it belongs to
                    ! remark: to_tile will be 0 (see dep_points_interp_driver_factory_mod.f90)
                    ! if point indx, indy falls outside the set of expected tiles, this will be used later as abort condition
                    indx = int((alpha(count) - mesh%tile(t)%alpha_0) / mesh%tile(t)%hx+1 - mesh%tile(t)%shift_i)
                    indy = int((beta(count)  - mesh%tile(t)%beta_0)  / mesh%tile(t)%hy+1 - mesh%tile(t)%shift_j)
                    to_tile = this%tile_map(panel_ind)%ind(indx,indy)

                    send_map(count) = to_tile
                    num_dp_on_tile(to_tile) = num_dp_on_tile(to_tile)+1

                end do
            end do
        end do
    end do

    !assign eta
    count = 1
    if(associated(eta_dp)) then
        do t = mesh%ts, mesh%te
            do k = mesh%tile(t)%ks, mesh%tile(t)%ke
                do j = mesh%tile(t)%js, mesh%tile(t)%je
                    do i = mesh%tile(t)%is, mesh%tile(t)%ie
                        eta(count) = eta_dp%tile(t)%p(i,j,k)
                        count = count+1
                    end do
                end do
            end do
        end do
    else
        do t = mesh%ts, mesh%te
            do k = mesh%tile(t)%ks, mesh%tile(t)%ke
                do j = mesh%tile(t)%js, mesh%tile(t)%je
                    do i = mesh%tile(t)%is, mesh%tile(t)%ie
                        eta(count) = mesh%tile(t)%get_eta(k)
                        count = count+1
                    end do
                end do
            end do
        end do
    end if

    ! 1.2 Check if any departure point is located outside expected set of tiles
    ! (these were previously assigned to tile 0)
    if(num_dp_on_tile(0) > 0) then
        call parcomm_global%abort("Departure points out of expected tiles. Check SL MAX_CFL")
    end if

    ! 2. Sort departure points by the tiles they belong to
    ! 2.1 Mark up send buffer and gather map by tiles
    msg_start(1) = 1
    msg_end(1)   = 3*num_dp_on_tile(1)
    this%gather_map%start_of_tile(1) = 1
    end_of_tile = num_dp_on_tile(1)
    do t = 2, ubound(num_dp_on_tile,1)
        msg_start(t) = msg_end(t-1)+1
        msg_end(t) = msg_end(t-1)+3*num_dp_on_tile(t)
        this%gather_map%start_of_tile(t) = end_of_tile+1
        end_of_tile = end_of_tile+num_dp_on_tile(t)
    end do

    ! 2.2 Place departure points coordinates in send buffer and
    !     place indices of points in gather map to distribute message
    !     with interpolated values that we will have in reply
    count = 0
    num_dp_on_tile_tmp(0:this%num_sends_dp) = 0
    do t = mesh%ts, mesh%te
        do k = mesh%tile(t)%ks, mesh%tile(t)%ke
            do j = mesh%tile(t)%js, mesh%tile(t)%je
                do i = mesh%tile(t)%is, mesh%tile(t)%ie

                    count = count +1

                    to_tile = send_map(count)
                    num_dp_on_tile_tmp(to_tile) = num_dp_on_tile_tmp(to_tile)+1
                    indx = msg_start(to_tile)+num_dp_on_tile_tmp(to_tile)-1

                    this%dp_send_buffer(indx) = alpha(count)
                    this%dp_send_buffer(indx+num_dp_on_tile(to_tile)) = beta(count)
                    this%dp_send_buffer(indx+2*num_dp_on_tile(to_tile)) = eta(count)

                    indx = this%gather_map%start_of_tile(to_tile)+num_dp_on_tile_tmp(to_tile)-1
                    this%gather_map%i(indx) = i
                    this%gather_map%j(indx) = j
                    this%gather_map%k(indx) = k
                    this%gather_map%t(indx) = t

                end do
            end do
        end do
    end do

    ! 2.3 Send departure points to tiles. TODO: avoid self-sends
    do t = 1, this%num_sends_dp
        to_tile = this%send_dp_to_tiles(t)
        call mpi_isend(this%dp_send_buffer(msg_start(t)),  &
                       3*num_dp_on_tile(t), MPI_DOUBLE,                   &
                       domain%partition%proc_map(to_tile), to_tile,   &
                       this%dp_comm, this%dp_send_request(t), ierr)
    end do

end subroutine send_dp_and_prepare_gather

subroutine recv_dp_interpolate_send_vals(this, q, domain)

    class(WD_interp_driver_t),  intent(inout) :: this
    class(stvec_flexible_t),    intent(inout) :: q
    type(domain_t),             intent(in)    :: domain

    integer(kind=4) :: msg_size, imsg, tile_ind, tag, source_rank, i, n_dp, st
    integer(kind=4) :: ierr, status(MPI_STATUS_SIZE)
    integer(kind=4) :: n_interp, n_fld, i_interp
    real(kind=8) :: x, y, z
    type(grid_field_t), pointer :: f

    st = 1
    do imsg = 1, this%num_recvs_dp

        call mpi_probe(MPI_ANY_SOURCE,MPI_ANY_TAG,this%dp_comm,status,ierr)

        tile_ind =  status(MPI_TAG)
        source_rank = status(MPI_SOURCE)
        call mpi_get_count(status,MPI_DOUBLE,msg_size,ierr)

        call mpi_recv(this%dp_recv_buffer,msg_size,MPI_DOUBLE, &
                      source_rank,tile_ind,this%dp_comm,status,ierr)


        n_dp = msg_size/3
        
        do i = 1, size(this%interpolator_array,1)
            call this%interpolator_array(i)%&
                 interpolator%calc_weights(this%dp_recv_buffer(1:n_dp),         &
                                           this%dp_recv_buffer(n_dp+1:2*n_dp),  &
                                           this%dp_recv_buffer(2*n_dp+1:3*n_dp),&
                                           tile_ind)
        end do

        n_fld = size(this%interp_fields,1)
        do i = 1, n_fld
            call q%get_field(f,this%interp_fields(i)%str)
            i_interp = this%interpolator_index(i)
            call this%interpolator_array(i_interp)%&
                 interpolator%interpolate(this%val_buffer(st+(i-1)*n_dp:st+i*n_dp-1),&
                                          f%tile(tile_ind))
        end do
        call mpi_isend(this%val_buffer(st),n_fld*n_dp,MPI_DOUBLE,source_rank, &
                       tile_ind,this%val_comm,this%val_send_request(imsg), ierr)

        st = st+n_fld*n_dp

    end do

end subroutine recv_dp_interpolate_send_vals

subroutine recv_interpolated_values(this,q_dep,domain)

    class(WD_interp_driver_t),  intent(inout) :: this
    class(stvec_flexible_t),    intent(inout) :: q_dep
    type(domain_t),             intent(in)    :: domain

    integer(kind=4) :: msg_size, expected_msg_count, imsg, tile_ind, tag, source_rank, n_dp, st
    integer(kind=4) :: i, j, k, t, count, count2, ifld, n_fld
    integer(kind=4) :: ierr, status(MPI_STATUS_SIZE)
    type(grid_field_t), pointer :: q
    real(kind=8) :: buffer(this%n_my_points*this%n_fld)

    n_fld = size(this%interp_fields,1)

    do imsg = 1, this%num_sends_dp
        call mpi_probe(MPI_ANY_SOURCE,MPI_ANY_TAG,this%val_comm,status,ierr)

        tile_ind =  status(MPI_TAG)
        source_rank = status(MPI_SOURCE)
        call mpi_get_count(status,MPI_DOUBLE,msg_size,ierr)

        call mpi_recv(buffer,msg_size,MPI_DOUBLE, &
                      source_rank,tile_ind,this%val_comm,status,ierr)

        tile_ind = findloc(this%send_dp_to_tiles,tile_ind,1)
        st = this%gather_map%start_of_tile(tile_ind)
        count2 = 1
        do ifld = 1, n_fld
            call q_dep%get_field(q,this%interp_fields(ifld)%str)
            do count = 1, msg_size/n_fld
                t = this%gather_map%t(st+count-1)
                i = this%gather_map%i(st+count-1)
                j = this%gather_map%j(st+count-1)
                k = this%gather_map%k(st+count-1)
                q%tile(t)%p(i,j,k) = buffer(count2)
                count2 = count2+1
            end do
        end do
    end do

end subroutine

subroutine cleanup_mpi_requests(this)

    class(WD_interp_driver_t),  intent(inout) :: this

    integer(kind=4) :: status_val(MPI_STATUS_SIZE,size(this%val_send_request,1))
    integer(kind=4) :: status_dp(MPI_STATUS_SIZE,size(this%dp_send_request,1))
    integer(kind=4) :: ierr

    call mpi_waitall(size(this%dp_send_request,1),this%dp_send_request,status_dp,ierr)
    call mpi_waitall(size(this%val_send_request,1),this%val_send_request,status_val,ierr)

end subroutine

end module