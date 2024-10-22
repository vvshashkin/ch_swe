module latlon_regrid_mod

use abstract_regrid_mod,          only : regrid_t, regrid_vec_t
use domain_mod,                   only : domain_t
use grid_field_mod,               only : grid_field_t, tile_field_t
use mesh_mod,                     only : mesh_t, tile_mesh_t
use parcomm_mod,                  only : parcomm_global
use halo_mod,                     only : halo_t, halo_vec_t
use abstract_interpolators2d_mod, only : interpolator2d_vec2vec_t
use mpi    

implicit none


type, public :: interp_t
    integer(kind=4)               :: nhor, nz
    real(kind=8),     allocatable :: wx(:,:), wy(:,:)
    integer(kind=4),  allocatable :: x_ind(:), y_ind(:)

    contains
    procedure, public :: init => init_interpolation, interpolate

end type interp_t

type, public :: real_2d_array_t
    real(kind=8), allocatable :: a(:,:)
end type
type, public :: real_3d_array_t
    real(kind=8), allocatable :: a(:,:,:)
end type

type, public :: indices_container_t
    integer(kind=4) :: n !size of arrays
    integer(kind=4), allocatable :: i(:), j(:)
end type

type, public, extends(regrid_t) :: latlon_regrid_t

    class(halo_t),                   allocatable :: halo
    type(grid_field_t)                           :: buff
    character(len=:),                allocatable :: mesh_name
    type(interp_t),                  allocatable :: latlon_interp(:)
    type(real_2d_array_t),           allocatable :: values(:)
    type(indices_container_t),       allocatable :: msg_to_latlon(:) !translate message 1d index to latlon ij indices

    integer(kind=4) :: nz, halo_width ! buffer size for recive at master

    contains

    procedure :: do_regrid     => do_latlon_scalar_regrid

end type latlon_regrid_t

type, public, extends(regrid_vec_t) :: latlon_regrid_vec_t

    class(interpolator2d_vec2vec_t), allocatable :: staggered2colocated
    class(halo_vec_t),               allocatable :: halo
    type(grid_field_t)                           :: buff_u, buff_v, buff_u_stag, buff_v_stag
    type(mesh_t),                    pointer     :: mesh, mesh_u, mesh_v
    type(interp_t),                  allocatable :: latlon_interp(:)
    type(real_3d_array_t),           allocatable :: values_uv(:)
    type(indices_container_t),       allocatable :: msg_to_latlon(:)
    !transform matrix from native uv components to latlon uv:
    type(real_3d_array_t),           allocatable :: q(:)

    integer(kind=4) :: max_points_at_tile, nz, halo_width ! buffer size for recive at master

    contains

    procedure :: do_regrid     => do_latlon_vector_regrid

end type latlon_regrid_vec_t
    
contains
    
subroutine do_latlon_scalar_regrid(this,fout,f,domain)

    class(latlon_regrid_t),  intent(inout) :: this
    real(kind=8),            intent(inout) :: fout(:,:,:)
    type(grid_field_t),      intent(in)    :: f
    type(domain_t),          intent(in)    :: domain
    
    integer(kind=4) :: t, ihor, k, i, j, ierr, proc_id, send_size
    integer(kind=4) :: request(domain%partition%ts:domain%partition%te), &
                       status(MPI_STATUS_SIZE, domain%partition%ts:domain%partition%te), &
                       status_recv(MPI_STATUS_SIZE)
    
    type(real_2d_array_t), allocatable :: recv_buffers(:)
    integer(kind=4) :: recv_request(domain%partition%Nt), &
                       recv_status(MPI_STATUS_SIZE), recv_from_tile(domain%partition%Nt)
    integer(kind=4) :: n_recv, i_recv, req_ind

    type(mesh_t), pointer :: mesh

    call domain%get_mesh(mesh,this%mesh_name)

    call this%buff%assign(f, mesh)
    call this%halo%get_halo_scalar(this%buff,domain,this%halo_width)
    
    do t = domain%partition%ts, domain%partition%te
        call this%latlon_interp(t)%interpolate(this%values(t)%a,this%buff%tile(t))

        if(domain%parcomm%myid == 0) then
            do k = 1, this%nz
                do ihor = 1, this%latlon_interp(t)%nhor
                    i = this%msg_to_latlon(t)%i(ihor)
                    j = this%msg_to_latlon(t)%j(ihor)
                    fout(i,j,k) = this%values(t)%a(ihor,k)
                end do
            end do
        else
            send_size = this%latlon_interp(t)%nhor*this%latlon_interp(t)%nz
            call mpi_isend(this%values(t)%a, send_size, MPI_DOUBLE, 0, t, &
                           domain%parcomm%comm_w, request(t), ierr)
        end if
    end do

    if(domain%parcomm%myid /= 0) then
        call mpi_waitall(size(request,1),request,status,ierr)
    else
        
        allocate(recv_buffers(domain%partition%Nt))

        n_recv = 0

        do t = 1, domain%partition%Nt

            if(t >= domain%partition%ts .and. t <= domain%partition%te) cycle

            proc_id = domain%partition%proc_map(t)

            allocate(recv_buffers(t)%a(this%msg_to_latlon(t)%n,this%nz))

            n_recv = n_recv+1

            call mpi_irecv(recv_buffers(t)%a, this%msg_to_latlon(t)%n*this%nz,&
                           MPI_DOUBLE, proc_id, t, domain%parcomm%comm_w,     &
                           recv_request(n_recv), ierr)

            recv_from_tile(n_recv) = t

        end do

        do i_recv = 1, n_recv

            call mpi_waitany(n_recv, recv_request, req_ind, status_recv, ierr)
            
            t = recv_from_tile(req_ind)

            do k = 1, this%nz
                do ihor = 1, this%msg_to_latlon(t)%n
                    i = this%msg_to_latlon(t)%i(ihor)
                    j = this%msg_to_latlon(t)%j(ihor)
                    fout(i,j,k) = recv_buffers(t)%a(ihor,k)
                end do
            end do
        end do

    end if

end subroutine do_latlon_scalar_regrid

subroutine do_latlon_vector_regrid(this,uout,vout,u,v,domain)

    class(latlon_regrid_vec_t),  intent(inout) :: this
    real(kind=8),                intent(inout) :: uout(:,:,:), vout(:,:,:)
    type(grid_field_t),          intent(in)    :: u, v
    type(domain_t),              intent(in)    :: domain
    
    integer(kind=4) :: t, ihor, k, i, j, ierr, proc_id, send_size
    integer(kind=4) :: request(domain%partition%ts:domain%partition%te), &
                       status(MPI_STATUS_SIZE, domain%partition%ts:domain%partition%te), &
                       status_recv(MPI_STATUS_SIZE)
    
    real(kind=8)    :: u_tmp, v_tmp
    type(real_3d_array_t), allocatable :: recv_buffers(:)

    integer(kind=4) :: recv_request(domain%partition%Nt), &
                       recv_status(MPI_STATUS_SIZE), recv_from_tile(domain%partition%Nt)
    integer(kind=4) :: n_recv, i_recv, req_ind

    if(.not. allocated(this%staggered2colocated)) then
        call this%buff_u%assign(u, this%mesh)
        call this%buff_v%assign(v, this%mesh)
    else
        call this%buff_u_stag%assign(u, this%mesh_u)
        call this%buff_v_stag%assign(v, this%mesh_v)
        call this%staggered2colocated%interp2d_vec2vec(this%buff_u, this%buff_v, &
                                                       this%buff_u_stag, this%buff_v_stag, &
                                                       domain)
    end if

    call this%halo%get_halo_vector(this%buff_u, this%buff_v, domain, this%halo_width)
    
    do t = domain%partition%ts, domain%partition%te

        call this%latlon_interp(t)%interpolate(this%values_uv(t)%a(:,:,1),this%buff_u%tile(t))
        call this%latlon_interp(t)%interpolate(this%values_uv(t)%a(:,:,2),this%buff_v%tile(t))

        !transform native components to lat-lon
        do k = 1, this%nz
            do ihor = 1, this%latlon_interp(t)%nhor
                u_tmp = this%q(t)%a(ihor,1,1)*this%values_uv(t)%a(ihor,k,1) + &
                        this%q(t)%a(ihor,1,2)*this%values_uv(t)%a(ihor,k,2)
                v_tmp = this%q(t)%a(ihor,2,1)*this%values_uv(t)%a(ihor,k,1) + &
                        this%q(t)%a(ihor,2,2)*this%values_uv(t)%a(ihor,k,2)
                this%values_uv(t)%a(ihor,k,1) = u_tmp
                this%values_uv(t)%a(ihor,k,2) = v_tmp
            end do
        end do

        if(domain%parcomm%myid == 0) then
            do k = 1, this%nz
                do ihor = 1, this%latlon_interp(t)%nhor
                    i = this%msg_to_latlon(t)%i(ihor)
                    j = this%msg_to_latlon(t)%j(ihor)
                    uout(i,j,k) = this%values_uv(t)%a(ihor,k,1)
                    vout(i,j,k) = this%values_uv(t)%a(ihor,k,2)
                end do
            end do
        else
            send_size = 2*this%latlon_interp(t)%nhor*this%latlon_interp(t)%nz
            call mpi_isend(this%values_uv(t)%a, send_size, MPI_DOUBLE, 0, t, &
                           domain%parcomm%comm_w, request(t), ierr)
        end if
    end do

    if(domain%parcomm%myid /= 0) then
        call mpi_waitall(size(request,1),request,status,ierr)
    else

        allocate(recv_buffers(domain%partition%Nt))

        n_recv = 0

        do t = 1, domain%partition%num_panels*domain%partition%num_tiles

            if(t >= domain%partition%ts .and. t <= domain%partition%te) cycle

            proc_id = domain%partition%proc_map(t)

            allocate(recv_buffers(t)%a(this%msg_to_latlon(t)%n,this%nz,2))

            n_recv = n_recv+1

            call mpi_irecv(recv_buffers(t)%a, 2*this%msg_to_latlon(t)%n*this%nz, MPI_DOUBLE, proc_id, &
                          t, domain%parcomm%comm_w, recv_request(n_recv), ierr)
            
            recv_from_tile(n_recv) = t

        end do

        do i_recv = 1, n_recv

            call mpi_waitany(n_recv, recv_request, req_ind, status_recv, ierr)
            
            t = recv_from_tile(req_ind)
            
            do k = 1, this%nz
                do ihor = 1, this%msg_to_latlon(t)%n
                    i = this%msg_to_latlon(t)%i(ihor)
                    j = this%msg_to_latlon(t)%j(ihor)
                    uout(i,j,k) = recv_buffers(t)%a(ihor,k,1)
                    vout(i,j,k) = recv_buffers(t)%a(ihor,k,2)
                end do
            end do

        end do
    end if

end subroutine do_latlon_vector_regrid

subroutine init_interpolation(this, alpha, beta, n_pts, mesh, interp_type)

    class(interp_t),   intent(inout) :: this
    real(kind=8),      intent(in)    :: alpha(:), beta(:)
    integer(kind=4),   intent(in)    :: n_pts
    type(tile_mesh_t), intent(in)    :: mesh
    character(len=*),  intent(in)    :: interp_type

    integer(kind=4) :: ihor
    real(kind=8) :: zdx, zdy

    this%nhor = n_pts
    this%nz = mesh%nz
    allocate(this%x_ind(n_pts), this%y_ind(n_pts))
    
    if(interp_type == "linear") then

        allocate(this%wx(0:1,n_pts), this%wy(0:1,n_pts))
        do ihor = 1, n_pts

            zdx = (alpha(ihor)- mesh%alpha_0) /mesh%hx+1-mesh%shift_i
            zdy = (beta(ihor) - mesh%beta_0)  /mesh%hy+1-mesh%shift_j
            this%x_ind(ihor) = max(1,min(mesh%nx-1,int(zdx)))
            this%y_ind(ihor) = max(1,min(mesh%ny-1,int(zdy)))
            zdx = zdx - this%x_ind(ihor)
            zdy = zdy - this%y_ind(ihor)
    
            this%wx(0,ihor) = 1.0_8-zdx
            this%wx(1,ihor) = zdx
            this%wy(0,ihor) = 1.0_8-zdy
            this%wy(1,ihor) = zdy
    
        end do

    else if (interp_type == "cubic") then

        allocate(this%wx(-1:2,n_pts), this%wy(-1:2,n_pts))

        do ihor = 1, n_pts

            zdx = (alpha(ihor)- mesh%alpha_0) /mesh%hx+1-mesh%shift_i
            zdy = (beta(ihor) - mesh%beta_0)  /mesh%hy+1-mesh%shift_j
            this%x_ind(ihor) = int(zdx)
            this%y_ind(ihor) = int(zdy)
            zdx = zdx - this%x_ind(ihor)
            zdy = zdy - this%y_ind(ihor)
    
            this%wx(-1,ihor) =-zdx*(zdx-1.0_8)*(zdx-2.0_8) / 6.0_8
            this%wx( 0,ihor) = (zdx+1.0_8)*(zdx-1.0_8)*(zdx-2.0_8) / 2.0_8
            this%wx( 1,ihor) =-(zdx+1.0_8)*zdx*(zdx-2.0_8) / 2.0_8
            this%wx( 2,ihor) = (zdx+1.0_8)*zdx*(zdx-1.0_8) / 6.0_8

            this%wy(-1,ihor) =-zdy*(zdy-1.0_8)*(zdy-2.0_8) / 6.0_8
            this%wy( 0,ihor) = (zdy+1.0_8)*(zdy-1.0_8)*(zdy-2.0_8) / 2.0_8
            this%wy( 1,ihor) =-(zdy+1.0_8)*zdy*(zdy-2.0_8) / 2.0_8
            this%wy( 2,ihor) = (zdy+1.0_8)*zdy*(zdy-1.0_8) / 6.0_8
    
        end do
    else
        call parcomm_global%abort("latlon regrid - unknown interpolation type: "// interp_type)

    end if

end subroutine

subroutine interpolate(this,fout,f)

    class(interp_t),    intent(in)      :: this
    real(kind=8),       intent(inout)   :: fout(this%nhor,this%nz)
    type(tile_field_t), intent(in) :: f

    integer(kind=4) :: i, j, k, ihor, i1, j1, i2, j2
    real(kind=8)    :: p

    i1 = lbound(this%wx,1)
    i2 = ubound(this%wx,1)
    j1 = lbound(this%wy,1)
    j2 = ubound(this%wy,1)

    do k = 1, this%nz
        do ihor = 1, this%nhor
            fout(ihor,k) = 0.0_8
            do j= j1, j2
                do i = i1, i2
                    p = f%p(this%x_ind(ihor)+i,this%y_ind(ihor)+j,k)
                    fout(ihor,k) = fout(ihor,k)+this%wx(i,ihor)*this%wy(j,ihor)*p
                end do
            end do
        end do
    end do

end subroutine interpolate

end module latlon_regrid_mod