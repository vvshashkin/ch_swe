module dp_interpolators_mod

use abstract_dp_interpolator_mod, only : dp_interpolator_t
use mesh_mod,                     only : mesh_t, tile_mesh_t
use grid_field_mod,               only : grid_field_t, tile_field_t
use exchange_abstract_mod,        only : exchange_t
use domain_mod,                   only : domain_t

implicit none

type, extends(dp_interpolator_t) :: nearest_neighbour_interpolator_t
    integer(kind=4), allocatable :: indx(:), indy(:), indz(:)
    integer(kind=4)              :: task_size = 0
    type(mesh_t),    pointer     :: source_mesh => null()

    contains
    procedure :: init => init_nn
    procedure :: calc_weights => calc_w_nn
    procedure :: interpolate => interpolate_nn
    procedure :: reallocate_if_required => reallocate_if_required_nn
end type

type, extends(dp_interpolator_t) :: bilinear_interpolator_t
    integer(kind=4),   allocatable :: indx(:), indy(:)
    real(kind=8),      allocatable :: zdx(:), zdy(:)
    class(exchange_t), allocatable :: exchange_halo
    type(mesh_t),      pointer     :: source_mesh => null()
    integer(kind=4)                :: task_size = 0

    contains
    procedure :: init => init_bilinear
    procedure :: calc_weights => calc_w_bilinear
    procedure :: interpolate => interpolate_bilinear
    procedure :: ext_halo => ext_halo_bilinear
    procedure :: reallocate_if_required => reallocate_if_required_bilinear
end type

type, extends(dp_interpolator_t) :: bicubic_Ah_interpolator_t
    integer(kind=4),   allocatable :: indx(:), indy(:)
    real(kind=8),      allocatable :: wx(:,:),wy(:,:)
    class(exchange_t), allocatable :: exchange_halo
    type(mesh_t),      pointer     :: source_mesh => null()
    integer(kind=4)                :: task_size = 0

    contains
    procedure :: init => init_bicubic_Ah
    procedure :: calc_weights => calc_w_bicubic_Ah
    procedure :: interpolate => interpolate_bicubic_Ah
    procedure :: ext_halo => ext_halo_bicubic_Ah
    procedure :: reallocate_if_required => reallocate_if_required_bicubic
end type

type, extends(dp_interpolator_t) :: trilinear_Ah_interpolator_t
    integer(kind=4),   allocatable :: indx(:), indy(:), indz(:)
    real(kind=8),      allocatable :: zdx(:), zdy(:), zdz(:)
    class(exchange_t), allocatable :: exchange_halo
    type(mesh_t),      pointer     :: source_mesh => null()
    integer(kind=4)                :: task_size = 0

    contains
    procedure :: init => init_trilinear
    procedure :: calc_weights => calc_w_trilinear
    procedure :: interpolate => interpolate_trilinear
    procedure :: ext_halo => ext_halo_trilinear
    procedure :: reallocate_if_required => reallocate_if_required_trilinear
end type

type, extends(dp_interpolator_t) :: tricubic_Ah_interpolator_t
    integer(kind=4),   allocatable :: indx(:), indy(:), indz(:)
    real(kind=8),      allocatable :: wx(:,:),wy(:,:), wz(:,:)
    class(exchange_t), allocatable :: exchange_halo
    type(mesh_t),      pointer     :: source_mesh
    integer(kind=4) :: task_size = 0

    contains
    procedure :: init => init_tricubic_Ah
    procedure :: calc_weights => calc_w_tricubic_Ah
    procedure :: interpolate => interpolate_tricubic_Ah
    procedure :: ext_halo => ext_halo_tricubic_Ah
    procedure :: reallocate_if_required => reallocate_if_required_tricubic_Ah
end type

contains

subroutine init_nn(this, mesh_name, domain)

    class(nearest_neighbour_interpolator_t), intent(inout) :: this
    character(len=*), intent(in) :: mesh_name
    type(domain_t),   intent(in) :: domain

    allocate(this%indx(this%min_alloc_size))
    allocate(this%indy(this%min_alloc_size))
    allocate(this%indz(this%min_alloc_size))

    call domain%get_mesh(this%source_mesh,mesh_name)

end subroutine

subroutine reallocate_if_required_nn(this, need_size)

    class(nearest_neighbour_interpolator_t), intent(inout) :: this
    integer(kind=4),                         intent(in)    :: need_size

    if(need_size <= size(this%indx,1)) return

    deallocate(this%indx, this%indy, this%indz)

    allocate(this%indx(need_size+this%extension_overhead))
    allocate(this%indy(need_size+this%extension_overhead))
    allocate(this%indz(need_size+this%extension_overhead))

end subroutine

subroutine calc_w_nn(this,alpha,beta,eta,tile_ind)

    class(nearest_neighbour_interpolator_t), intent(inout) :: this
    real(kind=8),     intent(in) :: alpha(:), beta(:), eta(:)
    integer(kind=4),  intent(in) :: tile_ind        

    integer(kind=4) :: i
    type(tile_mesh_t), pointer :: mesh

    mesh => this%source_mesh%tile(tile_ind)

    this%task_size = size(alpha,1)
    call this%reallocate_if_required(this%task_size)

    do i = 1, this%task_size
        this%indx(i) = int((alpha(i) - mesh%alpha_0) / mesh%hx +1 - mesh%shift_i)
        this%indy(i) = int((beta(i)  - mesh%beta_0)  / mesh%hy +1 - mesh%shift_j)
        this%indz(i) = int(eta(i)   / mesh%hz +1 - mesh%shift_k)
        this%indz(i) = max(1,min(mesh%nz,this%indz(i)))
    end do

end subroutine

subroutine interpolate_nn(this, qout, qin)

    class(nearest_neighbour_interpolator_t), intent(in) :: this
    real(kind=8),             intent(inout) :: qout(:)
    type(tile_field_t),       intent(in)    :: qin

    integer(kind=4) :: i

    do i = 1, this%task_size
        qout(i) = qin%p(this%indx(i),this%indy(i),this%indz(i))
    end do

end subroutine

subroutine init_bilinear(this,mesh_name,domain)

    class(bilinear_interpolator_t), intent(inout) :: this
    character(len=*), intent(in) :: mesh_name
    type(domain_t),   intent(in) :: domain

    allocate(this%indx(this%min_alloc_size))
    allocate(this%indy(this%min_alloc_size))
    allocate(this%zdx(this%min_alloc_size))
    allocate(this%zdy(this%min_alloc_size))

    call domain%get_mesh(this%source_mesh,mesh_name)

end subroutine

subroutine reallocate_if_required_bilinear(this,need_size)

    class(bilinear_interpolator_t), intent(inout) :: this
    integer(kind=4),                intent(in)    :: need_size

    if(need_size <= size(this%indx,1)) return

    deallocate(this%indx, this%indy, this%zdx, this%zdy)

    allocate(this%indx(need_size+this%extension_overhead))
    allocate(this%indy(need_size+this%extension_overhead))
    allocate(this%zdx (need_size+this%extension_overhead))
    allocate(this%zdy (need_size+this%extension_overhead))

end subroutine

subroutine calc_w_bilinear(this,alpha,beta,eta,tile_ind)

    class(bilinear_interpolator_t), intent(inout) :: this
    real(kind=8),     intent(in)    :: alpha(:), beta(:), eta(:)
    integer(kind=4),  intent(in)    :: tile_ind        

    integer(kind=4) :: i
    type(tile_mesh_t), pointer :: mesh

    mesh => this%source_mesh%tile(tile_ind)

    this%task_size = size(alpha,1)
    call this%reallocate_if_required(this%task_size)

    do i = 1, this%task_size
        this%zdx(i)  = (alpha(i) - mesh%alpha_0) / mesh%hx +1 - mesh%shift_i
        this%zdy(i)  = (beta(i)  - mesh%beta_0)  / mesh%hy +1 - mesh%shift_j
        this%indx(i) = min(int(this%zdx(i)),mesh%nx-1)
        this%indy(i) = min(int(this%zdy(i)),mesh%ny-1)
        this%zdx(i)  = this%zdx(i) - this%indx(i)
        this%zdy(i)  = this%zdy(i) - this%indy(i)
    end do

end subroutine

subroutine interpolate_bilinear(this, qout, qin)

    class(bilinear_interpolator_t), intent(in) :: this
    real(kind=8),             intent(inout) :: qout(:)
    type(tile_field_t),       intent(in)    :: qin

    integer(kind=4) :: i
    real(kind=8) :: p1, p2

    do i = 1, this%task_size
        p1 = qin%p(this%indx(i),this%indy(i)  ,1)+this%zdx(i)*(qin%p(this%indx(i)+1,this%indy(i)  ,1)-qin%p(this%indx(i),this%indy(i)  ,1))
        p2 = qin%p(this%indx(i),this%indy(i)+1,1)+this%zdx(i)*(qin%p(this%indx(i)+1,this%indy(i)+1,1)-qin%p(this%indx(i),this%indy(i)+1,1))
        qout(i) = p1 + this%zdy(i)*(p2-p1)
    end do

end subroutine

subroutine ext_halo_bilinear(this,q,domain)
    class(bilinear_interpolator_t), intent(inout) :: this
    type(grid_field_t),       intent(inout) :: q
    type(domain_t),           intent(in)    :: domain

    call this%exchange_halo%do(q,domain%parcomm)
end subroutine

subroutine init_bicubic_Ah(this, mesh_name, domain)

    class(bicubic_Ah_interpolator_t), intent(inout) :: this
    character(len=*), intent(in) :: mesh_name
    type(domain_t),   intent(in) :: domain

    allocate(this%indx(this%min_alloc_size))
    allocate(this%indy(this%min_alloc_size))
    allocate(this%wx(-1:2,this%min_alloc_size))
    allocate(this%wy(-1:2,this%min_alloc_size))

    call domain%get_mesh(this%source_mesh,mesh_name)

end subroutine

subroutine reallocate_if_required_bicubic(this,need_size)

    class(bicubic_Ah_interpolator_t), intent(inout) :: this
    integer(kind=4),                  intent(in)    :: need_size

    if(need_size <= size(this%indx,1)) return

    deallocate(this%indx, this%indy, this%wx, this%wy)

    allocate(this%indx(need_size+this%extension_overhead))
    allocate(this%indy(need_size+this%extension_overhead))
    allocate(this%wx (-1:2,need_size+this%extension_overhead))
    allocate(this%wy (-1:2,need_size+this%extension_overhead))

end subroutine

subroutine calc_w_bicubic_Ah(this,alpha,beta,eta,tile_ind)

    class(bicubic_Ah_interpolator_t), intent(inout) :: this
    real(kind=8),      intent(in)    :: alpha(:), beta(:), eta(:)
    integer(kind=4),   intent(in) :: tile_ind

    integer(kind=4) :: i
    type(tile_mesh_t), pointer :: mesh
    real(kind=8) :: p1, p2, p3, p4, zdx, zdy

    mesh => this%source_mesh%tile(tile_ind)

    this%task_size = size(alpha,1)
    call this%reallocate_if_required(this%task_size)

    do i = 1, this%task_size
        zdx  = (alpha(i) - mesh%alpha_0) / mesh%hx +1 - mesh%shift_i
        zdy  = (beta(i)  - mesh%beta_0)  / mesh%hy +1 - mesh%shift_j
        this%indx(i) = max(min(int(zdx),mesh%nx-2),2)
        this%indy(i) = max(min(int(zdy),mesh%ny-2),2)
        zdx  = zdx - this%indx(i)
        zdy  = zdy - this%indy(i)

        p1 = zdx+1.0_8
        p2 = zdx
        p3 = zdx-1.0_8
        p4 = zdx-2.0_8

        this%wx(-1,i) =-p2*p3*p4 / 6.0_8
        this%wx( 0,i) = p1*p3*p4 / 2.0_8
        this%wx( 1,i) =-p1*p2*p4 / 2.0_8
        this%wx( 2,i) = p1*p2*p3 / 6.0_8

        p1 = zdy+1.0_8
        p2 = zdy
        p3 = zdy-1.0_8
        p4 = zdy-2.0_8

        this%wy(-1,i) =-p2*p3*p4 / 6.0_8
        this%wy( 0,i) = p1*p3*p4 / 2.0_8
        this%wy( 1,i) =-p1*p2*p4 / 2.0_8
        this%wy( 2,i) = p1*p2*p3 / 6.0_8
    end do

end subroutine

subroutine interpolate_bicubic_Ah(this, qout, qin)

    class(bicubic_Ah_interpolator_t), intent(in) :: this
    real(kind=8),             intent(inout) :: qout(:)
    type(tile_field_t),       intent(in)    :: qin

    integer(kind=4) :: i, indx, indy, j
    real(kind=8) :: py1, py2, py3, py4

    do i = 1, this%task_size

        indx = this%indx(i)
        indy = this%indy(i)

        py1 = 0.0_8
        py2 = 0.0_8
        py3 = 0.0_8
        py4 = 0.0_8

        do j = -1,2
            py1 = py1 + this%wx(j,i)*qin%p(indx+j,indy-1,1)
            py2 = py2 + this%wx(j,i)*qin%p(indx+j,indy  ,1)
            py3 = py3 + this%wx(j,i)*qin%p(indx+j,indy+1,1)
            py4 = py4 + this%wx(j,i)*qin%p(indx+j,indy+2,1)
        end do

        qout(i) = this%wy(-1,i)*py1 + this%wy(0,i)*py2 + this%wy(1,i)*py3 + this%wy(2,i)*py4
    end do

end subroutine

subroutine ext_halo_bicubic_Ah(this,q,domain)
    class(bicubic_Ah_interpolator_t), intent(inout) :: this
    type(grid_field_t),       intent(inout) :: q
    type(domain_t),           intent(in)    :: domain

    call this%exchange_halo%do(q,domain%parcomm)
end subroutine

subroutine init_trilinear(this,mesh_name,domain)

    class(trilinear_Ah_interpolator_t), intent(inout) :: this
    character(len=*), intent(in) :: mesh_name
    type(domain_t),   intent(in) :: domain

    allocate(this%indx(this%min_alloc_size))
    allocate(this%indy(this%min_alloc_size))
    allocate(this%indz(this%min_alloc_size))
    allocate(this%zdx(this%min_alloc_size))
    allocate(this%zdy(this%min_alloc_size))
    allocate(this%zdz(this%min_alloc_size))

    call domain%get_mesh(this%source_mesh,mesh_name)

end subroutine

subroutine reallocate_if_required_trilinear(this,need_size)

    class(trilinear_Ah_interpolator_t), intent(inout) :: this
    integer(kind=4),                    intent(in)    :: need_size

    if(need_size <= size(this%indx,1)) return

    deallocate(this%indx, this%indy, this%indz, this%zdx, this%zdy, this%zdz)

    allocate(this%indx(need_size+this%extension_overhead))
    allocate(this%indy(need_size+this%extension_overhead))
    allocate(this%indz(need_size+this%extension_overhead))
    allocate(this%zdx (need_size+this%extension_overhead))
    allocate(this%zdy (need_size+this%extension_overhead))
    allocate(this%zdz (need_size+this%extension_overhead))

end subroutine

subroutine calc_w_trilinear(this,alpha,beta,eta,tile_ind)

    class(trilinear_Ah_interpolator_t), intent(inout) :: this
    real(kind=8),     intent(in)    :: alpha(:), beta(:), eta(:)
    integer(kind=4),  intent(in)    :: tile_ind        

    integer(kind=4) :: i
    type(tile_mesh_t), pointer :: mesh
    real(kind=8) :: eta_1, eta_N

    mesh => this%source_mesh%tile(tile_ind)
    eta_1 = mesh%get_eta(1)
    eta_N = mesh%get_eta(mesh%nz)

    this%task_size = size(alpha,1)
    call this%reallocate_if_required(this%task_size)

    do i = 1, this%task_size
        this%zdx(i)  = (alpha(i) - mesh%alpha_0) / mesh%hx +1 - mesh%shift_i
        this%zdy(i)  = (beta(i)  - mesh%beta_0)  / mesh%hy +1 - mesh%shift_j
        this%zdz(i)  = max(eta_1, min(eta_N, eta(i)) ) / mesh%hz +1 - mesh%shift_k
        this%indx(i) = min(int(this%zdx(i)),mesh%nx-1)
        this%indy(i) = min(int(this%zdy(i)),mesh%ny-1)
        this%indz(i) = min(int(this%zdz(i)),mesh%nz-1)
        this%zdx(i)  = this%zdx(i) - this%indx(i)
        this%zdy(i)  = this%zdy(i) - this%indy(i)
        this%zdz(i)  = this%zdz(i) - this%indz(i)
    end do

end subroutine

subroutine interpolate_trilinear(this, qout, qin)

    class(trilinear_Ah_interpolator_t), intent(in) :: this
    real(kind=8),             intent(inout) :: qout(:)
    type(tile_field_t),       intent(in)    :: qin

    integer(kind=4) :: i, ii, jj, kk
    real(kind=8) :: px1, px2, py1, py2

    do i = 1, this%task_size

        ii = this%indx(i)
        jj = this%indy(i)
        kk = this%indz(i)

        px1 = qin%p(ii,jj  ,kk) + this%zdx(i)*(qin%p(ii+1,jj  ,kk) - qin%p(ii,jj  ,kk))
        px2 = qin%p(ii,jj+1,kk) + this%zdx(i)*(qin%p(ii+1,jj+1,kk) - qin%p(ii,jj+1,kk))
        py1 = px1 + this%zdy(i)*(px2-px1)

        px1 = qin%p(ii,jj  ,kk+1) + this%zdx(i)*(qin%p(ii+1,jj  ,kk+1) - qin%p(ii,jj  ,kk+1))
        px2 = qin%p(ii,jj+1,kk+1) + this%zdx(i)*(qin%p(ii+1,jj+1,kk+1) - qin%p(ii,jj+1,kk+1))
        py2 = px1 + this%zdy(i)*(px2-px1)

        qout(i) = py1 + this%zdz(i)*(py2-py1)
    end do

end subroutine

subroutine ext_halo_trilinear(this,q,domain)
    class(trilinear_Ah_interpolator_t), intent(inout) :: this
    type(grid_field_t),       intent(inout) :: q
    type(domain_t),           intent(in)    :: domain

    call this%exchange_halo%do(q,domain%parcomm)
end subroutine

subroutine init_tricubic_Ah(this, mesh_name, domain)

    class(tricubic_Ah_interpolator_t), intent(inout) :: this
    character(len=*), intent(in) :: mesh_name
    type(domain_t),   intent(in) :: domain

    allocate(this%indx(this%min_alloc_size))
    allocate(this%indy(this%min_alloc_size))
    allocate(this%indz(this%min_alloc_size))
    allocate(this%wx(-1:2,this%min_alloc_size))
    allocate(this%wy(-1:2,this%min_alloc_size))
    allocate(this%wz(-1:2,this%min_alloc_size))

    call domain%get_mesh(this%source_mesh,mesh_name)

end subroutine

subroutine reallocate_if_required_tricubic_Ah(this,need_size)

    class(tricubic_Ah_interpolator_t), intent(inout) :: this
    integer(kind=4),                   intent(in)    :: need_size

    if(need_size <= size(this%indx,1)) return

    deallocate(this%indx, this%indy, this%indz, this%wx, this%wy, this%wz)

    allocate(this%indx(need_size+this%extension_overhead))
    allocate(this%indy(need_size+this%extension_overhead))
    allocate(this%indz(need_size+this%extension_overhead))
    allocate(this%wx(-1:2,need_size+this%extension_overhead))
    allocate(this%wy(-1:2,need_size+this%extension_overhead))
    allocate(this%wz(-1:2,need_size+this%extension_overhead))

end subroutine

subroutine calc_w_tricubic_Ah(this,alpha,beta,eta,tile_ind)

    class(tricubic_Ah_interpolator_t), intent(inout) :: this
    real(kind=8),    intent(in)    :: alpha(:), beta(:), eta(:)
    integer(kind=4), intent(in)    :: tile_ind

    type(tile_mesh_t), pointer :: mesh
    integer(kind=4) :: i
    real(kind=8) :: p1, p2, p3, p4, zdx, zdy, zdz, eta_1, eta_N

    mesh => this%source_mesh%tile(tile_ind)
    eta_1 = mesh%get_eta(1)
    eta_N = mesh%get_eta(mesh%nz)

    this%task_size = size(alpha,1)
    call this%reallocate_if_required(this%task_size)

    do i = 1, this%task_size
        zdx  = (alpha(i) - mesh%alpha_0) / mesh%hx +1 - mesh%shift_i
        zdy  = (beta(i)  - mesh%beta_0)  / mesh%hy +1 - mesh%shift_j
        zdz  = min(eta_N, max(eta_1, eta(i)) )  / mesh%hz +1 - mesh%shift_k

        this%indx(i) = max(min(int(zdx),mesh%nx-2),2)
        this%indy(i) = max(min(int(zdy),mesh%ny-2),2)
        this%indz(i) = max(min(int(zdz),mesh%nz-2),2)
        zdx  = zdx - this%indx(i)
        zdy  = zdy - this%indy(i)
        zdz  = zdz - this%indz(i)

        if(zdx < 0.0_8) then
            zdx = zdx+1.0_8
            this%wx(-1:2,i) = [1.0_8-zdx,zdx,0.0_8,0.0_8]
        else if(zdx > 1.0_8) then
            zdx = zdx-1.0_8
            this%wx(-1:2,i) = [0.0_8,0.0_8,1.0_8-zdx,zdx]
        else
        
            p1 = zdx+1.0_8
            p2 = zdx
            p3 = zdx-1.0_8
            p4 = zdx-2.0_8

            this%wx(-1,i) =-p2*p3*p4 / 6.0_8
            this%wx( 0,i) = p1*p3*p4 / 2.0_8
            this%wx( 1,i) =-p1*p2*p4 / 2.0_8
            this%wx( 2,i) = p1*p2*p3 / 6.0_8

        end if

        if(zdy < 0.0_8) then
            zdy = zdy+1.0_8
            this%wy(-1:2,i) = [1.0_8-zdy,zdy,0.0_8,0.0_8]
        else if(zdy > 1.0_8) then
            zdy = zdy-1.0_8
            this%wy(-1:2,i) = [0.0_8,0.0_8,1.0_8-zdy,zdy]
        else
            p1 = zdy+1.0_8
            p2 = zdy
            p3 = zdy-1.0_8
            p4 = zdy-2.0_8

            this%wy(-1,i) =-p2*p3*p4 / 6.0_8
            this%wy( 0,i) = p1*p3*p4 / 2.0_8
            this%wy( 1,i) =-p1*p2*p4 / 2.0_8
            this%wy( 2,i) = p1*p2*p3 / 6.0_8

        end if

        p1 = zdz+1.0_8
        p2 = zdz
        p3 = zdz-1.0_8
        p4 = zdz-2.0_8

        this%wz(-1,i) =-p2*p3*p4 / 6.0_8
        this%wz( 0,i) = p1*p3*p4 / 2.0_8
        this%wz( 1,i) =-p1*p2*p4 / 2.0_8
        this%wz( 2,i) = p1*p2*p3 / 6.0_8
    end do

end subroutine

subroutine interpolate_tricubic_Ah(this, qout, qin)

    class(tricubic_Ah_interpolator_t), intent(in) :: this
    real(kind=8),             intent(inout) :: qout(:)
    type(tile_field_t),       intent(in)    :: qin

    integer(kind=4) :: i, indx, indy, indz, k, j, l
    real(kind=8) :: py, pz

    do i = 1, this%task_size

        indx = this%indx(i)
        indy = this%indy(i)
        indz = this%indz(i)

        qout(i) = 0.0_8
        do k = -1,2
            pz = 0.0_8
            do j = -1,2
                py = 0.0_8
                do l = -1,2
                    py = py + this%wx(l,i)*qin%p(indx+l, indy+j, indz+k)
                end do
                pz = pz + this%wy(j,i)*py
            end do
            qout(i) = qout(i) + this%wz(k,i)*pz
        end do

    end do

end subroutine

subroutine ext_halo_tricubic_Ah(this,q,domain)
    class(tricubic_Ah_interpolator_t), intent(inout) :: this
    type(grid_field_t),       intent(inout) :: q
    type(domain_t),           intent(in)    :: domain

    call this%exchange_halo%do(q,domain%parcomm)
end subroutine

end module