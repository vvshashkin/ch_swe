module bicubic_interp_mod

use abstract_particle_interp_mod, only : particle_interp_t
use domain_mod,                   only : domain_t
use particles_mod,                only : particles_t
use particle_values_mod,          only : particle_values_t
use grid_field_mod,               only : tile_field_t
use mesh_mod,                     only : mesh_t, tile_mesh_t
use array_tools_mod,              only : derealloc

type, extends(particle_interp_t) :: bicubic_interp_t

    integer(kind=4), allocatable :: indx(:), indy(:)
    real(kind=8),    allocatable :: wx(:,:), wy(:,:)
    integer(kind=4) :: current_size = 0, task_size = 0

    character(len=:), allocatable :: mesh_name

    contains

    procedure :: calc_w
    procedure :: interp_2d_levs

end type

contains

subroutine calc_w(this, particles, t, domain)

    class(bicubic_interp_t),  intent(inout) :: this
    type(particles_t),        intent(in)    :: particles
    integer(kind=4),          intent(in)    :: t
    type(domain_t),           intent(in)    :: domain

    integer(kind=4) :: i
    integer(kind=4) :: indx, indy
    real(kind=8)    :: dx, dy

    type(mesh_t), pointer :: mesh
    type(tile_mesh_t), pointer :: tmesh
    real(kind=8) :: p1, p2, p3, p4

    call domain%get_mesh(mesh, this%mesh_name)
    tmesh => mesh%tile(t)

    if(particles%N > this%current_size) then
        call derealloc(this%indx,particles%N)
        call derealloc(this%indy,particles%N)

        deallocate(this%wx)
        deallocate(this%wy)
        allocate(this%wx(-1:2,particles%N))
        allocate(this%wy(-1:2,particles%N))

        this%current_size = particles%N
    end if

    this%task_size = particles%N

    do i = 1, particles%N

        dx   = (particles%coords%a(1,i) - tmesh%alpha_0) / tmesh%hx +1 - tmesh%shift_i
        dy   = (particles%coords%a(2,i)  - tmesh%beta_0)  / tmesh%hy +1 - tmesh%shift_j

        this%indx(i) = min(max(2,int(dx)),tmesh%nx-2)
        this%indy(i) = min(max(2,int(dy)),tmesh%ny-2)
        dx   = dx - this%indx(i)
        dy   = dy - this%indy(i)

        p1 = dx+1.0_8
        p2 = dx
        p3 = dx-1.0_8
        p4 = dx-2.0_8

        this%wx(-1,i) =-p2*p3*p4 / 6.0_8
        this%wx( 0,i) = p1*p3*p4 / 2.0_8
        this%wx( 1,i) =-p1*p2*p4 / 2.0_8
        this%wx( 2,i) = p1*p2*p3 / 6.0_8

        p1 = dy+1.0_8
        p2 = dy
        p3 = dy-1.0_8
        p4 = dy-2.0_8

        this%wy(-1,i) =-p2*p3*p4 / 6.0_8
        this%wy( 0,i) = p1*p3*p4 / 2.0_8
        this%wy( 1,i) =-p1*p2*p4 / 2.0_8
        this%wy( 2,i) = p1*p2*p3 / 6.0_8

    end do

end subroutine

subroutine interp_2d_levs(this, f_interp, fin, domain)

    class(bicubic_interp_t),  intent(in)    :: this
    type(particle_values_t),  intent(inout) :: f_interp
    type(tile_field_t),       intent(in)    :: fin
    type(domain_t),           intent(in)    :: domain

    integer(kind=4) :: i, k, i1, j1
    real(kind=8)    :: p


    do k = 1, f_interp%Nfld
        do i = 1, this%task_size
            
            f_interp%p(i,k) = 0.0_8

            do j1 = -1, 2
                p = 0.0_8
                do i1 = -1,2
                    p = p + this%wx(i1,i)*fin%p(this%indx(i)+i1,this%indy(i)+j1,k)
                end do
                f_interp%p(i,k) = f_interp%p(i,k) + this%wy(j1,i)*p
            end do

        end do
    end do

end subroutine

end module