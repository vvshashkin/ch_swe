module bilinear_interp_mod

use abstract_particle_interp_mod, only : particle_interp_t
use domain_mod,                   only : domain_t
use particles_mod,                only : particles_t
use particle_values_mod,          only : particle_values_t
use grid_field_mod,               only : tile_field_t
use mesh_mod,                     only : mesh_t, tile_mesh_t
use array_tools_mod,              only : derealloc

type, extends(particle_interp_t) :: bilinear_interp_t

    integer(kind=4), allocatable :: indx(:), indy(:)
    real(kind=8),    allocatable :: dx(:), dy(:)
    integer(kind=4) :: current_size = 0, task_size = 0

    character(len=:), allocatable :: mesh_name

    contains

    procedure :: calc_w
    procedure :: interp_2d_levs

end type

contains

subroutine calc_w(this, particles, t, domain)

    class(bilinear_interp_t), intent(inout) :: this
    type(particles_t),        intent(in)    :: particles
    integer(kind=4),          intent(in)    :: t
    type(domain_t),           intent(in)    :: domain

    integer(kind=4) :: i
    integer(kind=4) :: indx, indy
    real(kind=8)    :: dx, dy

    type(mesh_t), pointer :: mesh
    type(tile_mesh_t), pointer :: tmesh

    call domain%get_mesh(mesh, this%mesh_name)
    tmesh => mesh%tile(t)

    if(particles%N > this%current_size) then
        call derealloc(this%indx,particles%N)
        call derealloc(this%indy,particles%N)
        call derealloc(this%dx,particles%N)
        call derealloc(this%dy,particles%N)
        this%current_size = particles%N
    end if

    this%task_size = particles%N

    do i = 1, particles%N

        dx   = (particles%coords%a(1,i) - tmesh%alpha_0) / tmesh%hx +1 - tmesh%shift_i
        dy   = (particles%coords%a(2,i)  - tmesh%beta_0)  / tmesh%hy +1 - tmesh%shift_j

        this%indx(i) = min(int(dx),tmesh%nx-1)
        this%indy(i) = min(int(dy),tmesh%ny-1)
        this%dx(i)   = dx - this%indx(i)
        this%dy(i)   = dy - this%indy(i)

    end do

end subroutine

subroutine interp_2d_levs(this, f_interp, fin, domain)

    class(bilinear_interp_t), intent(in)    :: this
    type(particle_values_t),  intent(inout) :: f_interp
    type(tile_field_t),       intent(in)    :: fin
    type(domain_t),           intent(in)    :: domain

    integer(kind=4) :: i, k
    real(kind=8)    :: p1, p2


    associate(ix => this%indx, iy => this%indy, dx => this%dx, dy => this%dy)

        do k = 1, f_interp%Nfld
            do i = 1, this%task_size
            
                p1 = fin%p(ix(i),iy(i)  ,k)+&
                      dx(i)*(fin%p(ix(i)+1,iy(i)  ,k)-fin%p(ix(i),iy(i)  ,k))
                p2 = fin%p(ix(i),iy(i)+1,k)+&
                      dx(i)*(fin%p(ix(i)+1,iy(i)+1,k)-fin%p(ix(i),iy(i)+1,k))
                f_interp%p(i,k) = p1 + dy(i)*(p2-p1)

            end do
        end do

    end associate

end subroutine

end module