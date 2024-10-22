module particle_interp_factory_mod

use abstract_particle_interp_mod, only : particle_interp_t
use bilinear_interp_mod,          only : bilinear_interp_t
use bicubic_interp_mod,           only : bicubic_interp_t
use domain_mod,                   only : domain_t

implicit none

contains

subroutine create_particle_interp(particle_interp, particle_interp_name, &
                                  initial_task_size, mesh_name, domain)

    class(particle_interp_t), intent(out), allocatable :: particle_interp
    character(len=*), intent(in) :: particle_interp_name
    integer(kind=4),  intent(in) :: initial_task_size
    character(len=*), intent(in) :: mesh_name
    type(domain_t),   intent(in) :: domain

    select case (particle_interp_name)
    case("bicubic")
        call create_bicubic_particle_interp(particle_interp, initial_task_size, &
                                            mesh_name, domain)
    case("bilinear")
        call create_bilinear_particle_interp(particle_interp, initial_task_size, &
                                              mesh_name, domain)
    case default
        call domain%parcomm%abort("particle_interp_factory_mod, create_particle_interp unknown interp name: "// particle_interp_name)
    end select

end subroutine

subroutine create_bicubic_particle_interp(particle_interp, initial_task_size, &
                                          mesh_name, domain)

    class(particle_interp_t), intent(out), allocatable :: particle_interp
    integer(kind=4),  intent(in) :: initial_task_size
    character(len=*), intent(in) :: mesh_name
    type(domain_t),   intent(in) :: domain

    type(bicubic_interp_t), allocatable :: bicubic_interp

    allocate(bicubic_interp)

    allocate(bicubic_interp%indx(initial_task_size))
    allocate(bicubic_interp%indy(initial_task_size))
    allocate(bicubic_interp%wx(-1:2,initial_task_size))
    allocate(bicubic_interp%wy(-1:2,initial_task_size))
    bicubic_interp%current_size = initial_task_size
    bicubic_interp%task_size = 0

    bicubic_interp%mesh_name = mesh_name

    call move_alloc(bicubic_interp, particle_interp)

end subroutine

subroutine create_bilinear_particle_interp(particle_interp, initial_task_size, &
                                           mesh_name, domain)

    class(particle_interp_t), intent(out), allocatable :: particle_interp
    integer(kind=4),  intent(in) :: initial_task_size
    character(len=*), intent(in) :: mesh_name
    type(domain_t),   intent(in) :: domain

    type(bilinear_interp_t), allocatable :: bilinear_interp

    allocate(bilinear_interp)

    allocate(bilinear_interp%indx(initial_task_size))
    allocate(bilinear_interp%indy(initial_task_size))
    allocate(bilinear_interp%dx(initial_task_size))
    allocate(bilinear_interp%dy(initial_task_size))
    bilinear_interp%current_size = initial_task_size
    bilinear_interp%task_size = 0

    bilinear_interp%mesh_name = mesh_name

    call move_alloc(bilinear_interp, particle_interp)

end subroutine

end module particle_interp_factory_mod