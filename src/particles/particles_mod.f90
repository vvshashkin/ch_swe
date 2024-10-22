module particles_mod

use smart_array_mod, only : vec_array_r8_t, smart_array_i4_t
use string_mod,      only : string_t

implicit none

type, public :: particles_t

    integer(kind=4) :: n_coords, N
    type(string_t), allocatable :: coord_names(:)

    type(vec_array_r8_t)   :: coords
    type(smart_array_i4_t) :: face_ind

    contains

    procedure, public :: init
    procedure, public :: shrink_to_actuall_size
    procedure, public :: reset
    procedure, public :: add_particle
    procedure, public :: add_particles
    procedure, public :: add_particles_particles
    generic :: add => add_particle, add_particles, add_particles_particles

end type

contains

subroutine init(this, coord_names, size_alloc)

    class(particles_t), intent(out) :: this
    type(string_t),     intent(in)  :: coord_names(:)
    integer(kind=4),    intent(in)  :: size_alloc

    integer(kind=4) :: i

    this%N = 0

    this%coord_names = coord_names
    this%n_coords = size(coord_names,1)

    call this%coords%init(this%n_coords, size_alloc)
    
    call this%face_ind%init(size_alloc)
    
end subroutine

subroutine shrink_to_actuall_size(this)

    class(particles_t), intent(inout) :: this

    integer(kind=4) :: i

    call this%coords%shrink()
    call this%face_ind%shrink()

end subroutine

subroutine reset(this)

    class(particles_t), intent(inout) :: this

    integer(kind=4) :: i

    call this%coords%reset()
    call this%face_ind%reset()

    this%N = 0

end subroutine

subroutine add_particle(this, coords, face_ind, spare_size)

    class(particles_t), intent(inout) :: this
    real(kind=8),       intent(in)    :: coords(:)
    integer(kind=4),    intent(in)    :: face_ind

    integer(kind=4),    intent(in), optional :: spare_size

    integer(kind=4) :: spare_size_loc

    call this%coords%append(coords, spare_size)
    call this%face_ind%append(face_ind, spare_size)

    this%N = this%N+1

end subroutine

subroutine add_particles(this, coords, face_ind, spare_size)

    class(particles_t), intent(inout) :: this
    real(kind=8),       intent(in)    :: coords(:,:)
    integer(kind=4),    intent(in)    :: face_ind(:)

    integer(kind=4),    intent(in), optional :: spare_size

    integer(kind=4) :: spare_size_loc

    call this%coords%append(coords, spare_size)
    call this%face_ind%append(face_ind, spare_size)

    this%N = this%N+size(face_ind,1)

end subroutine

subroutine add_particles_particles(this, particles, spare_size)

    class(particles_t), intent(inout) :: this
    class(particles_t), intent(in)    :: particles

    integer(kind=4),    intent(in), optional :: spare_size

    integer(kind=4) :: spare_size_loc, n_coords, N

    N = particles%n
    n_coords = particles%n_coords

    call this%coords%append(particles%coords%a(1:n_coords,1:N), spare_size)
    call this%face_ind%append(particles%face_ind%a(1:N), spare_size)

    this%N = this%N+particles%N

end subroutine

end module