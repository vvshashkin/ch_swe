module particle_values_mod

use string_mod,    only : string_t
use particles_mod, only : particles_t

implicit none

type, public :: particle_values_t

    type(string_t), allocatable :: field_names(:)
    integer(kind=4) :: Nfld, N
    real(kind=8), allocatable :: p(:,:)

    contains
        procedure :: init

end type particle_values_t

contains

subroutine init(this,N,field_names)

    class(particle_values_t), intent(out) :: this
    integer(kind=4),          intent(in)  :: N
    type(string_t),           intent(in)  :: field_names(:)

    this%field_names = field_names
    this%Nfld = size(field_names,1)
    this%N = N

    allocate(this%p(N,this%Nfld))

end subroutine

end module particle_values_mod