module latlon_grid_generator_mod

use const_mod, only : pi

implicit none

type, public :: latlon_grid_generator_t

    integer(kind=4) :: Nlat=0, Nlon=0
    real(kind=8), allocatable :: phi(:), lam(:)
    real(kind=8)              :: rotation_matrix(3,3) = reshape([1._8,0._8,0._8, &
                                                                 0._8,1._8,0._8, &
                                                                 0._8,0._8,1._8],[3,3])
    logical :: is_rotated = .false.

    contains
        procedure :: init
        procedure :: get_cartesian_coords
        procedure :: get_basis_vectors

end type

contains

subroutine init(this, Nlat, Nlon, rotation_matrix)

    class(latlon_grid_generator_t), intent(out) :: this
    integer(kind=4),                intent(in)  :: Nlat, Nlon
    real(kind=8),        optional,  intent(in)  :: rotation_matrix


    integer(kind=4) :: i, j

    this%Nlat = Nlat
    this%Nlon = Nlon
    
    allocate(this%lam(Nlon))
    do i = 1, Nlon
        this%lam(i) = (i-1)*2.0_8*pi / Nlon
    end do

    allocate(this%phi(Nlat))
    do j = 1, Nlat
        this%phi(j) = -0.5_8*pi+(j-1)*pi / (Nlat-1.0_8)
    end do

    this%is_rotated = .false.
    if(present(rotation_matrix)) then
        this%is_rotated = .true.
        this%rotation_matrix = rotation_matrix
    end if

end subroutine

pure subroutine get_cartesian_coords(this, r, i, j)

    class(latlon_grid_generator_t), intent(in) :: this

    real(kind=8),    intent(out) :: r(3)
    integer(kind=4), intent(in)  :: i, j

    real(kind=8) :: x, y, z

    r(1) = cos(this%phi(j))*cos(this%lam(i))
    r(2) = cos(this%phi(j))*sin(this%lam(i))
    r(3) = sin(this%phi(j))

    if(this%is_rotated) then

        x = sum(this%rotation_matrix(1,1:3)*r(1:3))
        y = sum(this%rotation_matrix(2,1:3)*r(1:3))
        z = sum(this%rotation_matrix(3,1:3)*r(1:3))

        r = [x,y,z]
    end if

end subroutine

pure subroutine get_basis_vectors(this, iv, jv, i, j)

    class(latlon_grid_generator_t), intent(in) :: this

    real(kind=8),    intent(out) :: iv(3), jv(3)
    integer(kind=4), intent(in)  :: i, j

    real(kind=8) :: x,y,z

    iv(1:3) = [-sin(this%lam(i)),cos(this%lam(i)),0.0_8]
    jv(1:3) = [-sin(this%phi(j))*cos(this%lam(i)),-sin(this%phi(j))*sin(this%lam(i)),cos(this%phi(j))]

    if(this%is_rotated) then

        x = sum(this%rotation_matrix(1,1:3)*iv(1:3))
        y = sum(this%rotation_matrix(2,1:3)*iv(1:3))
        z = sum(this%rotation_matrix(3,1:3)*iv(1:3))

        iv(1:3) = [x,y,z]

        x = sum(this%rotation_matrix(1,1:3)*jv(1:3))
        y = sum(this%rotation_matrix(2,1:3)*jv(1:3))
        z = sum(this%rotation_matrix(3,1:3)*jv(1:3))

        jv(1:3) = [x,y,z]

    end if

end subroutine

end module latlon_grid_generator_mod