module sbp_diff_63_mod

use sbp_diff_mod,   only : sbp_diff_t
use grid_field_mod, only : grid_field_t, tile_field_t
use mesh_mod,       only : tile_mesh_t, mesh_t
use tile_mod,       only : tile_t
use parcomm_mod,    only : parcomm_global

implicit none

real(kind=8), parameter :: h11 = 0.13649d5 / 0.43200d5, &
                           h22 = 0.12013d5 / 0.8640d4,  &
                           h33 = 0.2711d4  / 0.4320d4,  &
                           h44 = 0.5359d4  / 0.4320d4,  &
                           h55 = 0.7877d4  / 0.8640d4,  &
                           h66 = 0.43801d5 / 0.43200d5

real(kind=8), parameter :: c56 =5591070156686698065364559.d0/7931626489314500743872000.d0

real(kind=8), parameter :: c11 = -1.d0 / 2.d0
real(kind=8), parameter :: c12 = -0.953d3 / 0.16200d5 + c56
real(kind=8), parameter :: c13 = 0.715489d6 / 0.259200d6 - (4.d0 * c56)
real(kind=8), parameter :: c14 = -0.62639d5 / 0.14400d5 + (6.d0 * c56)
real(kind=8), parameter :: c15 = 0.147127d6 / 0.51840d5 - (4.d0 * c56)
real(kind=8), parameter :: c16 = -0.89387d5 / 0.129600d6 + c56
real(kind=8), parameter :: c23 = -0.57139d5 / 0.8640d4 + (10.d0 * c56)
real(kind=8), parameter :: c24 = 0.745733d6 / 0.51840d5 - (20.d0 * c56)
real(kind=8), parameter :: c25 = -0.18343d5 / 0.1728d4 + (15.d0 * c56)
real(kind=8), parameter :: c26 = 0.240569d6 / 0.86400d5 - (4.d0 * c56)
real(kind=8), parameter :: c34 = -0.176839d6 / 0.12960d5 + (20.d0 * c56)
real(kind=8), parameter :: c35 = 0.242111d6 / 0.17280d5 - (20.d0 * c56)
real(kind=8), parameter :: c36 = -0.182261d6 / 0.43200d5 + (6.d0 * c56)
real(kind=8), parameter :: c45 = -0.165041d6 / 0.25920d5 + (10.d0 * c56)
real(kind=8), parameter :: c46 = 0.710473d6 / 0.259200d6 - (4.d0 * c56)
real(kind=8), parameter :: c47 = 1.d0 / 6.d1
real(kind=8), parameter :: c57 = -3.D0 / 2.d1
real(kind=8), parameter :: c58 = 1.d0 / 6.d1
real(kind=8), parameter :: c67 = 3.d0 / 4.d0
real(kind=8), parameter :: c68 = -3.d0 / 2.d1
real(kind=8), parameter :: c69 = 1.d0 / 6.d1

real(kind=8), parameter :: q11 =  c11/h11, q12 =  c12/h11, q13 =  c13/h11, &
                           q14 =  c14/h11, q15 =  c15/h11, q16 =  c16/h11, &
                           q21 = -c12/h22, q23 =  c23/h22, q24 =  c24/h22, &
                           q25 =  c25/h22, q26 =  c26/h22,                 &
                           q31 = -c13/h33, q32 = -c23/h33, q34 =  c34/h33, &
                           q35 =  c35/h33, q36 =  c36/h33,                 &
                           q41 = -c14/h44, q42 = -c24/h44, q43 = -c34/h44, &
                           q45 =  c45/h44, q46 =  c46/h44, q47 =  c47/h44, &
                           q51 = -c15/h55, q52 = -c25/h55, q53 = -c35/h55, &
                           q54 = -c45/h55, q56 =  c56/h55, q57 =  c57/h55, &
                           q58 =  c58/h55,                                 &
                           q61 = -c16/h66, q62 = -c26/h66, q63 = -c36/h66, &
                           q64 = -c46/h66, q65 = -c56/h66, q67 =  c67/h66, &
                           q68 =  c68/h66, q69 =  c69/h66
type, public, extends(sbp_diff_t) :: sbp_diff_63_t
contains
    procedure, public :: apply_tile
    procedure, public :: add_SAT_correction
end type sbp_diff_63_t

contains

subroutine apply_tile(this, f_out, f_in, mesh, direction)
    class(sbp_diff_63_t), intent(in)    :: this
    type(tile_field_t),      intent(inout) :: f_out
    type(tile_field_t),      intent(in)    :: f_in
    type(tile_mesh_t),       intent(in)    :: mesh
    character(len=*),        intent(in)    :: direction

    integer(kind=4) :: k, i, j

    select case(direction)
    case("x")
        do k = mesh%ks, mesh%ke
            do j = mesh%js, mesh%je
                if (mesh%is<=1 .and. mesh%ie>=1) then
                    f_out%p(1,j,k) = q11*f_in%p(1,j,k)+q12*f_in%p(2,j,k)+q13*f_in%p(3,j,k) + &
                                   + q14*f_in%p(4,j,k)+q15*f_in%p(5,j,k)+q16*f_in%p(6,j,k)
                end if
                if (mesh%is<=2 .and. mesh%ie>=2) then
                    f_out%p(2,j,k) = q21*f_in%p(1,j,k)                  +q23*f_in%p(3,j,k) + &
                                   + q24*f_in%p(4,j,k)+q25*f_in%p(5,j,k)+q26*f_in%p(6,j,k)
                end if
                if (mesh%is<=3 .and. mesh%ie>=3) then
                    f_out%p(3,j,k) = q31*f_in%p(1,j,k)+q32*f_in%p(2,j,k)                   + &
                                   + q34*f_in%p(4,j,k)+q35*f_in%p(5,j,k)+q36*f_in%p(6,j,k)
                end if
                if (mesh%is<=4 .and. mesh%ie>=4) then
                    f_out%p(4,j,k) = q41*f_in%p(1,j,k)+q42*f_in%p(2,j,k)+q43*f_in%p(3,j,k) + &
                                                      +q45*f_in%p(5,j,k)+q46*f_in%p(6,j,k) + &
                                   + q47*f_in%p(7,j,k)
                end if
                if (mesh%is<=5 .and. mesh%ie>=5) then
                    f_out%p(5,j,k) = q51*f_in%p(1,j,k)+q52*f_in%p(2,j,k)+q53*f_in%p(3,j,k) + &
                                   + q54*f_in%p(4,j,k)                  +q56*f_in%p(6,j,k) + &
                                   + q57*f_in%p(7,j,k)+q58*f_in%p(8,j,k)
                end if
                if (mesh%is<=6 .and. mesh%ie>=6) then
                    f_out%p(6,j,k) = q61*f_in%p(1,j,k)+q62*f_in%p(2,j,k)+q63*f_in%p(3,j,k) + &
                                   + q64*f_in%p(4,j,k)+q65*f_in%p(5,j,k)                   + &
                                   + q67*f_in%p(7,j,k)+q68*f_in%p(8,j,k)+q69*f_in%p(9,j,k)
                end if
                do i = max(mesh%is, 7), min(mesh%ie, mesh%nx-6)
                    f_out%p(i,j,k) = (    -f_in%p(i-3,j,k) + &
                                         9*f_in%p(i-2,j,k) + &
                                       -45*f_in%p(i-1,j,k) + &
                                        45*f_in%p(i+1,j,k) + &
                                        -9*f_in%p(i+2,j,k) + &
                                           f_in%p(i+3,j,k)) / 60.0_8
                end do
                if (mesh%is<=mesh%nx-5 .and. mesh%ie>=mesh%nx-5) then
                    f_out%p(mesh%nx-5,j,k) = -q61*f_in%p(mesh%nx  ,j,k)-q62*f_in%p(mesh%nx-1,j,k)-q63*f_in%p(mesh%nx-2,j,k) &
                                             -q64*f_in%p(mesh%nx-3,j,k)-q65*f_in%p(mesh%nx-4,j,k)                           &
                                             -q67*f_in%p(mesh%nx-6,j,k)-q68*f_in%p(mesh%nx-7,j,k)-q69*f_in%p(mesh%nx-8,j,k)
                end if
                if (mesh%is<=mesh%nx-4 .and. mesh%ie>=mesh%nx-4) then
                    f_out%p(mesh%nx-4,j,k) = -q51*f_in%p(mesh%nx  ,j,k)-q52*f_in%p(mesh%nx-1,j,k)-q53*f_in%p(mesh%nx-2,j,k) &
                                             -q54*f_in%p(mesh%nx-3,j,k)                          -q56*f_in%p(mesh%nx-5,j,k) &
                                             -q57*f_in%p(mesh%nx-6,j,k)-q58*f_in%p(mesh%nx-7,j,k)
                end if
                if (mesh%is<=mesh%nx-3 .and. mesh%ie>=mesh%nx-3) then
                    f_out%p(mesh%nx-3,j,k) = -q41*f_in%p(mesh%nx  ,j,k)-q42*f_in%p(mesh%nx-1,j,k)-q43*f_in%p(mesh%nx-2,j,k) &
                                                                       -q45*f_in%p(mesh%nx-4,j,k)-q46*f_in%p(mesh%nx-5,j,k) &
                                             -q47*f_in%p(mesh%nx-6,j,k)
                end if
                if (mesh%is<=mesh%nx-2 .and. mesh%ie>=mesh%nx-2) then
                    f_out%p(mesh%nx-2,j,k) = -q31*f_in%p(mesh%nx  ,j,k)-q32*f_in%p(mesh%nx-1,j,k)                           &
                                             -q34*f_in%p(mesh%nx-3,j,k)-q35*f_in%p(mesh%nx-4,j,k)-q36*f_in%p(mesh%nx-5,j,k)
                end if
                if (mesh%is<=mesh%nx-1 .and. mesh%ie>=mesh%nx-1) then
                    f_out%p(mesh%nx-1,j,k) = -q21*f_in%p(mesh%nx  ,j,k)                          -q23*f_in%p(mesh%nx-2,j,k) &
                                             -q24*f_in%p(mesh%nx-3,j,k)-q25*f_in%p(mesh%nx-4,j,k)-q26*f_in%p(mesh%nx-5,j,k)
                end if
                if (mesh%is<=mesh%nx .and. mesh%ie>=mesh%nx) then
                    f_out%p(mesh%nx  ,j,k) = -q11*f_in%p(mesh%nx  ,j,k)-q12*f_in%p(mesh%nx-1,j,k)-q13*f_in%p(mesh%nx-2,j,k) &
                                             -q14*f_in%p(mesh%nx-3,j,k)-q15*f_in%p(mesh%nx-4,j,k)-q16*f_in%p(mesh%nx-5,j,k)
                end if
            end do
        end do
    case("y")
        do k = mesh%ks, mesh%ke
            if (mesh%js<=1 .and. mesh%je>=1) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,1,k) = q11*f_in%p(i,1,k)+q12*f_in%p(i,2,k)+q13*f_in%p(i,3,k) + &
                                   + q14*f_in%p(i,4,k)+q15*f_in%p(i,5,k)+q16*f_in%p(i,6,k)
                end do
            end if
            if (mesh%js<=2 .and. mesh%je>=2) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,2,k) = q21*f_in%p(i,1,k)                  +q23*f_in%p(i,3,k) + &
                                   + q24*f_in%p(i,4,k)+q25*f_in%p(i,5,k)+q26*f_in%p(i,6,k)
                end do
            end if
            if (mesh%js<=3 .and. mesh%je>=3) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,3,k) = q31*f_in%p(i,1,k)+q32*f_in%p(i,2,k)                   + &
                                   + q34*f_in%p(i,4,k)+q35*f_in%p(i,5,k)+q36*f_in%p(i,6,k)
                end do
            end if
            if (mesh%js<=4 .and. mesh%je>=4) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,4,k) = q41*f_in%p(i,1,k)+q42*f_in%p(i,2,k)+q43*f_in%p(i,3,k) + &
                                                      +q45*f_in%p(i,5,k)+q46*f_in%p(i,6,k) + &
                                   + q47*f_in%p(i,7,k)
                end do
            end if
            if (mesh%js<=5 .and. mesh%je>=5) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,5,k) = q51*f_in%p(i,1,k)+q52*f_in%p(i,2,k)+q53*f_in%p(i,3,k) + &
                                   + q54*f_in%p(i,4,k)                  +q56*f_in%p(i,6,k) + &
                                   + q57*f_in%p(i,7,k)+q58*f_in%p(i,8,k)
                end do
            end if
            if (mesh%js<=6 .and. mesh%je>=6) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,6,k) = q61*f_in%p(i,1,k)+q62*f_in%p(i,2,k)+q63*f_in%p(i,3,k) + &
                                   + q64*f_in%p(i,4,k)+q65*f_in%p(i,5,k)                   + &
                                   + q67*f_in%p(i,7,k)+q68*f_in%p(i,8,k)+q69*f_in%p(i,9,k)
                end do
            end if
            do j = max(mesh%js, 7), min(mesh%je, mesh%ny-6)
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,k) = (    -f_in%p(i,j-3,k) + &
                                         9*f_in%p(i,j-2,k) + &
                                       -45*f_in%p(i,j-1,k) + &
                                        45*f_in%p(i,j+1,k) + &
                                        -9*f_in%p(i,j+2,k) + &
                                           f_in%p(i,j+3,k)) / 60.0_8
                end do
            end do
            if (mesh%js<=mesh%ny-5 .and. mesh%je>=mesh%ny-5) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-5,k) = -q61*f_in%p(i,mesh%ny  ,k)-q62*f_in%p(i,mesh%ny-1,k)-q63*f_in%p(i,mesh%ny-2,k) &
                                             -q64*f_in%p(i,mesh%ny-3,k)-q65*f_in%p(i,mesh%ny-4,k)                           &
                                             -q67*f_in%p(i,mesh%ny-6,k)-q68*f_in%p(i,mesh%ny-7,k)-q69*f_in%p(i,mesh%ny-8,k)
                end do
            end if
            if (mesh%js<=mesh%ny-4 .and. mesh%je>=mesh%ny-4) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-4,k) = -q51*f_in%p(i,mesh%ny  ,k)-q52*f_in%p(i,mesh%ny-1,k)-q53*f_in%p(i,mesh%ny-2,k) &
                                             -q54*f_in%p(i,mesh%ny-3,k)                          -q56*f_in%p(i,mesh%ny-5,k) &
                                             -q57*f_in%p(i,mesh%ny-6,k)-q58*f_in%p(i,mesh%ny-7,k)
                end do
            end if
            if (mesh%js<=mesh%ny-3 .and. mesh%je>=mesh%ny-3) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-3,k) = -q41*f_in%p(i,mesh%ny  ,k)-q42*f_in%p(i,mesh%ny-1,k)-q43*f_in%p(i,mesh%ny-2,k) &
                                                                       -q45*f_in%p(i,mesh%ny-4,k)-q46*f_in%p(i,mesh%ny-5,k) &
                                             -q47*f_in%p(i,mesh%ny-6,k)
                end do
            end if
            if (mesh%js<=mesh%ny-2 .and. mesh%je>=mesh%ny-2) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-2,k) = -q31*f_in%p(i,mesh%ny  ,k)-q32*f_in%p(i,mesh%ny-1,k)                           &
                                             -q34*f_in%p(i,mesh%ny-3,k)-q35*f_in%p(i,mesh%ny-4,k)-q36*f_in%p(i,mesh%ny-5,k)
                end do
            end if
            if (mesh%js<=mesh%ny-1 .and. mesh%je>=mesh%ny-1) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-1,k) = -q21*f_in%p(i,mesh%ny  ,k)                          -q23*f_in%p(i,mesh%ny-2,k) &
                                             -q24*f_in%p(i,mesh%ny-3,k)-q25*f_in%p(i,mesh%ny-4,k)-q26*f_in%p(i,mesh%ny-5,k)
                end do
            end if
            if (mesh%js<=mesh%ny .and. mesh%je>=mesh%ny) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny  ,k) = -q11*f_in%p(i,mesh%ny  ,k)-q12*f_in%p(i,mesh%ny-1,k)-q13*f_in%p(i,mesh%ny-2,k) &
                                             -q14*f_in%p(i,mesh%ny-3,k)-q15*f_in%p(i,mesh%ny-4,k)-q16*f_in%p(i,mesh%ny-5,k)
                end do
            end if
        end do
    case("z")
        if (mesh%ks<=1 .and. mesh%ke>=1) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,1) = q11*f_in%p(i,j,1)+q12*f_in%p(i,j,2)+q13*f_in%p(i,j,3) + &
                                   + q14*f_in%p(i,j,4)+q15*f_in%p(i,j,5)+q16*f_in%p(i,j,6)
                end do
            end do
        end if
        if (mesh%ks<=2 .and. mesh%ke>=2) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,2) = q21*f_in%p(i,j,1)                  +q23*f_in%p(i,j,3) + &
                                   + q24*f_in%p(i,j,4)+q25*f_in%p(i,j,5)+q26*f_in%p(i,j,6)
                end do
            end do
        end if
        if (mesh%ks<=3 .and. mesh%ke>=3) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,3) = q31*f_in%p(i,j,1)+q32*f_in%p(i,j,2)                   + &
                                   + q34*f_in%p(i,j,4)+q35*f_in%p(i,j,5)+q36*f_in%p(i,j,6)
                end do
            end do
        end if
        if (mesh%ks<=4 .and. mesh%ke>=4) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,4) = q41*f_in%p(i,j,1)+q42*f_in%p(i,j,2)+q43*f_in%p(i,j,3) + &
                                                      +q45*f_in%p(i,j,5)+q46*f_in%p(i,j,6) + &
                                   + q47*f_in%p(i,j,7)
                end do
            end do
        end if
        if (mesh%ks<=5 .and. mesh%ke>=5) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,5) = q51*f_in%p(i,j,1)+q52*f_in%p(i,j,2)+q53*f_in%p(i,j,3) + &
                                   + q54*f_in%p(i,j,4)                  +q56*f_in%p(i,j,6) + &
                                   + q57*f_in%p(i,j,7)+q58*f_in%p(i,j,8)
                end do
            end do
        end if
        if (mesh%ks<=6 .and. mesh%ke>=6) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,6) = q61*f_in%p(i,j,1)+q62*f_in%p(i,j,2)+q63*f_in%p(i,j,3) + &
                                   + q64*f_in%p(i,j,4)+q65*f_in%p(i,j,5)                   + &
                                   + q67*f_in%p(i,j,7)+q68*f_in%p(i,j,8)+q69*f_in%p(i,j,9)
                end do
            end do
        end if
        do k = max(mesh%ks, 7), min(mesh%ke, mesh%nz-6)
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,k) = (    -f_in%p(i,j,k-3) + &
                                         9*f_in%p(i,j,k-2) + &
                                       -45*f_in%p(i,j,k-1) + &
                                        45*f_in%p(i,j,k+1) + &
                                        -9*f_in%p(i,j,k+2) + &
                                           f_in%p(i,j,k+3)) / 60.0_8
                end do
            end do
        end do
        if (mesh%ks<=mesh%nz-5 .and. mesh%ke>=mesh%nz-5) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-5) = -q61*f_in%p(i,j,mesh%nz  )-q62*f_in%p(i,j,mesh%nz-1)-q63*f_in%p(i,j,mesh%nz-2) &
                                             -q64*f_in%p(i,j,mesh%nz-3)-q65*f_in%p(i,j,mesh%nz-4)                           &
                                             -q67*f_in%p(i,j,mesh%nz-6)-q68*f_in%p(i,j,mesh%nz-7)-q69*f_in%p(i,j,mesh%nz-8)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz-4 .and. mesh%ke>=mesh%nz-4) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-4) = -q51*f_in%p(i,j,mesh%nz  )-q52*f_in%p(i,j,mesh%nz-1)-q53*f_in%p(i,j,mesh%nz-2) &
                                             -q54*f_in%p(i,j,mesh%nz-3)                          -q56*f_in%p(i,j,mesh%nz-5) &
                                             -q57*f_in%p(i,j,mesh%nz-6)-q58*f_in%p(i,j,mesh%nz-7)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz-3 .and. mesh%ke>=mesh%nz-3) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-3) = -q41*f_in%p(i,j,mesh%nz  )-q42*f_in%p(i,j,mesh%nz-1)-q43*f_in%p(i,j,mesh%nz-2) &
                                                                       -q45*f_in%p(i,j,mesh%nz-4)-q46*f_in%p(i,j,mesh%nz-5) &
                                             -q47*f_in%p(i,j,mesh%nz-6)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz-2 .and. mesh%ke>=mesh%nz-2) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-2) = -q31*f_in%p(i,j,mesh%nz  )-q32*f_in%p(i,j,mesh%nz-1)                           &
                                             -q34*f_in%p(i,j,mesh%nz-3)-q35*f_in%p(i,j,mesh%nz-4)-q36*f_in%p(i,j,mesh%nz-5)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz-1 .and. mesh%ke>=mesh%nz-1) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-1) = -q21*f_in%p(i,j,mesh%nz  )                          -q23*f_in%p(i,j,mesh%nz-2) &
                                             -q24*f_in%p(i,j,mesh%nz-3)-q25*f_in%p(i,j,mesh%nz-4)-q26*f_in%p(i,j,mesh%nz-5)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz .and. mesh%ke>=mesh%nz) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz  ) = -q11*f_in%p(i,j,mesh%nz  )-q12*f_in%p(i,j,mesh%nz-1)-q13*f_in%p(i,j,mesh%nz-2) &
                                             -q14*f_in%p(i,j,mesh%nz-3)-q15*f_in%p(i,j,mesh%nz-4)-q16*f_in%p(i,j,mesh%nz-5)
                end do
            end do
        end if
    case default
        call parcomm_global%abort("Wrong direction in sbp_diff_63_mod")
    end select

end subroutine apply_tile

subroutine add_SAT_correction(this, f_out, f_in, mesh, direction)
    class(sbp_diff_63_t), intent(in)    :: this
    type(tile_field_t),   intent(inout) :: f_out
    type(tile_field_t),   intent(in)    :: f_in
    type(tile_mesh_t),    intent(in)    :: mesh
    character(len=*),     intent(in)    :: direction

    call parcomm_global%abort("add_SAT_correction not implemented for sbp_diff_63_t")
end subroutine add_SAT_correction

end module sbp_diff_63_mod
