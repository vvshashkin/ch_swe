module sbp_diff_c2i_63_mod

use sbp_diff_mod,   only : sbp_diff_t
use grid_field_mod, only : grid_field_t, tile_field_t
use mesh_mod,       only : tile_mesh_t, mesh_t
use tile_mod,       only : tile_t
use parcomm_mod,    only : parcomm_global

use sbp_diff_i2c_63_mod, only : c34, c55

implicit none

! dot prodict matrix coeffs
real(kind=8), parameter :: h11 =  95.0_8 / 288.0_8, &
                           h22 =  317._8 / 240._8,  &
                           h33 =  23.0_8 /  30.0_8,  &
                           h44 = 793.0_8 / 720.0_8, &
                           h55 = 157.0_8 / 160.0_8

! boundary projection operator coeffs
real(kind=8), parameter :: R11 =  35.0_8 / 16.0_8, &
                           R12 = -35.0_8 / 16.0_8, &
                           R13 =  21.0_8 / 16.0_8, &
                           R14 = - 5.0_8 / 16.0_8

real(kind=8), parameter :: q11 = -60711983.0_8/21888000.0_8+312623.0_8/456000.0_8*c55+8639.0_8/36480.0_8*c34, &
                           q12 = 101173243.0_8/17510400.0_8-312623.0_8/182400.0_8*c55-8639.0_8/ 9728.0_8*c34, &
                           q13 = -2593653.0_8/486400.0_8+8639.0_8/7296.0_8*c34, &
                           q14 = -8639.0_8/14592.0_8*c34+7121893.0_8/1751040.0_8+312623.0_8/91200.0_8*c55, &
                           q15 = -312623.0_8/91200.0_8*c55-5209847.0_8/2188800.0_8, &
                           q16 = 18712829.0_8/29184000.0_8+312623.0_8/304000.0_8*c55+8639.0_8/145920.0_8*c34, &
                           q21 = -53376169.0_8/43822080.0_8-312623.0_8/456480.0_8*c55-43195.0_8/182592.0_8*c34, &
                           q22 =  7190801.0_8/5842944.0_8+312623.0_8/182592.0_8*c55+215975.0_8/243456.0_8*c34, & 
                           q23 = -215975.0_8/182592.0_8*c34+2846555.0_8/2191104.0_8, &
                           q24 = -27181195.0_8/8764416.0_8-312623.0_8/91296.0_8*c55+215975.0_8/365184.0_8*c34, & 
                           q25 =  2407853.0_8/973824.0_8+312623.0_8/91296.0_8*c55, & 
                           q26 = -59866697.0_8/87644160.0_8-312623.0_8/304320.0_8*c55-43195.0_8/730368.0_8*c34, &
                           q31 = 1807.0_8/1920.0_8+312623.0_8/176640.0_8*c55+43195.0_8/70656.0_8*c34, &
                           q32 = -6326795.0_8/2260992.0_8-312623.0_8/70656.0_8*c55-215975.0_8/94208.0_8*c34, &
                           q33 =  215975.0_8/70656.0_8*c34-221735.0_8/188416.0_8, &
                           q34 =  8940511.0_8/1130496.0_8+312623.0_8/35328.0_8*c55-215975.0_8/141312.0_8*c34, &
                           q35 =  -3843253.0_8/565248.0_8-312623.0_8/35328.0_8*c55, &
                           q36 =  7252803.0_8/3768320.0_8+312623.0_8/117760.0_8*c55+43195.0_8/282624.0_8*c34, &
                           q41 = -288299.0_8/599040.0_8-312623.0_8/380640.0_8*c55-43195.0_8/152256.0_8*c34, &
                           q42 = 230891.0_8/239616.0_8+312623.0_8/152256.0_8*c55+215975.0_8/203008.0_8*c34, &
                           q43 = -215975.0_8/152256.0_8*c34, &
                           q44 =  215975.0_8/304512.0_8*c34-312623.0_8/76128.0_8*c55-355691.0_8/119808.0_8, &
                           q45 =  418091.0_8/119808.0_8+312623.0_8/76128.0_8*c55, &
                           q46 =  -400619.0_8/399360.0_8-312623.0_8/253760.0_8*c55-43195.0_8/609024.0_8*c34, &
                           q51 = 9606527.0_8/65111040.0_8+312623.0_8/1356480.0_8*c55+43195.0_8/542592.0_8*c34, &
                           q52 = -4598783.0_8/17362944.0_8-312623.0_8/542592.0_8*c55-215975.0_8/723456.0_8*c34, &
                           q53 = -4811905.0_8/13022208.0_8+215975.0_8/542592.0_8*c34, &
                           q54 = 5665537.0_8/26044416.0_8-215975.0_8/1085184.0_8*c34+312623.0_8/271296.0_8*c55, &
                           q55 =  -312623.0_8/271296.0_8*c55, &
                           q56 =  68894207.0_8/260444160.0_8+312623.0_8/904320.0_8*c55+43195.0_8/2170368.0_8*c34, &
                           q57 =  3.0_8/628.0_8

real(kind=8), parameter :: qc1 =   75.0_8 /  64.0_8, &
                           qc2 = - 25.0_8 / 384.0_8, &
                           qc3 =    3.0_8 / 640.0_8

type, public, extends(sbp_diff_t) :: sbp_diff_c2i_63_t
contains
    procedure, public :: apply_tile
    procedure, public :: add_SAT_correction
end type sbp_diff_c2i_63_t

contains

subroutine apply_tile(this, f_out, f_in, mesh, direction)
    class(sbp_diff_c2i_63_t), intent(in)    :: this
    type(tile_field_t),       intent(inout) :: f_out
    type(tile_field_t),       intent(in)    :: f_in
    type(tile_mesh_t),        intent(in)    :: mesh
    character(len=*),         intent(in)    :: direction

    integer(kind=4) :: k, i, j

    select case(direction)
    case("x")
        do k = mesh%ks, mesh%ke
            do j = mesh%js, mesh%je
                if (mesh%is<=1 .and. mesh%ie>=1) then
                    f_out%p(1,j,k) = q11*f_in%p(1,j,k)+q12*f_in%p(2,j,k)+q13*f_in%p(3,j,k)+q14*f_in%p(4,j,k)+q15*f_in%p(5,j,k)+q16*f_in%p(6,j,k)
                end if
                if (mesh%is<=2 .and. mesh%ie>=2) then
                    f_out%p(2,j,k) = q21*f_in%p(1,j,k)+q22*f_in%p(2,j,k)+q23*f_in%p(3,j,k)+q24*f_in%p(4,j,k)+q25*f_in%p(5,j,k)+q26*f_in%p(6,j,k)
                end if
                if (mesh%is<=3 .and. mesh%ie>=3) then
                    f_out%p(3,j,k) = q31*f_in%p(1,j,k)+q32*f_in%p(2,j,k)+q33*f_in%p(3,j,k)+q34*f_in%p(4,j,k)+q35*f_in%p(5,j,k)+q36*f_in%p(6,j,k)
                end if
                if (mesh%is<=4 .and. mesh%ie>=4) then
                    f_out%p(4,j,k) = q41*f_in%p(1,j,k)+q42*f_in%p(2,j,k)+q43*f_in%p(3,j,k)+q44*f_in%p(4,j,k)+q45*f_in%p(5,j,k)+q46*f_in%p(6,j,k)
                end if
                if (mesh%is<=5 .and. mesh%ie>=5) then
                    f_out%p(5,j,k) = q51*f_in%p(1,j,k)+q52*f_in%p(2,j,k)+q53*f_in%p(3,j,k)+q54*f_in%p(4,j,k)+q55*f_in%p(5,j,k)+q56*f_in%p(6,j,k)+q57*f_in%p(7,j,k)
                end if
                do i = max(mesh%is, 6), min(mesh%ie, mesh%nx-5)
                    f_out%p(i,j,k) = qc1*(f_in%p(i  ,j,k) - f_in%p(i-1,j,k)) &
                                   + qc2*(f_in%p(i+1,j,k) - f_in%p(i-2,j,k)) &
                                   + qc3*(f_in%p(i+2,j,k) - f_in%p(i-3,j,k))
                end do
                if (mesh%is<=mesh%nx-4 .and. mesh%ie>=mesh%nx-4) then
                    f_out%p(mesh%nx-4,j,k) = -q51*f_in%p(mesh%nx-1,j,k)-q52*f_in%p(mesh%nx-2,j,k)-q53*f_in%p(mesh%nx-3,j,k)-q54*f_in%p(mesh%nx-4,j,k)-q55*f_in%p(mesh%nx-5,j,k)-q56*f_in%p(mesh%nx-6,j,k)-q57*f_in%p(mesh%nx-7,j,k)
                end if
                if (mesh%is<=mesh%nx-3 .and. mesh%ie>=mesh%nx-3) then
                    f_out%p(mesh%nx-3,j,k) = -q41*f_in%p(mesh%nx-1,j,k)-q42*f_in%p(mesh%nx-2,j,k)-q43*f_in%p(mesh%nx-3,j,k)-q44*f_in%p(mesh%nx-4,j,k)-q45*f_in%p(mesh%nx-5,j,k)-q46*f_in%p(mesh%nx-6,j,k)
                end if
                if (mesh%is<=mesh%nx-2 .and. mesh%ie>=mesh%nx-2) then
                    f_out%p(mesh%nx-2,j,k) = -q31*f_in%p(mesh%nx-1,j,k)-q32*f_in%p(mesh%nx-2,j,k)-q33*f_in%p(mesh%nx-3,j,k)-q34*f_in%p(mesh%nx-4,j,k)-q35*f_in%p(mesh%nx-5,j,k)-q36*f_in%p(mesh%nx-6,j,k)
                end if
                if (mesh%is<=mesh%nx-1 .and. mesh%ie>=mesh%nx-1) then
                    f_out%p(mesh%nx-1,j,k) = -q21*f_in%p(mesh%nx-1,j,k)-q22*f_in%p(mesh%nx-2,j,k)-q23*f_in%p(mesh%nx-3,j,k)-q24*f_in%p(mesh%nx-4,j,k)-q25*f_in%p(mesh%nx-5,j,k)-q26*f_in%p(mesh%nx-6,j,k)
                end if
                if (mesh%is<=mesh%nx .and. mesh%ie>=mesh%nx) then
                    f_out%p(mesh%nx,j,k) = -q11*f_in%p(mesh%nx-1,j,k)-q12*f_in%p(mesh%nx-2,j,k)-q13*f_in%p(mesh%nx-3,j,k)-q14*f_in%p(mesh%nx-4,j,k)-q15*f_in%p(mesh%nx-5,j,k)-q16*f_in%p(mesh%nx-6,j,k)
                end if
            end do
        end do
    case("y")
        do k = mesh%ks, mesh%ke
            if (mesh%js<=1 .and. mesh%je>=1) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,1,k) = q11*f_in%p(i,1,k)+q12*f_in%p(i,2,k)+q13*f_in%p(i,3,k)+q14*f_in%p(i,4,k)+q15*f_in%p(i,5,k)+q16*f_in%p(i,6,k)
                end do
            end if
            if (mesh%js<=2 .and. mesh%je>=2) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,2,k) = q21*f_in%p(i,1,k)+q22*f_in%p(i,2,k)+q23*f_in%p(i,3,k)+q24*f_in%p(i,4,k)+q25*f_in%p(i,5,k)+q26*f_in%p(i,6,k)
                end do
            end if
            if (mesh%js<=3 .and. mesh%je>=3) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,3,k) = q31*f_in%p(i,1,k)+q32*f_in%p(i,2,k)+q33*f_in%p(i,3,k)+q34*f_in%p(i,4,k)+q35*f_in%p(i,5,k)+q36*f_in%p(i,6,k)
                end do
            end if
            if (mesh%js<=4 .and. mesh%je>=4) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,4,k) = q41*f_in%p(i,1,k)+q42*f_in%p(i,2,k)+q43*f_in%p(i,3,k)+q44*f_in%p(i,4,k)+q45*f_in%p(i,5,k)+q46*f_in%p(i,6,k)
                end do
            end if
            if (mesh%js<=5 .and. mesh%je>=5) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,5,k) = q51*f_in%p(i,1,k)+q52*f_in%p(i,2,k)+q53*f_in%p(i,3,k)+q54*f_in%p(i,4,k)+q55*f_in%p(i,5,k)+q56*f_in%p(i,6,k)+q57*f_in%p(i,7,k)
                end do
            end if
            do j = max(mesh%js, 6), min(mesh%je, mesh%ny-5)
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,k) = qc1*(f_in%p(i,j  ,k) - f_in%p(i,j-1,k)) &
                                   + qc2*(f_in%p(i,j+1,k) - f_in%p(i,j-2,k)) &
                                   + qc3*(f_in%p(i,j+2,k) - f_in%p(i,j-3,k))
                end do
            end do
            if (mesh%js<=mesh%ny-4 .and. mesh%je>=mesh%ny-4) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-4,k) = -q51*f_in%p(i,mesh%ny-1,k)-q52*f_in%p(i,mesh%ny-2,k)-q53*f_in%p(i,mesh%ny-3,k)-q54*f_in%p(i,mesh%ny-4,k)-q55*f_in%p(i,mesh%ny-5,k)-q56*f_in%p(i,mesh%ny-6,k)-q57*f_in%p(i,mesh%ny-7,k)
                end do
            end if
            if (mesh%js<=mesh%ny-3 .and. mesh%je>=mesh%ny-3) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-3,k) = -q41*f_in%p(i,mesh%ny-1,k)-q42*f_in%p(i,mesh%ny-2,k)-q43*f_in%p(i,mesh%ny-3,k)-q44*f_in%p(i,mesh%ny-4,k)-q45*f_in%p(i,mesh%ny-5,k)-q46*f_in%p(i,mesh%ny-6,k)
                end do
            end if
            if (mesh%js<=mesh%ny-2 .and. mesh%je>=mesh%ny-2) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-2,k) = -q31*f_in%p(i,mesh%ny-1,k)-q32*f_in%p(i,mesh%ny-2,k)-q33*f_in%p(i,mesh%ny-3,k)-q34*f_in%p(i,mesh%ny-4,k)-q35*f_in%p(i,mesh%ny-5,k)-q36*f_in%p(i,mesh%ny-6,k)
                end do
            end if
            if (mesh%js<=mesh%ny-1 .and. mesh%je>=mesh%ny-1) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-1,k) = -q21*f_in%p(i,mesh%ny-1,k)-q22*f_in%p(i,mesh%ny-2,k)-q23*f_in%p(i,mesh%ny-3,k)-q24*f_in%p(i,mesh%ny-4,k)-q25*f_in%p(i,mesh%ny-5,k)-q26*f_in%p(i,mesh%ny-6,k)
                end do
            end if
            if (mesh%js<=mesh%ny .and. mesh%je>=mesh%ny) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny,k) = -q11*f_in%p(i,mesh%ny-1,k)-q12*f_in%p(i,mesh%ny-2,k)-q13*f_in%p(i,mesh%ny-3,k)-q14*f_in%p(i,mesh%ny-4,k)-q15*f_in%p(i,mesh%ny-5,k)-q16*f_in%p(i,mesh%ny-6,k)
                end do
            end if
        end do
    case("z")
        if (mesh%ks<=1 .and. mesh%ke>=1) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,1) = q11*f_in%p(i,j,1)+q12*f_in%p(i,j,2)+q13*f_in%p(i,j,3)+q14*f_in%p(i,j,4)+q15*f_in%p(i,j,5)+q16*f_in%p(i,j,6)
                end do
            end do
        end if
        if (mesh%ks<=2 .and. mesh%ke>=2) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,2) = q21*f_in%p(i,j,1)+q22*f_in%p(i,j,2)+q23*f_in%p(i,j,3)+q24*f_in%p(i,j,4)+q25*f_in%p(i,j,5)+q26*f_in%p(i,j,6)
                end do
            end do
        end if
        if (mesh%ks<=3 .and. mesh%ke>=3) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,3) = q31*f_in%p(i,j,1)+q32*f_in%p(i,j,2)+q33*f_in%p(i,j,3)+q34*f_in%p(i,j,4)+q35*f_in%p(i,j,5)+q36*f_in%p(i,j,6)
                end do
            end do
        end if
        if (mesh%ks<=4 .and. mesh%ke>=4) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,4) = q41*f_in%p(i,j,1)+q42*f_in%p(i,j,2)+q43*f_in%p(i,j,3)+q44*f_in%p(i,j,4)+q45*f_in%p(i,j,5)+q46*f_in%p(i,j,6)
                end do
            end do
        end if
        if (mesh%ks<=5 .and. mesh%ke>=5) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,5) = q51*f_in%p(i,j,1)+q52*f_in%p(i,j,2)+q53*f_in%p(i,j,3)+q54*f_in%p(i,j,4)+q55*f_in%p(i,j,5)+q56*f_in%p(i,j,6)+q57*f_in%p(i,j,7)
                end do
            end do
        end if
        do k = max(mesh%ks, 6), min(mesh%ke, mesh%nz-5)
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,k) = qc1*(f_in%p(i,j,k  ) - f_in%p(i,j,k-1)) &
                                   + qc2*(f_in%p(i,j,k+1) - f_in%p(i,j,k-2)) &
                                   + qc3*(f_in%p(i,j,k+2) - f_in%p(i,j,k-3))
                end do
            end do
        end do
        if (mesh%ks<=mesh%nz-4 .and. mesh%ke>=mesh%nz-4) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-4) = -q51*f_in%p(i,j,mesh%nz-1)-q52*f_in%p(i,j,mesh%nz-2)-q53*f_in%p(i,j,mesh%nz-3)-q54*f_in%p(i,j,mesh%nz-4)-q55*f_in%p(i,j,mesh%nz-5)-q56*f_in%p(i,j,mesh%nz-6)-q57*f_in%p(i,j,mesh%nz-7)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz-3 .and. mesh%ke>=mesh%nz-3) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-3) = -q41*f_in%p(i,j,mesh%nz-1)-q42*f_in%p(i,j,mesh%nz-2)-q43*f_in%p(i,j,mesh%nz-3)-q44*f_in%p(i,j,mesh%nz-4)-q45*f_in%p(i,j,mesh%nz-5)-q46*f_in%p(i,j,mesh%nz-6)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz-2 .and. mesh%ke>=mesh%nz-2) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-2) = -q31*f_in%p(i,j,mesh%nz-1)-q32*f_in%p(i,j,mesh%nz-2)-q33*f_in%p(i,j,mesh%nz-3)-q34*f_in%p(i,j,mesh%nz-4)-q35*f_in%p(i,j,mesh%nz-5)-q36*f_in%p(i,j,mesh%nz-6)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz-1 .and. mesh%ke>=mesh%nz-1) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-1) = -q21*f_in%p(i,j,mesh%nz-1)-q22*f_in%p(i,j,mesh%nz-2)-q23*f_in%p(i,j,mesh%nz-3)-q24*f_in%p(i,j,mesh%nz-4)-q25*f_in%p(i,j,mesh%nz-5)-q26*f_in%p(i,j,mesh%nz-6)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz .and. mesh%ke>=mesh%nz) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz) = -q11*f_in%p(i,j,mesh%nz-1)-q12*f_in%p(i,j,mesh%nz-2)-q13*f_in%p(i,j,mesh%nz-3)-q14*f_in%p(i,j,mesh%nz-4)-q15*f_in%p(i,j,mesh%nz-5)-q16*f_in%p(i,j,mesh%nz-6)
                end do
            end do
        end if
    case default
        call parcomm_global%abort("Wrong direction in sbp_diff_c2i_63_mod::apply_tile")
    end select
                
end subroutine apply_tile


subroutine add_SAT_correction(this, f_out, f_in, mesh, direction)
    class(sbp_diff_c2i_63_t), intent(in)    :: this
    type(tile_field_t),       intent(inout) :: f_out
    type(tile_field_t),       intent(in)    :: f_in
    type(tile_mesh_t),        intent(in)    :: mesh
    character(len=*),         intent(in)    :: direction

    integer(kind=4) :: k, i, j

    select case(direction)
    case('x')
        if (mesh%is<=1 .and. mesh%ie>=1) then
            do k = mesh%ks, mesh%ke
                do j = mesh%js, mesh%je
                    f_out%p(1,j,k) = f_out%p(1,j,k) &
                                   + 0.5_8/h11*(R11*(f_in%p(1,j,k)-f_in%p( 0,j,k)) + &
                                                R12*(f_in%p(2,j,k)-f_in%p(-1,j,k)) + &
                                                R13*(f_in%p(3,j,k)-f_in%p(-2,j,k)) + &
                                                R14*(f_in%p(4,j,k)-f_in%p(-3,j,k)) )
                end do
            end do
        end if
        if (mesh%is<=mesh%nx .and. mesh%ie>=mesh%nx) then
            do k = mesh%ks, mesh%ke
                do j = mesh%js, mesh%je
                    f_out%p(mesh%nx,j,k) = f_out%p(mesh%nx,j,k) &
                                          + 0.5_8/h11*(R11*(f_in%p(mesh%nx  ,j,k)-f_in%p(mesh%nx-1,j,k)) + &
                                                       R12*(f_in%p(mesh%nx+1,j,k)-f_in%p(mesh%nx-2,j,k)) + &
                                                       R13*(f_in%p(mesh%nx+2,j,k)-f_in%p(mesh%nx-3,j,k)) + &
                                                       R14*(f_in%p(mesh%nx+3,j,k)-f_in%p(mesh%nx-4,j,k)) )
                end do
            end do
        end if
    case('y')
        if (mesh%js<=1 .and. mesh%je>=1) then
            do k = mesh%ks, mesh%ke
                do i = mesh%is, mesh%ie
                    f_out%p(i,1,k) = f_out%p(i,1,k) &
                                   + 0.5_8/h11*(R11*(f_in%p(i,1,k)-f_in%p(i, 0,k)) + &
                                                R12*(f_in%p(i,2,k)-f_in%p(i,-1,k)) + &
                                                R13*(f_in%p(i,3,k)-f_in%p(i,-2,k)) + &
                                                R14*(f_in%p(i,4,k)-f_in%p(i,-3,k)) )
                end do
            end do
        end if
        if (mesh%js<=mesh%ny .and. mesh%je>=mesh%ny) then
            do k = mesh%ks, mesh%ke
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny,k) = f_out%p(i,mesh%ny,k) &
                                          + 0.5_8/h11*(R11*(f_in%p(i,mesh%ny  ,k)-f_in%p(i,mesh%ny-1,k)) + &
                                                       R12*(f_in%p(i,mesh%ny+1,k)-f_in%p(i,mesh%ny-2,k)) + &
                                                       R13*(f_in%p(i,mesh%ny+2,k)-f_in%p(i,mesh%ny-3,k)) + &
                                                       R14*(f_in%p(i,mesh%ny+3,k)-f_in%p(i,mesh%ny-4,k)) )
                end do
            end do
        end if
    case default
        call parcomm_global%abort("Wrong direction in sbp_diff_c2i_63_mod::add_SAT_correction")
    end select

end subroutine add_SAT_correction

end module sbp_diff_c2i_63_mod