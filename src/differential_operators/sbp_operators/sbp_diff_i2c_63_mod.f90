module sbp_diff_i2c_63_mod

use sbp_diff_mod,   only : sbp_diff_t
use grid_field_mod, only : grid_field_t, tile_field_t
use mesh_mod,       only : tile_mesh_t, mesh_t
use tile_mod,       only : tile_t
use parcomm_mod,    only : parcomm_global


implicit none
! dot prodict matrix coeffs
real(kind=8), parameter :: h11 = 325363.0_8 / 276480.0_8, &
                           h22 = 144001.0_8 / 276480.0_8, &
                           h33 =  43195.0_8 /  27648.0_8, &
                           h44 =  86857.0_8 / 138240.0_8, &
                           h55 = 312623.0_8 / 276480.0_8, &
                           h66 = 271229.0_8 / 276480.0_8

! boundary projection operator coeffs
real(kind=8), parameter :: R11 =  35.0_8 / 16.0_8, &
                           R12 = -35.0_8 / 16.0_8, &
                           R13 =  21.0_8 / 16.0_8, &
                           R14 = - 5.0_8 / 16.0_8
                           
! c34 and c55 are free parameters left after imposing summation-by-pats property.
! Current values obtained using minimization of ep4^T*Hp*ep4 + eu4^T*Hu*eu4
! where ep4 and eu4 are the fourth order polynomial differentiation error vectors
! i.e. e4p = Du2p * xu^4 - 4 * xp^3, e4u = Dp2u * xp^4 - 4 * xu^3 
real(kind=8), parameter :: c34 =  0.72279320467934028075_8 ! 805495023252648197134120997227394580304849/1114419751096024376224265917430133834276000
real(kind=8), parameter :: c55 = -0.77959568679678007297_8 ! -1006064115833286443955180421621841266681577/1290494717803051014490968617311302621777024
! real(kind=8), parameter :: c34 = 0.491150983242187_8, c55 = -0.724044733065893_8                           

real(kind=8), parameter :: q11 = -84440017.0_8/78087120.0_8-312623.0_8/1626815.0_8*c55-43195.0_8/650726.0_8*c34, &
                           q12 = 53376169.0_8/39043560.0_8+1250492.0_8/1626815.0_8*c55+86390.0_8/325363.0_8*c34, &
                           q13 = -997464.0_8/1626815.0_8-1875738.0_8/1626815.0_8*c55-129585.0_8/325363.0_8*c34, &
                           q14 = 17586239.0_8/39043560.0_8+1250492.0_8/1626815.0_8*c55+86390.0_8/325363.0_8*c34, &
                           q15 = -9606527.0_8/78087120.0_8-312623.0_8/1626815.0_8*c55-43195.0_8/650726.0_8*c34, &
                           q21 = 14948357.0_8/27648192.0_8+312623.0_8/288002.0_8*c55+12225.0_8/21736.0_8*c34, &
                           q22 = -7190801.0_8/2304016.0_8-625246.0_8/144001.0_8*c55-12225.0_8/5434.0_8*c34, &
                           q23 = 18980385.0_8/4608032.0_8+937869.0_8/144001.0_8*c55+36675.0_8/10868.0_8*c34, &
                           q24 = -14084351.0_8/6912048.0_8-625246.0_8/144001.0_8*c55-12225.0_8/5434.0_8*c34, &
                           q25 = 4598783.0_8/9216064.0_8+312623.0_8/288002.0_8*c55+12225.0_8/21736.0_8*c34, &
                           q31 = 1974879.0_8/6911200.0_8-1.0_8/4.0_8*c34, &
                           q32 = c34-569311.0_8/518340.0_8, &
                           q33 = -3.0_8/2.0_8*c34+399123.0_8/691120.0_8, &
                           q34 = c34, &
                           q35 = 962381.0_8/4146720.0_8-1.0_8/4.0_8*c34, &
                           q36 = 648.0_8/215975.0_8, &
                           q41 = 215975.0_8/694856.0_8*c34-27315065.0_8/16676544.0_8-312623.0_8/173714.0_8*c55, &
                           q42 = 27181195.0_8/4169136.0_8+625246.0_8/86857.0_8*c55-215975.0_8/173714.0_8*c34, &
                           q43 = -26821533.0_8/2779424.0_8-937869.0_8/86857.0_8*c55+647925.0_8/347428.0_8*c34, &
                           q44 = -215975.0_8/173714.0_8*c34+625246.0_8/86857.0_8*c55+21697151.0_8/4169136.0_8, &
                           q45 = -5665537.0_8/16676544.0_8+215975.0_8/694856.0_8*c34-312623.0_8/173714.0_8*c55, &
                           q46 = -9000.0_8/86857.0_8, &
                           q47 = 648.0_8/86857.0_8, &
                           q51 = c55+5209847.0_8/7502952.0_8, &
                           q52 = -7223559.0_8/2500984.0_8-4.0_8*c55, &
                           q53 = 11529759.0_8/2500984.0_8+6.0_8*c55, &
                           q54 = -25503551.0_8/7502952.0_8-4.0_8*c55, &
                           q55 = c55, &
                           q56 = 324000.0_8/312623.0_8, &
                           q57 = -18000.0_8/312623.0_8, &
                           q58 = 1296.0_8/312623.0_8, &
                           q61 = -18712829.0_8/86793280.0_8-937869.0_8/2712290.0_8*c55-43195.0_8/2169832.0_8*c34, &
                           q62 = 59866697.0_8/65094960.0_8+1875738.0_8/1356145.0_8*c55+43195.0_8/542458.0_8*c34, &
                           q63 = -65275227.0_8/43396640.0_8-2813607.0_8/1356145.0_8*c55-129585.0_8/1084916.0_8*c34, &
                           q64 = 24437759.0_8/21698320.0_8+1875738.0_8/1356145.0_8*c55+43195.0_8/542458.0_8*c34, &
                           q65 = -68894207.0_8/260379840.0_8-937869.0_8/2712290.0_8*c55-43195.0_8/2169832.0_8*c34, &
                           q66 = -324000.0_8/271229.0_8, &
                           q67 = 324000.0_8/271229.0_8, &
                           q68 = -18000.0_8/271229.0_8, &
                           q69 = 1296.0_8/271229.0_8

real(kind=8), parameter :: qc1 =   75.0_8 /  64.0_8, &
                           qc2 = - 25.0_8 / 384.0_8, &
                           qc3 =    3.0_8 / 640.0_8

type, public, extends(sbp_diff_t) :: sbp_diff_i2c_63_t
contains
    procedure, public :: apply_tile
    procedure, public :: add_SAT_correction
end type sbp_diff_i2c_63_t

contains

subroutine apply_tile(this, f_out, f_in, mesh, direction)
    class(sbp_diff_i2c_63_t), intent(in)    :: this
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
                    f_out%p(1,j,k) = q11*f_in%p(1,j,k)+q12*f_in%p(2,j,k)+q13*f_in%p(3,j,k)+q14*f_in%p(4,j,k)+q15*f_in%p(5,j,k)
                end if
                if (mesh%is<=2 .and. mesh%ie>=2) then
                    f_out%p(2,j,k) = q21*f_in%p(1,j,k)+q22*f_in%p(2,j,k)+q23*f_in%p(3,j,k)+q24*f_in%p(4,j,k)+q25*f_in%p(5,j,k)
                end if
                if (mesh%is<=3 .and. mesh%ie>=3) then
                    f_out%p(3,j,k) = q31*f_in%p(1,j,k)+q32*f_in%p(2,j,k)+q33*f_in%p(3,j,k)+q34*f_in%p(4,j,k)+q35*f_in%p(5,j,k)+q36*f_in%p(6,j,k)
                end if
                if (mesh%is<=4 .and. mesh%ie>=4) then
                    f_out%p(4,j,k) = q41*f_in%p(1,j,k)+q42*f_in%p(2,j,k)+q43*f_in%p(3,j,k)+q44*f_in%p(4,j,k)+q45*f_in%p(5,j,k)+q46*f_in%p(6,j,k)+q47*f_in%p(7,j,k)
                end if
                if (mesh%is<=5 .and. mesh%ie>=5) then
                    f_out%p(5,j,k) = q51*f_in%p(1,j,k)+q52*f_in%p(2,j,k)+q53*f_in%p(3,j,k)+q54*f_in%p(4,j,k)+q55*f_in%p(5,j,k)+q56*f_in%p(6,j,k)+q57*f_in%p(7,j,k)+q58*f_in%p(8,j,k)
                end if
                if (mesh%is<=6 .and. mesh%ie>=6) then
                    f_out%p(6,j,k) = q61*f_in%p(1,j,k)+q62*f_in%p(2,j,k)+q63*f_in%p(3,j,k)+q64*f_in%p(4,j,k)+q65*f_in%p(5,j,k)+q66*f_in%p(6,j,k)+q67*f_in%p(7,j,k)+q68*f_in%p(8,j,k)+q69*f_in%p(9,j,k)
                end if
                do i = max(mesh%is, 7), min(mesh%ie, mesh%nx-6)
                    f_out%p(i,j,k) = qc1*(f_in%p(i+1,j,k) - f_in%p(i  ,j,k)) &
                                   + qc2*(f_in%p(i+2,j,k) - f_in%p(i-1,j,k)) &
                                   + qc3*(f_in%p(i+3,j,k) - f_in%p(i-2,j,k))
                end do
                if (mesh%is<=mesh%nx-5 .and. mesh%ie>=mesh%nx-5) then
                    f_out%p(mesh%nx-5,j,k) = -q61*f_in%p(mesh%nx+1,j,k)-q62*f_in%p(mesh%nx,j,k)-q63*f_in%p(mesh%nx-1,j,k)-q64*f_in%p(mesh%nx-2,j,k)-q65*f_in%p(mesh%nx-3,j,k)-q66*f_in%p(mesh%nx-4,j,k)-q67*f_in%p(mesh%nx-5,j,k)-q68*f_in%p(mesh%nx-6,j,k)-q69*f_in%p(mesh%nx-7,j,k)
                end if
                if (mesh%is<=mesh%nx-4 .and. mesh%ie>=mesh%nx-4) then
                    f_out%p(mesh%nx-4,j,k) = -q51*f_in%p(mesh%nx+1,j,k)-q52*f_in%p(mesh%nx,j,k)-q53*f_in%p(mesh%nx-1,j,k)-q54*f_in%p(mesh%nx-2,j,k)-q55*f_in%p(mesh%nx-3,j,k)-q56*f_in%p(mesh%nx-4,j,k)-q57*f_in%p(mesh%nx-5,j,k)-q58*f_in%p(mesh%nx-6,j,k)
                end if
                if (mesh%is<=mesh%nx-3 .and. mesh%ie>=mesh%nx-3) then
                    f_out%p(mesh%nx-3,j,k) = -q41*f_in%p(mesh%nx+1,j,k)-q42*f_in%p(mesh%nx,j,k)-q43*f_in%p(mesh%nx-1,j,k)-q44*f_in%p(mesh%nx-2,j,k)-q45*f_in%p(mesh%nx-3,j,k)-q46*f_in%p(mesh%nx-4,j,k)-q47*f_in%p(mesh%nx-5,j,k)
                end if
                if (mesh%is<=mesh%nx-2 .and. mesh%ie>=mesh%nx-2) then
                    f_out%p(mesh%nx-2,j,k) = -q31*f_in%p(mesh%nx+1,j,k)-q32*f_in%p(mesh%nx,j,k)-q33*f_in%p(mesh%nx-1,j,k)-q34*f_in%p(mesh%nx-2,j,k)-q35*f_in%p(mesh%nx-3,j,k)-q36*f_in%p(mesh%nx-4,j,k)
                end if
                if (mesh%is<=mesh%nx-1 .and. mesh%ie>=mesh%nx-1) then
                    f_out%p(mesh%nx-1,j,k) = -q21*f_in%p(mesh%nx+1,j,k)-q22*f_in%p(mesh%nx,j,k)-q23*f_in%p(mesh%nx-1,j,k)-q24*f_in%p(mesh%nx-2,j,k)-q25*f_in%p(mesh%nx-3,j,k)
                end if
                if (mesh%is<=mesh%nx .and. mesh%ie>=mesh%nx) then
                    f_out%p(mesh%nx,j,k) = -q11*f_in%p(mesh%nx+1,j,k)-q12*f_in%p(mesh%nx,j,k)-q13*f_in%p(mesh%nx-1,j,k)-q14*f_in%p(mesh%nx-2,j,k)-q15*f_in%p(mesh%nx-3,j,k)
                end if
            end do
        end do
    case("y")
        do k = mesh%ks, mesh%ke
            if (mesh%js<=1 .and. mesh%je>=1) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,1,k) = q11*f_in%p(i,1,k)+q12*f_in%p(i,2,k)+q13*f_in%p(i,3,k)+q14*f_in%p(i,4,k)+q15*f_in%p(i,5,k)
                end do
            end if
            if (mesh%js<=2 .and. mesh%je>=2) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,2,k) = q21*f_in%p(i,1,k)+q22*f_in%p(i,2,k)+q23*f_in%p(i,3,k)+q24*f_in%p(i,4,k)+q25*f_in%p(i,5,k)
                end do
            end if
            if (mesh%js<=3 .and. mesh%je>=3) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,3,k) = q31*f_in%p(i,1,k)+q32*f_in%p(i,2,k)+q33*f_in%p(i,3,k)+q34*f_in%p(i,4,k)+q35*f_in%p(i,5,k)+q36*f_in%p(i,6,k)
                end do
            end if
            if (mesh%js<=4 .and. mesh%je>=4) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,4,k) = q41*f_in%p(i,1,k)+q42*f_in%p(i,2,k)+q43*f_in%p(i,3,k)+q44*f_in%p(i,4,k)+q45*f_in%p(i,5,k)+q46*f_in%p(i,6,k)+q47*f_in%p(i,7,k)
                end do
            end if
            if (mesh%js<=5 .and. mesh%je>=5) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,5,k) = q51*f_in%p(i,1,k)+q52*f_in%p(i,2,k)+q53*f_in%p(i,3,k)+q54*f_in%p(i,4,k)+q55*f_in%p(i,5,k)+q56*f_in%p(i,6,k)+q57*f_in%p(i,7,k)+q58*f_in%p(i,8,k)
                end do
            end if
            if (mesh%js<=6 .and. mesh%je>=6) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,6,k) = q61*f_in%p(i,1,k)+q62*f_in%p(i,2,k)+q63*f_in%p(i,3,k)+q64*f_in%p(i,4,k)+q65*f_in%p(i,5,k)+q66*f_in%p(i,6,k)+q67*f_in%p(i,7,k)+q68*f_in%p(i,8,k)+q69*f_in%p(i,9,k)
                end do
            end if
            do j = max(mesh%js, 7), min(mesh%je, mesh%ny-6)
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,k) = qc1*(f_in%p(i,j+1,k) - f_in%p(i,j  ,k)) &
                                   + qc2*(f_in%p(i,j+2,k) - f_in%p(i,j-1,k)) &
                                   + qc3*(f_in%p(i,j+3,k) - f_in%p(i,j-2,k))
                end do 
            end do
            if (mesh%js<=mesh%ny-5 .and. mesh%je>=mesh%ny-5) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-5,k) = -q61*f_in%p(i,mesh%ny+1,k)-q62*f_in%p(i,mesh%ny,k)-q63*f_in%p(i,mesh%ny-1,k)-q64*f_in%p(i,mesh%ny-2,k)-q65*f_in%p(i,mesh%ny-3,k)-q66*f_in%p(i,mesh%ny-4,k)-q67*f_in%p(i,mesh%ny-5,k)-q68*f_in%p(i,mesh%ny-6,k)-q69*f_in%p(i,mesh%ny-7,k)
                end do
            end if
            if (mesh%js<=mesh%ny-4 .and. mesh%je>=mesh%ny-4) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-4,k) = -q51*f_in%p(i,mesh%ny+1,k)-q52*f_in%p(i,mesh%ny,k)-q53*f_in%p(i,mesh%ny-1,k)-q54*f_in%p(i,mesh%ny-2,k)-q55*f_in%p(i,mesh%ny-3,k)-q56*f_in%p(i,mesh%ny-4,k)-q57*f_in%p(i,mesh%ny-5,k)-q58*f_in%p(i,mesh%ny-6,k)
                end do
            end if
            if (mesh%js<=mesh%ny-3 .and. mesh%je>=mesh%ny-3) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-3,k) = -q41*f_in%p(i,mesh%ny+1,k)-q42*f_in%p(i,mesh%ny,k)-q43*f_in%p(i,mesh%ny-1,k)-q44*f_in%p(i,mesh%ny-2,k)-q45*f_in%p(i,mesh%ny-3,k)-q46*f_in%p(i,mesh%ny-4,k)-q47*f_in%p(i,mesh%ny-5,k)
                end do
            end if
            if (mesh%js<=mesh%ny-2 .and. mesh%je>=mesh%ny-2) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-2,k) = -q31*f_in%p(i,mesh%ny+1,k)-q32*f_in%p(i,mesh%ny,k)-q33*f_in%p(i,mesh%ny-1,k)-q34*f_in%p(i,mesh%ny-2,k)-q35*f_in%p(i,mesh%ny-3,k)-q36*f_in%p(i,mesh%ny-4,k)
                end do
            end if
            if (mesh%js<=mesh%ny-1 .and. mesh%je>=mesh%ny-1) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-1,k) = -q21*f_in%p(i,mesh%ny+1,k)-q22*f_in%p(i,mesh%ny,k)-q23*f_in%p(i,mesh%ny-1,k)-q24*f_in%p(i,mesh%ny-2,k)-q25*f_in%p(i,mesh%ny-3,k)
                end do
            end if
            if (mesh%js<=mesh%ny .and. mesh%je>=mesh%ny) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny,k) = -q11*f_in%p(i,mesh%ny+1,k)-q12*f_in%p(i,mesh%ny,k)-q13*f_in%p(i,mesh%ny-1,k)-q14*f_in%p(i,mesh%ny-2,k)-q15*f_in%p(i,mesh%ny-3,k)
                end do
            end if
        end do
    case("z")
        if (mesh%ks<=1 .and. mesh%ke>=1) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,1) = q11*f_in%p(i,j,1)+q12*f_in%p(i,j,2)+q13*f_in%p(i,j,3)+q14*f_in%p(i,j,4)+q15*f_in%p(i,j,5)
                end do
            end do
        end if
        if (mesh%ks<=2 .and. mesh%ke>=2) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,2) = q21*f_in%p(i,j,1)+q22*f_in%p(i,j,2)+q23*f_in%p(i,j,3)+q24*f_in%p(i,j,4)+q25*f_in%p(i,j,5)
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
                    f_out%p(i,j,4) = q41*f_in%p(i,j,1)+q42*f_in%p(i,j,2)+q43*f_in%p(i,j,3)+q44*f_in%p(i,j,4)+q45*f_in%p(i,j,5)+q46*f_in%p(i,j,6)+q47*f_in%p(i,j,7)
                end do
            end do
        end if
        if (mesh%ks<=5 .and. mesh%ke>=5) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,5) = q51*f_in%p(i,j,1)+q52*f_in%p(i,j,2)+q53*f_in%p(i,j,3)+q54*f_in%p(i,j,4)+q55*f_in%p(i,j,5)+q56*f_in%p(i,j,6)+q57*f_in%p(i,j,7)+q58*f_in%p(i,j,8)
                end do
            end do
        end if
        if (mesh%ks<=6 .and. mesh%ke>=6) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,6) = q61*f_in%p(i,j,1)+q62*f_in%p(i,j,2)+q63*f_in%p(i,j,3)+q64*f_in%p(i,j,4)+q65*f_in%p(i,j,5)+q66*f_in%p(i,j,6)+q67*f_in%p(i,j,7)+q68*f_in%p(i,j,8)+q69*f_in%p(i,j,9)
                end do
            end do
        end if
        do k = max(mesh%ks, 7), min(mesh%ke, mesh%nz-6)
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,k) = qc1*(f_in%p(i,j,k+1) - f_in%p(i,j,k  )) &
                                   + qc2*(f_in%p(i,j,k+2) - f_in%p(i,j,k-1)) & 
                                   + qc3*(f_in%p(i,j,k+3) - f_in%p(i,j,k-2))
                end do
            end do
        end do
        if (mesh%ks<=mesh%nz-5 .and. mesh%ke>=mesh%nz-5) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-5) = -q61*f_in%p(i,j,mesh%nz+1)-q62*f_in%p(i,j,mesh%nz)-q63*f_in%p(i,j,mesh%nz-1)-q64*f_in%p(i,j,mesh%nz-2)-q65*f_in%p(i,j,mesh%nz-3)-q66*f_in%p(i,j,mesh%nz-4)-q67*f_in%p(i,j,mesh%nz-5)-q68*f_in%p(i,j,mesh%nz-6)-q69*f_in%p(i,j,mesh%nz-7)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz-4 .and. mesh%ke>=mesh%nz-4) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-4) = -q51*f_in%p(i,j,mesh%nz+1)-q52*f_in%p(i,j,mesh%nz)-q53*f_in%p(i,j,mesh%nz-1)-q54*f_in%p(i,j,mesh%nz-2)-q55*f_in%p(i,j,mesh%nz-3)-q56*f_in%p(i,j,mesh%nz-4)-q57*f_in%p(i,j,mesh%nz-5)-q58*f_in%p(i,j,mesh%nz-6)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz-3 .and. mesh%ke>=mesh%nz-3) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-3) = -q41*f_in%p(i,j,mesh%nz+1)-q42*f_in%p(i,j,mesh%nz)-q43*f_in%p(i,j,mesh%nz-1)-q44*f_in%p(i,j,mesh%nz-2)-q45*f_in%p(i,j,mesh%nz-3)-q46*f_in%p(i,j,mesh%nz-4)-q47*f_in%p(i,j,mesh%nz-5)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz-2 .and. mesh%ke>=mesh%nz-2) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-2) = -q31*f_in%p(i,j,mesh%nz+1)-q32*f_in%p(i,j,mesh%nz)-q33*f_in%p(i,j,mesh%nz-1)-q34*f_in%p(i,j,mesh%nz-2)-q35*f_in%p(i,j,mesh%nz-3)-q36*f_in%p(i,j,mesh%nz-4)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz-1 .and. mesh%ke>=mesh%nz-1) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-1) = -q21*f_in%p(i,j,mesh%nz+1)-q22*f_in%p(i,j,mesh%nz)-q23*f_in%p(i,j,mesh%nz-1)-q24*f_in%p(i,j,mesh%nz-2)-q25*f_in%p(i,j,mesh%nz-3)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz .and. mesh%ke>=mesh%nz) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz) = -q11*f_in%p(i,j,mesh%nz+1)-q12*f_in%p(i,j,mesh%nz)-q13*f_in%p(i,j,mesh%nz-1)-q14*f_in%p(i,j,mesh%nz-2)-q15*f_in%p(i,j,mesh%nz-3)
                end do
            end do
        end if
    case default
        call parcomm_global%abort("Wrong direction in sbp_diff_i2c_42_mod")
    end select

end subroutine apply_tile

subroutine add_SAT_correction(this, f_out, f_in, mesh, direction)
    class(sbp_diff_i2c_63_t), intent(in)    :: this
    type(tile_field_t),       intent(inout) :: f_out
    type(tile_field_t),       intent(in)    :: f_in
    type(tile_mesh_t),        intent(in)    :: mesh
    character(len=*),         intent(in)    :: direction

    integer(kind=4) :: k, i, j

    select case(direction)
    case("x")
        if (mesh%is<=1 .and. mesh%ie>=1) then
            do k = mesh%ks, mesh%ke
                do j = mesh%js, mesh%je
                    f_out%p(1,j,k) = f_out%p(1,j,k) &
                                   + 0.5_8*R11/h11*(f_in%p(1,j,k)-f_in%p(0,j,k))
                end do
            end do
        end if
        if (mesh%is<=2 .and. mesh%ie>=2) then
            do k = mesh%ks, mesh%ke
                do j = mesh%js, mesh%je
                    f_out%p(2,j,k) = f_out%p(2,j,k) &
                                   + 0.5_8*R12/h22*(f_in%p(1,j,k)-f_in%p(0,j,k))
                end do
            end do
        end if
        if (mesh%is<=3 .and. mesh%ie>=3) then
            do k = mesh%ks, mesh%ke
                do j = mesh%js, mesh%je
                    f_out%p(3,j,k) = f_out%p(3,j,k) &
                                   + 0.5_8*R13/h33*(f_in%p(1,j,k)-f_in%p(0,j,k))
                end do
            end do
        end if
        if (mesh%is<=4 .and. mesh%ie>=4) then
            do k = mesh%ks, mesh%ke
                do j = mesh%js, mesh%je
                    f_out%p(4,j,k) = f_out%p(4,j,k) &
                                   + 0.5_8*R14/h44*(f_in%p(1,j,k)-f_in%p(0,j,k))
                end do
            end do
        end if
        if (mesh%is<=mesh%nx .and. mesh%ie>=mesh%nx) then
            do k = mesh%ks, mesh%ke
                do j = mesh%js, mesh%je
                    f_out%p(mesh%nx,j,k) = f_out%p(mesh%nx,j,k) &
                                   + 0.5_8*R11/h11*(f_in%p(mesh%nx+2,j,k)-f_in%p(mesh%nx+1,j,k))
                end do
            end do
        end if
        if (mesh%is<=mesh%nx-1 .and. mesh%ie>=mesh%nx-1) then
            do k = mesh%ks, mesh%ke
                do j = mesh%js, mesh%je
                    f_out%p(mesh%nx-1,j,k) = f_out%p(mesh%nx-1,j,k) &
                                   + 0.5_8*R12/h22*(f_in%p(mesh%nx+2,j,k)-f_in%p(mesh%nx+1,j,k))
                end do
            end do
        end if
        if (mesh%is<=mesh%nx-2 .and. mesh%ie>=mesh%nx-2) then
            do k = mesh%ks, mesh%ke
                do j = mesh%js, mesh%je
                    f_out%p(mesh%nx-2,j,k) = f_out%p(mesh%nx-2,j,k) &
                                   + 0.5_8*R13/h33*(f_in%p(mesh%nx+2,j,k)-f_in%p(mesh%nx+1,j,k))
                end do
            end do
        end if
        if (mesh%is<=mesh%nx-3 .and. mesh%ie>=mesh%nx-3) then
            do k = mesh%ks, mesh%ke
                do j = mesh%js, mesh%je
                    f_out%p(mesh%nx-3,j,k) = f_out%p(mesh%nx-3,j,k) &
                                   + 0.5_8*R14/h44*(f_in%p(mesh%nx+2,j,k)-f_in%p(mesh%nx+1,j,k))
                end do
            end do
        end if
    case("y")
        if (mesh%js<=1 .and. mesh%je>=1) then
            do k = mesh%ks, mesh%ke
                do i = mesh%is, mesh%ie
                    f_out%p(i,1,k) = f_out%p(i,1,k) &
                                   + 0.5_8*R11/h11*(f_in%p(i,1,k)-f_in%p(i,0,k))
                end do
            end do
        end if
        if (mesh%js<=2 .and. mesh%je>=2) then
            do k = mesh%ks, mesh%ke
                do i = mesh%is, mesh%ie
                    f_out%p(i,2,k) = f_out%p(i,2,k) &
                                   + 0.5_8*R12/h22*(f_in%p(i,1,k)-f_in%p(i,0,k))
                end do
            end do
        end if
        if (mesh%js<=3 .and. mesh%je>=3) then
            do k = mesh%ks, mesh%ke
                do i = mesh%is, mesh%ie
                    f_out%p(i,3,k) = f_out%p(i,3,k) &
                                   + 0.5_8*R13/h33*(f_in%p(i,1,k)-f_in%p(i,0,k))
                end do
            end do
        end if
        if (mesh%js<=4 .and. mesh%je>=4) then
            do k = mesh%ks, mesh%ke
                do i = mesh%is, mesh%ie
                    f_out%p(i,4,k) = f_out%p(i,4,k) &
                                   + 0.5_8*R14/h44*(f_in%p(i,1,k)-f_in%p(i,0,k))
                end do
            end do
        end if
        if (mesh%js<=mesh%ny .and. mesh%je>=mesh%ny) then
            do k = mesh%ks, mesh%ke
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny,k) = f_out%p(i,mesh%ny,k) &
                                   + 0.5_8*R11/h11*(f_in%p(i,mesh%ny+2,k)-f_in%p(i,mesh%ny+1,k))
                end do
            end do
        end if
        if (mesh%js<=mesh%ny-1 .and. mesh%je>=mesh%ny-1) then
            do k = mesh%ks, mesh%ke
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-1,k) = f_out%p(i,mesh%ny-1,k) &
                                   + 0.5_8*R12/h22*(f_in%p(i,mesh%ny+2,k)-f_in%p(i,mesh%ny+1,k))
                end do
            end do
        end if
        if (mesh%js<=mesh%ny-2 .and. mesh%je>=mesh%ny-2) then
            do k = mesh%ks, mesh%ke
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-2,k) = f_out%p(i,mesh%ny-2,k) &
                                   + 0.5_8*R13/h33*(f_in%p(i,mesh%ny+2,k)-f_in%p(i,mesh%ny+1,k))
                end do
            end do
        end if
        if (mesh%js<=mesh%ny-3 .and. mesh%je>=mesh%ny-3) then
            do k = mesh%ks, mesh%ke
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-3,k) = f_out%p(i,mesh%ny-3,k) &
                                   + 0.5_8*R14/h44*(f_in%p(i,mesh%ny+2,k)-f_in%p(i,mesh%ny+1,k))
                end do
            end do
        end if
    case default
        call parcomm_global%abort("Wrong direction in add_SAT_correction in sbp_diff_i2c_63_mod")
    end select

end subroutine add_SAT_correction
end module sbp_diff_i2c_63_mod
