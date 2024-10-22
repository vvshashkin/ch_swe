module sbp_interp_i2c_63_mod

use sbp_interp_mod, only : sbp_interp_t
use grid_field_mod, only : grid_field_t, tile_field_t
use mesh_mod,       only : tile_mesh_t, mesh_t
use tile_mod,       only : tile_t
use parcomm_mod,    only : parcomm_global

implicit none

! c42, c43, c52, c53, c62, c64 are free parameters left after imposing
! SBP-preserving interpolation propertie. Current values obtained from
! successive minimization of the ep3^T*Hp*ep3 + eu3^T*Hu*eu3
! and  ep4^T*Hp*ep4 + eu4^T*Hu*eu4 functionals, 
! where eu3, ep3 and eu4, ep4 are the 3rd and 4th order polynomial interpolation error vectors
! i.e. e3u = Ip2u * xp^3 - xu^3, e3p = Iu2p * xu^3 - xp^3, 
! e4u = Ip2u * xp^4 - xu^4, e4p = Iu2p * xu^4 - xp^4

real(kind=8), parameter :: c42 =-0.3332211159670528_8,& !-0.14404519153844240456_8, & ! -926292970421086411583174718907177714793225786763096390206657584/6430571965145256254841052294837532198444291389539149195700255513
                           c43 = 0.3310769312612241_8,&!!0.032475525312611441344_8, & ! 290554716700722460974523199025570338981553700029589292713683817/8946882734115139137170159714556566536965970628924033663582964192
                           c52 =-0.07099703081266314_8,&!-0.080803485810270262054_8, & !-2137409560167516896903956685900349819708778665363881855560396625/26451947446750477499563970224543051267837607720695559552004539208, &
                           c53 =-0.2916164053358880_8, &!-0.14653181730277016544_8, & ! -674095987925168472730843075283023612247217159194874438740472431/4600338686391387391228516560790095872667410038381836443826876384
                           c62 = 0.05753938634775091_8, & !0.040209557114396073708_8, & ! 14329012774676524375994568354759783135527368981660032821701549/356358383503467209802084561878312385801436827978869550322912860
                           c64 =-0.1230378129758785_8!-0.044522001860911594323_8  !-95194731680956888263984476758232028156063030229847764949053621/2138150301020803258812507371269874314808620967873217301937477160

real(kind=8), parameter :: q11 = 4474753.0_8/7808712.0_8+312623.0_8/650726.0_8*c53+937869.0_8/650726.0_8*c52-813687.0_8/1301452.0_8*c64+2441061.0_8/1301452.0_8*c62+86857.0_8/325363.0_8*c42+86857.0_8/976089.0_8*c43, &
                           q12 = 136944.0_8/325363.0_8-173714.0_8/325363.0_8*c42-1627374.0_8/325363.0_8*c62-937869.0_8/325363.0_8*c52, &
                           q13 = -848457.0_8/2602904.0_8-173714.0_8/325363.0_8*c43+2441061.0_8/650726.0_8*c64+2441061.0_8/650726.0_8*c62-937869.0_8/325363.0_8*c53, &
                           q14 = 2331127.0_8/3904356.0_8+173714.0_8/325363.0_8*c42+694856.0_8/976089.0_8*c43-1627374.0_8/325363.0_8*c64+1250492.0_8/325363.0_8*c53+937869.0_8/325363.0_8*c52, &
                           q15 = -10145.0_8/38278.0_8-937869.0_8/650726.0_8*c53-937869.0_8/650726.0_8*c52+2441061.0_8/1301452.0_8*c64-813687.0_8/1301452.0_8*c62-86857.0_8/325363.0_8*c42-86857.0_8/325363.0_8*c43, &
                           q21 = -373145.0_8/354464.0_8-86857.0_8/144001.0_8*c43+4068435.0_8/1152008.0_8*c64-12205305.0_8/1152008.0_8*c62-1250492.0_8/432003.0_8*c53-260571.0_8/144001.0_8*c42-1250492.0_8/144001.0_8*c52, &
                           q22 = 273888.0_8/144001.0_8+521142.0_8/144001.0_8*c42+4068435.0_8/144001.0_8*c62+2500984.0_8/144001.0_8*c52, &
                           q23 = 520551.0_8/209456.0_8+521142.0_8/144001.0_8*c43-12205305.0_8/576004.0_8*c64-12205305.0_8/576004.0_8*c62+2500984.0_8/144001.0_8*c53, &
                           q24 = -1142117.0_8/288002.0_8-521142.0_8/144001.0_8*c42-694856.0_8/144001.0_8*c43+4068435.0_8/144001.0_8*c64-10003936.0_8/432003.0_8*c53-2500984.0_8/144001.0_8*c52, &
                           q25 = 7516251.0_8/4608032.0_8+260571.0_8/144001.0_8*c43-12205305.0_8/1152008.0_8*c64+4068435.0_8/1152008.0_8*c62+1250492.0_8/144001.0_8*c53+260571.0_8/144001.0_8*c42+1250492.0_8/144001.0_8*c52, &
                           q31 = 930131.0_8/6911200.0_8+312623.0_8/431950.0_8*c53-271229.0_8/345560.0_8*c64+86857.0_8/431950.0_8*c43+813687.0_8/345560.0_8*c62+937869.0_8/431950.0_8*c52+260571.0_8/431950.0_8*c42, &
                           q32 = -22824.0_8/215975.0_8-260571.0_8/215975.0_8*c42-271229.0_8/43195.0_8*c62-937869.0_8/215975.0_8*c52, &
                           q33 = -424911.0_8/3455600.0_8-937869.0_8/215975.0_8*c53+813687.0_8/172780.0_8*c64-260571.0_8/215975.0_8*c43+813687.0_8/172780.0_8*c62, &
                           q34 = 330902.0_8/215975.0_8+1250492.0_8/215975.0_8*c53+937869.0_8/215975.0_8*c52-271229.0_8/43195.0_8*c64+260571.0_8/215975.0_8*c42+347428.0_8/215975.0_8*c43, &
                           q35 = -615889.0_8/1382240.0_8-937869.0_8/431950.0_8*c53+813687.0_8/345560.0_8*c64-260571.0_8/431950.0_8*c43-271229.0_8/345560.0_8*c62-937869.0_8/431950.0_8*c52-260571.0_8/431950.0_8*c42, &
                           q36 = 324.0_8/43195.0_8, &
                           q41 = -17737.0_8/4169136.0_8-1.0_8/2.0_8*c42-1.0_8/6.0_8*c43, &
                           q42 = c42, &
                           q43 = c43, &
                           q44 = 415759.0_8/1042284.0_8-c42-4.0_8/3.0_8*c43, &
                           q45 = 1031359.0_8/1389712.0_8+1.0_8/2.0_8*c42+1.0_8/2.0_8*c43, &
                           q46 = -13500.0_8/86857.0_8, &
                           q47 = 1620.0_8/86857.0_8, &
                           q51 = 44783.0_8/5001968.0_8-1.0_8/6.0_8*c53-1.0_8/2.0_8*c52, &
                           q52 = c52, &
                           q53 = c53, &
                           q54 = -199149.0_8/1250492.0_8-4.0_8/3.0_8*c53-c52, &
                           q55 = 3541941.0_8/5001968.0_8+1.0_8/2.0_8*c53+1.0_8/2.0_8*c52, &
                           q56 = 162000.0_8/312623.0_8, &
                           q57 = -27000.0_8/312623.0_8, &
                           q58 = 3240.0_8/312623.0_8, &
                           q61 = -123231.0_8/8679328.0_8+1.0_8/8.0_8*c64-3.0_8/8.0_8*c62, &
                           q62 = c62, &
                           q63 = 30309.0_8/619952.0_8-3.0_8/4.0_8*c64-3.0_8/4.0_8*c62, &
                           q64 = c64, &
                           q65 = -1229447.0_8/8679328.0_8-3.0_8/8.0_8*c64+1.0_8/8.0_8*c62, &
                           q66 = 162000.0_8/271229.0_8, &
                           q67 = 162000.0_8/271229.0_8, &
                           q68 = -27000.0_8/271229.0_8, &
                           q69 = 3240.0_8/271229.0_8                           

real(kind=8), parameter :: qc1 = 150.0_8 / 256.0_8, &
                           qc2 = -25.0_8 / 256.0_8, &
                           qc3 =   3.0_8 / 256.0_8

type, public, extends(sbp_interp_t) :: sbp_interp_i2c_63_t
contains
    procedure, public :: apply_tile
end type sbp_interp_i2c_63_t

contains

subroutine apply_tile(this, f_out, f_in, mesh, direction)

    class(sbp_interp_i2c_63_t), intent(in)    :: this
    type(tile_field_t),         intent(inout) :: f_out
    type(tile_field_t),         intent(in)    :: f_in
    type(tile_mesh_t),          intent(in)    :: mesh
    character(len=*),           intent(in)    :: direction

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
                    f_out%p(i,j,k) = qc1*(f_in%p(i+1,j,k) + f_in%p(i  ,j,k)) &
                                   + qc2*(f_in%p(i+2,j,k) + f_in%p(i-1,j,k)) &
                                   + qc3*(f_in%p(i+3,j,k) + f_in%p(i-2,j,k))
                end do
                if (mesh%is<=mesh%nx-5 .and. mesh%ie>=mesh%nx-5) then
                    f_out%p(mesh%nx-5,j,k) = q61*f_in%p(mesh%nx+1,j,k)+q62*f_in%p(mesh%nx,j,k)+q63*f_in%p(mesh%nx-1,j,k)+q64*f_in%p(mesh%nx-2,j,k)+q65*f_in%p(mesh%nx-3,j,k) &
                                           + q66*f_in%p(mesh%nx-4,j,k)+q67*f_in%p(mesh%nx-5,j,k)+q68*f_in%p(mesh%nx-6,j,k)+q69*f_in%p(mesh%nx-7,j,k)
                end if
                if (mesh%is<=mesh%nx-4 .and. mesh%ie>=mesh%nx-4) then
                    f_out%p(mesh%nx-4,j,k) = q51*f_in%p(mesh%nx+1,j,k)+q52*f_in%p(mesh%nx,j,k)+q53*f_in%p(mesh%nx-1,j,k)+q54*f_in%p(mesh%nx-2,j,k)+q55*f_in%p(mesh%nx-3,j,k) &
                                           + q56*f_in%p(mesh%nx-4,j,k)+q57*f_in%p(mesh%nx-5,j,k)+q58*f_in%p(mesh%nx-6,j,k)
                end if
                if (mesh%is<=mesh%nx-3 .and. mesh%ie>=mesh%nx-3) then
                    f_out%p(mesh%nx-3,j,k) = q41*f_in%p(mesh%nx+1,j,k)+q42*f_in%p(mesh%nx,j,k)+q43*f_in%p(mesh%nx-1,j,k)+q44*f_in%p(mesh%nx-2,j,k)+q45*f_in%p(mesh%nx-3,j,k) &
                                           + q46*f_in%p(mesh%nx-4,j,k)+q47*f_in%p(mesh%nx-5,j,k)
                end if
                if (mesh%is<=mesh%nx-2 .and. mesh%ie>=mesh%nx-2) then
                    f_out%p(mesh%nx-2,j,k) = q31*f_in%p(mesh%nx+1,j,k)+q32*f_in%p(mesh%nx,j,k)+q33*f_in%p(mesh%nx-1,j,k)+q34*f_in%p(mesh%nx-2,j,k)+q35*f_in%p(mesh%nx-3,j,k)+q36*f_in%p(mesh%nx-4,j,k)
                end if
                if (mesh%is<=mesh%nx-1 .and. mesh%ie>=mesh%nx-1) then
                    f_out%p(mesh%nx-1,j,k) = q21*f_in%p(mesh%nx+1,j,k)+q22*f_in%p(mesh%nx,j,k)+q23*f_in%p(mesh%nx-1,j,k)+q24*f_in%p(mesh%nx-2,j,k)+q25*f_in%p(mesh%nx-3,j,k)
                end if
                if (mesh%is<=mesh%nx .and. mesh%ie>=mesh%nx) then
                    f_out%p(mesh%nx,j,k) = q11*f_in%p(mesh%nx+1,j,k)+q12*f_in%p(mesh%nx,j,k)+q13*f_in%p(mesh%nx-1,j,k)+q14*f_in%p(mesh%nx-2,j,k)+q15*f_in%p(mesh%nx-3,j,k)
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
                    f_out%p(i,5,k) = q51*f_in%p(i,1,k)+q52*f_in%p(i,2,k)+q53*f_in%p(i,3,k)+q54*f_in%p(i,4,k)+q55*f_in%p(i,5,k)+q56*f_in%p(i,6,k)+q57*f_in%p(i,7,k) &
                                   + q58*f_in%p(i,8,k)
                end do
            end if
            if (mesh%js<=6 .and. mesh%je>=6) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,6,k) = q61*f_in%p(i,1,k)+q62*f_in%p(i,2,k)+q63*f_in%p(i,3,k)+q64*f_in%p(i,4,k)+q65*f_in%p(i,5,k)+q66*f_in%p(i,6,k)+q67*f_in%p(i,7,k) &
                                   + q68*f_in%p(i,8,k)+q69*f_in%p(i,9,k)
                end do
            end if
            do j = max(mesh%js, 7), min(mesh%je, mesh%ny-6)
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,k) = qc1*(f_in%p(i,j+1,k) + f_in%p(i,j  ,k)) &
                                   + qc2*(f_in%p(i,j+2,k) + f_in%p(i,j-1,k)) &
                                   + qc3*(f_in%p(i,j+3,k) + f_in%p(i,j-2,k))
                end do
            end do
            if (mesh%js<=mesh%ny-5 .and. mesh%je>=mesh%ny-5) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-5,k) = q61*f_in%p(i,mesh%ny+1,k)+q62*f_in%p(i,mesh%ny,k)+q63*f_in%p(i,mesh%ny-1,k)+q64*f_in%p(i,mesh%ny-2,k)+q65*f_in%p(i,mesh%ny-3,k) &
                                           + q66*f_in%p(i,mesh%ny-4,k)+q67*f_in%p(i,mesh%ny-5,k)+q68*f_in%p(i,mesh%ny-6,k)+q69*f_in%p(i,mesh%ny-7,k)
                end do
            end if
            if (mesh%js<=mesh%ny-4 .and. mesh%je>=mesh%ny-4) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-4,k) = q51*f_in%p(i,mesh%ny+1,k)+q52*f_in%p(i,mesh%ny,k)+q53*f_in%p(i,mesh%ny-1,k)+q54*f_in%p(i,mesh%ny-2,k)+q55*f_in%p(i,mesh%ny-3,k) &
                                           + q56*f_in%p(i,mesh%ny-4,k)+q57*f_in%p(i,mesh%ny-5,k)+q58*f_in%p(i,mesh%ny-6,k)
                end do
            end if
            if (mesh%js<=mesh%ny-3 .and. mesh%je>=mesh%ny-3) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-3,k) = q41*f_in%p(i,mesh%ny+1,k)+q42*f_in%p(i,mesh%ny,k)+q43*f_in%p(i,mesh%ny-1,k)+q44*f_in%p(i,mesh%ny-2,k)+q45*f_in%p(i,mesh%ny-3,k) &
                                           + q46*f_in%p(i,mesh%ny-4,k)+q47*f_in%p(i,mesh%ny-5,k)
                end do
            end if
            if (mesh%js<=mesh%ny-2 .and. mesh%je>=mesh%ny-2) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-2,k) = q31*f_in%p(i,mesh%ny+1,k)+q32*f_in%p(i,mesh%ny,k)+q33*f_in%p(i,mesh%ny-1,k)+q34*f_in%p(i,mesh%ny-2,k)+q35*f_in%p(i,mesh%ny-3,k)+q36*f_in%p(i,mesh%ny-4,k)
                end do
            end if
            if (mesh%js<=mesh%ny-1 .and. mesh%je>=mesh%ny-1) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-1,k) = q21*f_in%p(i,mesh%ny+1,k)+q22*f_in%p(i,mesh%ny,k)+q23*f_in%p(i,mesh%ny-1,k)+q24*f_in%p(i,mesh%ny-2,k)+q25*f_in%p(i,mesh%ny-3,k)
                end do
            end if
            if (mesh%js<=mesh%ny .and. mesh%je>=mesh%ny) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny,k) = q11*f_in%p(i,mesh%ny+1,k)+q12*f_in%p(i,mesh%ny,k)+q13*f_in%p(i,mesh%ny-1,k)+q14*f_in%p(i,mesh%ny-2,k)+q15*f_in%p(i,mesh%ny-3,k)
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
                    f_out%p(i,j,k) = qc1*(f_in%p(i,j,k+1) + f_in%p(i,j,k  )) &
                                   + qc2*(f_in%p(i,j,k+2) + f_in%p(i,j,k-1)) &
                                   + qc3*(f_in%p(i,j,k+3) + f_in%p(i,j,k-2))
                end do
            end do
        end do
        if (mesh%ks<=mesh%nz-5 .and. mesh%ke>=mesh%nz-5) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-5) = q61*f_in%p(i,j,mesh%nz+1)+q62*f_in%p(i,j,mesh%nz)+q63*f_in%p(i,j,mesh%nz-1)+q64*f_in%p(i,j,mesh%nz-2)+q65*f_in%p(i,j,mesh%nz-3)+q66*f_in%p(i,j,mesh%nz-4)+q67*f_in%p(i,j,mesh%nz-5)+q68*f_in%p(i,j,mesh%nz-6)+q69*f_in%p(i,j,mesh%nz-7)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz-4 .and. mesh%ke>=mesh%nz-4) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-4) = q51*f_in%p(i,j,mesh%nz+1)+q52*f_in%p(i,j,mesh%nz)+q53*f_in%p(i,j,mesh%nz-1)+q54*f_in%p(i,j,mesh%nz-2)+q55*f_in%p(i,j,mesh%nz-3)+q56*f_in%p(i,j,mesh%nz-4)+q57*f_in%p(i,j,mesh%nz-5)+q58*f_in%p(i,j,mesh%nz-6)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz-3 .and. mesh%ke>=mesh%nz-3) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-3) = q41*f_in%p(i,j,mesh%nz+1)+q42*f_in%p(i,j,mesh%nz)+q43*f_in%p(i,j,mesh%nz-1)+q44*f_in%p(i,j,mesh%nz-2)+q45*f_in%p(i,j,mesh%nz-3)+q46*f_in%p(i,j,mesh%nz-4)+q47*f_in%p(i,j,mesh%nz-5)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz-2 .and. mesh%ke>=mesh%nz-2) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-2) = q31*f_in%p(i,j,mesh%nz+1)+q32*f_in%p(i,j,mesh%nz)+q33*f_in%p(i,j,mesh%nz-1)+q34*f_in%p(i,j,mesh%nz-2)+q35*f_in%p(i,j,mesh%nz-3)+q36*f_in%p(i,j,mesh%nz-4)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz-1 .and. mesh%ke>=mesh%nz-1) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-1) = q21*f_in%p(i,j,mesh%nz+1)+q22*f_in%p(i,j,mesh%nz)+q23*f_in%p(i,j,mesh%nz-1)+q24*f_in%p(i,j,mesh%nz-2)+q25*f_in%p(i,j,mesh%nz-3)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz .and. mesh%ke>=mesh%nz) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz) = q11*f_in%p(i,j,mesh%nz+1)+q12*f_in%p(i,j,mesh%nz)+q13*f_in%p(i,j,mesh%nz-1)+q14*f_in%p(i,j,mesh%nz-2)+q15*f_in%p(i,j,mesh%nz-3)
                end do
            end do
        end if
    case default
        call parcomm_global%abort("Wrong direction in sbp_interp_i2c_63_mod")
    end select

end subroutine apply_tile

end module sbp_interp_i2c_63_mod