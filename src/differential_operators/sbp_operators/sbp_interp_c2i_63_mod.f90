module sbp_interp_c2i_63_mod

use sbp_interp_mod, only : sbp_interp_t
use grid_field_mod, only : grid_field_t, tile_field_t
use mesh_mod,       only : tile_mesh_t, mesh_t
use tile_mod,       only : tile_t
use parcomm_mod,    only : parcomm_global

use sbp_interp_i2c_63_mod, only : c42, c43, c52, c53, c62, c64

implicit none

real(kind=8), parameter :: q11 = 4474753.0_8/2188800.0_8+312623.0_8/182400.0_8*c53+312623.0_8/60800.0_8*c52-271229.0_8/121600.0_8*c64+813687.0_8/121600.0_8*c62+86857.0_8/91200.0_8*c42+86857.0_8/273600.0_8*c43, &
                           q12 = -970177.0_8/583680.0_8-86857.0_8/91200.0_8*c43+271229.0_8/48640.0_8*c64-813687.0_8/48640.0_8*c62-312623.0_8/68400.0_8*c53-86857.0_8/30400.0_8*c42-312623.0_8/22800.0_8*c52, &
                           q13 = 930131.0_8/1459200.0_8+312623.0_8/91200.0_8*c53-271229.0_8/72960.0_8*c64+86857.0_8/91200.0_8*c43+271229.0_8/24320.0_8*c62+312623.0_8/30400.0_8*c52+86857.0_8/30400.0_8*c42, &
                           q14 = -17737.0_8/2188800.0_8-86857.0_8/91200.0_8*c42-86857.0_8/273600.0_8*c43, &
                           q15 = 2357.0_8/76800.0_8-312623.0_8/547200.0_8*c53-312623.0_8/182400.0_8*c52, &
                           q16 = -41077.0_8/972800.0_8+271229.0_8/729600.0_8*c64-271229.0_8/243200.0_8*c62, &
                           q21 = 3.0_8/8.0_8-86857.0_8/182592.0_8*c42-271229.0_8/60864.0_8*c62-312623.0_8/121728.0_8*c52, &
                           q22 = 3.0_8/4.0_8+86857.0_8/60864.0_8*c42+1356145.0_8/121728.0_8*c62+312623.0_8/45648.0_8*c52, &
                           q23 = -1.0_8/8.0_8-86857.0_8/60864.0_8*c42-1356145.0_8/182592.0_8*c62-312623.0_8/60864.0_8*c52, &
                           q24 = 86857.0_8/182592.0_8*c42, &
                           q25 = 312623.0_8/365184.0_8*c52, &
                           q26 = 271229.0_8/365184.0_8*c62, &
                           q31 = -94273.0_8/188416.0_8-86857.0_8/105984.0_8*c43+271229.0_8/47104.0_8*c64+271229.0_8/47104.0_8*c62-312623.0_8/70656.0_8*c53, &
                           q32 = 636229.0_8/376832.0_8+86857.0_8/35328.0_8*c43-1356145.0_8/94208.0_8*c64-1356145.0_8/94208.0_8*c62+312623.0_8/26496.0_8*c53, &
                           q33 = -141637.0_8/565248.0_8-312623.0_8/35328.0_8*c53+1356145.0_8/141312.0_8*c64-86857.0_8/35328.0_8*c43+1356145.0_8/141312.0_8*c62, &
                           q34 = 86857.0_8/105984.0_8*c43, &
                           q35 = 312623.0_8/211968.0_8*c53, &
                           q36 = 70721.0_8/1130496.0_8-271229.0_8/282624.0_8*c64-271229.0_8/282624.0_8*c62, &
                           q41 = 2331127.0_8/3654144.0_8+86857.0_8/152256.0_8*c42+86857.0_8/114192.0_8*c43-271229.0_8/50752.0_8*c64+312623.0_8/76128.0_8*c53+312623.0_8/101504.0_8*c52, &
                           q42 = -1142117.0_8/609024.0_8-86857.0_8/50752.0_8*c42-86857.0_8/38064.0_8*c43+1356145.0_8/101504.0_8*c64-312623.0_8/28548.0_8*c53-312623.0_8/38064.0_8*c52, &
                           q43 = 12727.0_8/5856.0_8+312623.0_8/38064.0_8*c53+312623.0_8/50752.0_8*c52-1356145.0_8/152256.0_8*c64+86857.0_8/50752.0_8*c42+86857.0_8/38064.0_8*c43, &
                           q44 = 415759.0_8/1827072.0_8-86857.0_8/152256.0_8*c42-86857.0_8/114192.0_8*c43, &
                           q45 = -66383.0_8/406016.0_8-312623.0_8/228384.0_8*c53-312623.0_8/304512.0_8*c52, &
                           q46 = 271229.0_8/304512.0_8*c64, &
                           q51 = -172465.0_8/542592.0_8-312623.0_8/180864.0_8*c53-312623.0_8/180864.0_8*c52+271229.0_8/120576.0_8*c64-271229.0_8/361728.0_8*c62-86857.0_8/271296.0_8*c42-86857.0_8/271296.0_8*c43, &
                           q52 = 835139.0_8/964608.0_8+86857.0_8/90432.0_8*c43-1356145.0_8/241152.0_8*c64+1356145.0_8/723456.0_8*c62+312623.0_8/67824.0_8*c53+86857.0_8/90432.0_8*c42+312623.0_8/67824.0_8*c52, &
                           q53 = -3079445.0_8/4340736.0_8-312623.0_8/90432.0_8*c53+1356145.0_8/361728.0_8*c64-86857.0_8/90432.0_8*c43-1356145.0_8/1085184.0_8*c62-312623.0_8/90432.0_8*c52-86857.0_8/90432.0_8*c42, &
                           q54 = 1031359.0_8/2170368.0_8+86857.0_8/271296.0_8*c42+86857.0_8/271296.0_8*c43, &
                           q55 = 131183.0_8/160768.0_8+312623.0_8/542592.0_8*c53+312623.0_8/542592.0_8*c52, &
                           q56 = -1229447.0_8/8681472.0_8-271229.0_8/723456.0_8*c64+271229.0_8/2170368.0_8*c62, &
                           q57 = 15.0_8/1256.0_8

real(kind=8), parameter :: qc1 = 150.0_8 / 256.0_8, &
                           qc2 = -25.0_8 / 256.0_8, &
                           qc3 =   3.0_8 / 256.0_8

type, public, extends(sbp_interp_t) :: sbp_interp_c2i_63_t
contains
    procedure, public :: apply_tile
end type sbp_interp_c2i_63_t

contains
    
subroutine apply_tile(this, f_out, f_in, mesh, direction)
    class(sbp_interp_c2i_63_t), intent(in)    :: this
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
                    f_out%p(i,j,k) = qc1*(f_in%p(i  ,j,k) + f_in%p(i-1,j,k)) &
                                   + qc2*(f_in%p(i+1,j,k) + f_in%p(i-2,j,k)) &
                                   + qc3*(f_in%p(i+2,j,k) + f_in%p(i-3,j,k))
                end do
                if (mesh%is<=mesh%nx-4 .and. mesh%ie>=mesh%nx-4) then
                    f_out%p(mesh%nx-4,j,k) = q51*f_in%p(mesh%nx-1,j,k)+q52*f_in%p(mesh%nx-2,j,k)+q53*f_in%p(mesh%nx-3,j,k)+q54*f_in%p(mesh%nx-4,j,k)+q55*f_in%p(mesh%nx-5,j,k)+q56*f_in%p(mesh%nx-6,j,k)+q57*f_in%p(mesh%nx-7,j,k)
                end if
                if (mesh%is<=mesh%nx-3 .and. mesh%ie>=mesh%nx-3) then
                    f_out%p(mesh%nx-3,j,k) = q41*f_in%p(mesh%nx-1,j,k)+q42*f_in%p(mesh%nx-2,j,k)+q43*f_in%p(mesh%nx-3,j,k)+q44*f_in%p(mesh%nx-4,j,k)+q45*f_in%p(mesh%nx-5,j,k)+q46*f_in%p(mesh%nx-6,j,k)
                end if
                if (mesh%is<=mesh%nx-2 .and. mesh%ie>=mesh%nx-2) then
                    f_out%p(mesh%nx-2,j,k) = q31*f_in%p(mesh%nx-1,j,k)+q32*f_in%p(mesh%nx-2,j,k)+q33*f_in%p(mesh%nx-3,j,k)+q34*f_in%p(mesh%nx-4,j,k)+q35*f_in%p(mesh%nx-5,j,k)+q36*f_in%p(mesh%nx-6,j,k)
                end if
                if (mesh%is<=mesh%nx-1 .and. mesh%ie>=mesh%nx-1) then
                    f_out%p(mesh%nx-1,j,k) = q21*f_in%p(mesh%nx-1,j,k)+q22*f_in%p(mesh%nx-2,j,k)+q23*f_in%p(mesh%nx-3,j,k)+q24*f_in%p(mesh%nx-4,j,k)+q25*f_in%p(mesh%nx-5,j,k)+q26*f_in%p(mesh%nx-6,j,k)
                end if
                if (mesh%is<=mesh%nx .and. mesh%ie>=mesh%nx) then
                    f_out%p(mesh%nx  ,j,k) = q11*f_in%p(mesh%nx-1,j,k)+q12*f_in%p(mesh%nx-2,j,k)+q13*f_in%p(mesh%nx-3,j,k)+q14*f_in%p(mesh%nx-4,j,k)+q15*f_in%p(mesh%nx-5,j,k)+q16*f_in%p(mesh%nx-6,j,k)
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
                    f_out%p(i,j,k) = qc1*(f_in%p(i,j  ,k) + f_in%p(i,j-1,k)) &
                                   + qc2*(f_in%p(i,j+1,k) + f_in%p(i,j-2,k)) &
                                   + qc3*(f_in%p(i,j+2,k) + f_in%p(i,j-3,k))
                end do
            end do
            if (mesh%js<=mesh%ny-4 .and. mesh%je>=mesh%ny-4) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-4,k) = q51*f_in%p(i,mesh%ny-1,k)+q52*f_in%p(i,mesh%ny-2,k)+q53*f_in%p(i,mesh%ny-3,k)+q54*f_in%p(i,mesh%ny-4,k)+q55*f_in%p(i,mesh%ny-5,k)+q56*f_in%p(i,mesh%ny-6,k)+q57*f_in%p(i,mesh%ny-7,k)
                end do
            end if
            if (mesh%js<=mesh%ny-3 .and. mesh%je>=mesh%ny-3) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-3,k) = q41*f_in%p(i,mesh%ny-1,k)+q42*f_in%p(i,mesh%ny-2,k)+q43*f_in%p(i,mesh%ny-3,k)+q44*f_in%p(i,mesh%ny-4,k)+q45*f_in%p(i,mesh%ny-5,k)+q46*f_in%p(i,mesh%ny-6,k)
                end do
            end if
            if (mesh%js<=mesh%ny-2 .and. mesh%je>=mesh%ny-2) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-2,k) = q31*f_in%p(i,mesh%ny-1,k)+q32*f_in%p(i,mesh%ny-2,k)+q33*f_in%p(i,mesh%ny-3,k)+q34*f_in%p(i,mesh%ny-4,k)+q35*f_in%p(i,mesh%ny-5,k)+q36*f_in%p(i,mesh%ny-6,k)
                end do
            end if
            if (mesh%js<=mesh%ny-1 .and. mesh%je>=mesh%ny-1) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-1,k) = q21*f_in%p(i,mesh%ny-1,k)+q22*f_in%p(i,mesh%ny-2,k)+q23*f_in%p(i,mesh%ny-3,k)+q24*f_in%p(i,mesh%ny-4,k)+q25*f_in%p(i,mesh%ny-5,k)+q26*f_in%p(i,mesh%ny-6,k)
                end do
            end if
            if (mesh%js<=mesh%ny .and. mesh%je>=mesh%ny) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny  ,k) = q11*f_in%p(i,mesh%ny-1,k)+q12*f_in%p(i,mesh%ny-2,k)+q13*f_in%p(i,mesh%ny-3,k)+q14*f_in%p(i,mesh%ny-4,k)+q15*f_in%p(i,mesh%ny-5,k)+q16*f_in%p(i,mesh%ny-6,k)
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
                    f_out%p(i,j,k) = qc1*(f_in%p(i,j,k  ) + f_in%p(i,j,k-1)) &
                                   + qc2*(f_in%p(i,j,k+1) + f_in%p(i,j,k-2)) & 
                                   + qc3*(f_in%p(i,j,k+2) + f_in%p(i,j,k-3))
                end do
            end do
        end do
        if (mesh%ks<=mesh%nz-4 .and. mesh%ke>=mesh%nz-4) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-4) = q51*f_in%p(i,j,mesh%nz-1)+q52*f_in%p(i,j,mesh%nz-2)+q53*f_in%p(i,j,mesh%nz-3)+q54*f_in%p(i,j,mesh%nz-4)+q55*f_in%p(i,j,mesh%nz-5)+q56*f_in%p(i,j,mesh%nz-6)+q57*f_in%p(i,j,mesh%nz-7)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz-3 .and. mesh%ke>=mesh%nz-3) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-3) = q41*f_in%p(i,j,mesh%nz-1)+q42*f_in%p(i,j,mesh%nz-2)+q43*f_in%p(i,j,mesh%nz-3)+q44*f_in%p(i,j,mesh%nz-4)+q45*f_in%p(i,j,mesh%nz-5)+q46*f_in%p(i,j,mesh%nz-6)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz-2 .and. mesh%ke>=mesh%nz-2) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-2) = q31*f_in%p(i,j,mesh%nz-1)+q32*f_in%p(i,j,mesh%nz-2)+q33*f_in%p(i,j,mesh%nz-3)+q34*f_in%p(i,j,mesh%nz-4)+q35*f_in%p(i,j,mesh%nz-5)+q36*f_in%p(i,j,mesh%nz-6)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz-1 .and. mesh%ke>=mesh%nz-1) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-1) = q21*f_in%p(i,j,mesh%nz-1)+q22*f_in%p(i,j,mesh%nz-2)+q23*f_in%p(i,j,mesh%nz-3)+q24*f_in%p(i,j,mesh%nz-4)+q25*f_in%p(i,j,mesh%nz-5)+q26*f_in%p(i,j,mesh%nz-6)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz .and. mesh%ke>=mesh%nz) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-1) = q11*f_in%p(i,j,mesh%nz-1)+q12*f_in%p(i,j,mesh%nz-2)+q13*f_in%p(i,j,mesh%nz-3)+q14*f_in%p(i,j,mesh%nz-4)+q15*f_in%p(i,j,mesh%nz-5)+q16*f_in%p(i,j,mesh%nz-6)
                end do
            end do
        end if
    case default
        call parcomm_global%abort("Wrong direction in sbp_interp_c2i_63_mod")
    end select

end subroutine apply_tile

end module sbp_interp_c2i_63_mod