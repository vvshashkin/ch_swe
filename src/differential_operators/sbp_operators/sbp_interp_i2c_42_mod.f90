module sbp_interp_i2c_42_mod

use sbp_interp_mod, only : sbp_interp_t
use grid_field_mod, only : grid_field_t, tile_field_t
use mesh_mod,       only : tile_mesh_t, mesh_t
use tile_mod,       only : tile_t
use parcomm_mod,    only : parcomm_global


implicit none

! c13, c14 are free parameters left after imposing
! SBP-preserving interpolation propertie. Current values obtained from
! minimization of the ep2^T*Hp*ep2 + eu2^T*Hu*eu2 functional, 
! where eu2 and ep2 are the 2nd order polynomial interpolation error vectors
! i.e. eu2 = Ip2u * xp^2 - xu^2, ep2 = Iu2p * xu^2 - xp^2

real(kind=8), parameter :: c13 = 0.1264926207034815_8, &!10456777.0_8/180130106.0_8, 
                           c14 =-0.029892648949074313_8 !-89562251.0_8/1080780636.0_8

real(kind=8), parameter :: q11 = 1.0_8/2.0_8 + c13 + 2.0_8*c14, &
                           q12 = 1.0_8/2.0_8 - 2.0_8*c13 -3.0_8*c14, &
                           q13 = c13, &
                           q14 = c14, &
                           q21 = -8.0_8/63.0_8 - 52.0_8/21.0_8*c13 - 104.0_8/21.0_8*c14, &
                           q22 = 29.0_8/42.0_8 + 104.0_8/21.0_8*c13 + 52.0_8/7.0_8*c14, &
                           q23 = -52.0_8/21.0_8*c13 + 1.0_8/2.0_8, &
                           q24 = -52.0_8/21.0_8*c14 - 4.0_8/63.0_8, &
                           q31 = 26.0_8/25.0_8*c13 + 52.0_8/25.0_8*c14 - 1.0_8/25.0_8, &
                           q32 = -1.0_8/50.0_8 - 52.0_8/25.0_8*c13 - 78.0_8/25.0_8*c14, &
                           q33 = 3.0_8/5.0_8 + 26.0_8/25.0_8*c13, &
                           q34 = 13.0_8/25.0_8 + 26.0_8/25.0_8*c14, &
                           q35 = -3.0_8/50.0_8


!Optimized vector->scalar x^2 interpolation
! real(kind=8), parameter :: q11 =  0.4892964885611128_8,   q12 = 0.48123107353632727_8, &
!                            q13 =  0.0696483872440123_8,   q14 = -0.040175949341450676_8,  &
!                            q21 = -0.07913099957454726_8,  q22 = 0.6127279550862234_8, &
!                            q23 =  0.5119370885511979_8,   q24 = -0.045534044062872786_8, &
!                            q31 = -0.0869982983587545_8,   q32 = 0.1691760860218048_8, &
!                            q33 =  0.3626427230326682_8,   q34 = 0.6151794893042863_8, &
!                            q35 = -0.06_8, &
!                            q41 =  0.018680545032460714_8, q42 = -0.10869571330417979_8, &
!                            q43 =  0.09884979151099299_8,  q44 = 0.4911653767607321_8, &
!                            q45 =  0.5625_8, q46 = -0.0625_8

!Optimized vector <-> scalar interp in x^2
! real(kind=8), parameter :: q11 =  0.5396575759839645_8,     q12 = 0.4475591035840635_8, &
!                            q13 = -0.014090935120015191_8,   q14 = 0.026874255551989024_8,  &
!                            q21 = -0.13589939229291748_8,    q22 = 0.6508224599464276_8, &
!                            q23 = 0.6060532569859003_8,      q24 = -0.12097632463940905_8, &
!                            q31 = -0.14875399135119022_8,    q32 = 0.2102338641077255_8, &
!                            q33 = 0.46579424583813406_8,     q34 = 0.5327258814053357_8, &
!                            q35 = -0.06_8, &
!                            q41 = 0.07812389082006573_8,     q42 = -0.14831895644807316_8, &
!                            q43 = -0.00023375956403533802_8, q44 = 0.5704288251920487_8, &
!                            q45 = 0.5625_8, q46 = -0.0625_8

type, public, extends(sbp_interp_t) :: sbp_interp_i2c_42_t
contains
    procedure, public :: apply_tile
end type sbp_interp_i2c_42_t

contains

subroutine apply_tile(this, f_out, f_in, mesh, direction)
    class(sbp_interp_i2c_42_t), intent(in)    :: this
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
                    f_out%p(1,j,k) = q11*f_in%p(1,j,k)+q12*f_in%p(2,j,k)+q13*f_in%p(3,j,k)+q14*f_in%p(4,j,k)
                end if
                if (mesh%is<=2 .and. mesh%ie>=2) then
                    f_out%p(2,j,k) = q21*f_in%p(1,j,k)+q22*f_in%p(2,j,k)+q23*f_in%p(3,j,k)+q24*f_in%p(4,j,k)
                end if
                if (mesh%is<=3 .and. mesh%ie>=3) then
                    f_out%p(3,j,k) = q31*f_in%p(1,j,k)+q32*f_in%p(2,j,k)+q33*f_in%p(3,j,k)+q34*f_in%p(4,j,k)+q35*f_in%p(5,j,k)
                end if
                do i = max(mesh%is, 4), min(mesh%ie, mesh%nx-3)
                    f_out%p(i,j,k) = (9.0_8*(f_in%p(i+1,j,k) + f_in%p(i  ,j,k)) &
                                           -(f_in%p(i+2,j,k) + f_in%p(i-1,j,k))) / 16.0_8
                end do
                if (mesh%is<=mesh%nx-2 .and. mesh%ie>=mesh%nx-2) then
                    f_out%p(mesh%nx-2,j,k) = q31*f_in%p(mesh%nx+1,j,k)+q32*f_in%p(mesh%nx,j,k)+q33*f_in%p(mesh%nx-1,j,k)+q34*f_in%p(mesh%nx-2,j,k)+q35*f_in%p(mesh%nx-3,j,k)
                end if
                if (mesh%is<=mesh%nx-1 .and. mesh%ie>=mesh%nx-1) then
                    f_out%p(mesh%nx-1,j,k) = q21*f_in%p(mesh%nx+1,j,k)+q22*f_in%p(mesh%nx,j,k)+q23*f_in%p(mesh%nx-1,j,k)+q24*f_in%p(mesh%nx-2,j,k)
                end if
                if (mesh%is<=mesh%nx .and. mesh%ie>=mesh%nx) then
                    f_out%p(mesh%nx,j,k) = q11*f_in%p(mesh%nx+1,j,k)+q12*f_in%p(mesh%nx,j,k)+q13*f_in%p(mesh%nx-1,j,k)+q14*f_in%p(mesh%nx-2,j,k)
                end if
            end do
        end do
    case("y")
        do k = mesh%ks, mesh%ke
            if (mesh%js<=1 .and. mesh%je>=1) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,1,k) = q11*f_in%p(i,1,k)+q12*f_in%p(i,2,k)+q13*f_in%p(i,3,k)+q14*f_in%p(i,4,k)
                end do
            end if
            if (mesh%js<=2 .and. mesh%je>=2) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,2,k) = q21*f_in%p(i,1,k)+q22*f_in%p(i,2,k)+q23*f_in%p(i,3,k)+q24*f_in%p(i,4,k)
                end do
            end if
            if (mesh%js<=3 .and. mesh%je>=3) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,3,k) = q31*f_in%p(i,1,k)+q32*f_in%p(i,2,k)+q33*f_in%p(i,3,k)+q34*f_in%p(i,4,k)+q35*f_in%p(i,5,k)
                end do
            end if
            do j = max(mesh%js, 4), min(mesh%je, mesh%ny-3)
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,k) = (9.0_8*(f_in%p(i,j+1,k) + f_in%p(i,j  ,k)) &
                                           -(f_in%p(i,j+2,k) + f_in%p(i,j-1,k))) / 16.0_8
                end do
            end do
            if (mesh%js<=mesh%ny-2 .and. mesh%je>=mesh%ny-2) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-2,k) = q31*f_in%p(i,mesh%ny+1,k)+q32*f_in%p(i,mesh%ny,k)+q33*f_in%p(i,mesh%ny-1,k)+q34*f_in%p(i,mesh%ny-2,k)+q35*f_in%p(i,mesh%ny-3,k)
                end do
            end if
            if (mesh%js<=mesh%ny-1 .and. mesh%je>=mesh%ny-1) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-1,k) = q21*f_in%p(i,mesh%ny+1,k)+q22*f_in%p(i,mesh%ny,k)+q23*f_in%p(i,mesh%ny-1,k)+q24*f_in%p(i,mesh%ny-2,k)
                end do
            end if
            if (mesh%js<=mesh%ny .and. mesh%je>=mesh%ny) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny,k) = q11*f_in%p(i,mesh%ny+1,k)+q12*f_in%p(i,mesh%ny,k)+q13*f_in%p(i,mesh%ny-1,k)+q14*f_in%p(i,mesh%ny-2,k)
                end do
            end if
        end do
    case("z")
        if (mesh%ks<=1 .and. mesh%ke>=1) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,1) = q11*f_in%p(i,j,1)+q12*f_in%p(i,j,2)+q13*f_in%p(i,j,3)+q14*f_in%p(i,j,4)
                end do
            end do
        end if
        if (mesh%ks<=2 .and. mesh%ke>=2) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,2) = q21*f_in%p(i,j,1)+q22*f_in%p(i,j,2)+q23*f_in%p(i,j,3)+q24*f_in%p(i,j,4)
                end do
            end do
        end if
        if (mesh%ks<=3 .and. mesh%ke>=3) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,3) = q31*f_in%p(i,j,1)+q32*f_in%p(i,j,2)+q33*f_in%p(i,j,3)+q34*f_in%p(i,j,4)+q35*f_in%p(i,j,5)
                end do
            end do
        end if
        do k = max(mesh%ks, 4), min(mesh%ke, mesh%nz-3)
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,k) = (9.0_8*(f_in%p(i,j,k+1) + f_in%p(i,j,k  )) &
                                           -(f_in%p(i,j,k+2) + f_in%p(i,j,k-1))) / 16.0_8
                end do
            end do
        end do
        if (mesh%ks<=mesh%nz-2 .and. mesh%ke>=mesh%nz-2) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-2) = q31*f_in%p(i,j,mesh%nz+1)+q32*f_in%p(i,j,mesh%nz)+q33*f_in%p(i,j,mesh%nz-1)+q34*f_in%p(i,j,mesh%nz-2)+q35*f_in%p(i,j,mesh%nz-3)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz-1 .and. mesh%ke>=mesh%nz-1) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-1) = q21*f_in%p(i,j,mesh%nz+1)+q22*f_in%p(i,j,mesh%nz)+q23*f_in%p(i,j,mesh%nz-1)+q24*f_in%p(i,j,mesh%nz-2)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz .and. mesh%ke>=mesh%nz) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz) = q11*f_in%p(i,j,mesh%nz+1)+q12*f_in%p(i,j,mesh%nz)+q13*f_in%p(i,j,mesh%nz-1)+q14*f_in%p(i,j,mesh%nz-2)
                end do
            end do
        end if
    case default
        call parcomm_global%abort("Wrong direction in sbp_interp_i2c_42_mod")
    end select

end subroutine apply_tile

end module sbp_interp_i2c_42_mod
