module sbp_interp_c2i_42_mod

use sbp_interp_mod,   only : sbp_interp_t
use grid_field_mod, only : grid_field_t, tile_field_t
use mesh_mod,       only : tile_mesh_t, mesh_t
use tile_mod,       only : tile_t
use parcomm_mod,    only : parcomm_global

use sbp_interp_i2c_42_mod, only : c13, c14

implicit none

real(kind=8), parameter :: q11 = 39.0_8/28.0_8+39.0_8/14.0_8*c13+39.0_8/7.0_8*c14, &
                           q12 = -2.0_8/7.0_8-39.0_8/7.0_8*c13-78.0_8/7.0_8*c14, &
                           q13 = 39.0_8/14.0_8*c13+39.0_8/7.0_8*c14-3.0_8/28.0_8, &
                           q21 = 13.0_8/27.0_8-52.0_8/27.0_8*c13-26.0_8/9.0_8*c14, &
                           q22 = 29.0_8/54.0_8+104.0_8/27.0_8*c13+52.0_8/9.0_8*c14, &
                           q23 = -1.0_8/54.0_8-52.0_8/27.0_8*c13-26.0_8/9.0_8*c14, &
                           q31 = 13.0_8/12.0_8*c13, &
                           q32 = -13.0_8/6.0_8*c13+7.0_8/16.0_8, &
                           q33 = 5.0_8/8.0_8+13.0_8/12.0_8*c13, &
                           q34 = -1.0_8/16.0_8, &
                           q41 = 78.0_8/71.0_8*c14, &
                           q42 = -156.0_8/71.0_8*c14-4.0_8/71.0_8, &
                           q43 =  39.0_8/71.0_8+78.0_8/71.0_8*c14, &
                           q44 = 81.0_8/142.0_8, & 
                           q45 = -9.0_8/142.0_8

!optimized only vector to scalar interp
! real(kind=8), parameter :: q11 =  1.3630402181345285_8,  q12 = -0.17804474904273135_8,  &
!                            q13 = -0.23303115631809243_8, q14 =  0.048035687226327554_8, &
!                            q21 =  0.46340770044238916_8, q22 =  0.4765661872892848_8,   &
!                            q23 =  0.15664452409426372_8, q24 = -0.09661841182593758_8,  &
!                            q31 =  0.07545241951434666_8, q32 =  0.44794495248229815_8,  &
!                            q33 =  0.37775283649236274_8, q34 = 0.09884979151099299_8,   &
!                            q41 = -0.04413695843145285_8, q42 = -0.04040344754874627_8,  &
!                            q43 =  0.649837488701711_8,   q44 =  0.4980831989686297_8,   &
!                            q45 = -0.06338028169014084_8
!optimized both vector to scalar and scalar to vector interp
! real(kind=8), parameter :: q11 =  1.5033318188124725_8,   q12 = -0.3057736326590643_8,     &
!                            q13 = -0.3984481911192596_8,   q14 =  0.2008900049658833_8,     &
!                            q21 =  0.43098284048835733_8,  q22 =  0.5061952466249993_8,     &
!                            q23 =  0.194660985284931_8,    q24 = -0.13183907239828724_8,    &
!                            q31 = -0.01526517971334979_8,  q32 =  0.5302965998626628_8,     &
!                            q33 =  0.485202339414723_8,    q34 = -0.00023375956403533802_8, &
!                            q41 =  0.029523830043030195_8, q42 = -0.1073451894687714_8,     &
!                            q43 =  0.5627386071183124_8,   q44 =  0.5784630339975705_8,     &
!                            q45 = -0.06338028169014084_8

type, public, extends(sbp_interp_t) :: sbp_interp_c2i_42_t
contains
    procedure, public :: apply_tile
end type sbp_interp_c2i_42_t

contains

subroutine apply_tile(this, f_out, f_in, mesh, direction)
    class(sbp_interp_c2i_42_t), intent(in)    :: this
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
                    f_out%p(1,j,k) = q11*f_in%p(1,j,k)+q12*f_in%p(2,j,k)+q13*f_in%p(3,j,k)
                end if
                if (mesh%is<=2 .and. mesh%ie>=2) then
                    f_out%p(2,j,k) = q21*f_in%p(1,j,k)+q22*f_in%p(2,j,k)+q23*f_in%p(3,j,k)
                end if
                if (mesh%is<=3 .and. mesh%ie>=3) then
                    f_out%p(3,j,k) = q31*f_in%p(1,j,k)+q32*f_in%p(2,j,k)+q33*f_in%p(3,j,k)+q34*f_in%p(4,j,k)
                end if
                if (mesh%is<=4 .and. mesh%ie>=4) then
                    f_out%p(4,j,k) = q41*f_in%p(1,j,k)+q42*f_in%p(2,j,k)+q43*f_in%p(3,j,k)+q44*f_in%p(4,j,k)+q45*f_in%p(5,j,k)
                end if
                do i = max(mesh%is, 5), min(mesh%ie, mesh%nx-4)
                    f_out%p(i,j,k) = (9.0_8*(f_in%p(i  ,j,k) + f_in%p(i-1,j,k)) &
                                           -(f_in%p(i+1,j,k) + f_in%p(i-2,j,k)))/16.0_8
                end do
                if (mesh%is<=mesh%nx-3 .and. mesh%ie>=mesh%nx-3) then
                    f_out%p(mesh%nx-3,j,k) = q41*f_in%p(mesh%nx-1,j,k)+q42*f_in%p(mesh%nx-2,j,k)+q43*f_in%p(mesh%nx-3,j,k)+q44*f_in%p(mesh%nx-4,j,k)+q45*f_in%p(mesh%nx-5,j,k)
                end if
                if (mesh%is<=mesh%nx-2 .and. mesh%ie>=mesh%nx-2) then
                    f_out%p(mesh%nx-2,j,k) = q31*f_in%p(mesh%nx-1,j,k)+q32*f_in%p(mesh%nx-2,j,k)+q33*f_in%p(mesh%nx-3,j,k)+q34*f_in%p(mesh%nx-4,j,k)
                end if
                if (mesh%is<=mesh%nx-1 .and. mesh%ie>=mesh%nx-1) then
                    f_out%p(mesh%nx-1,j,k) = q21*f_in%p(mesh%nx-1,j,k)+q22*f_in%p(mesh%nx-2,j,k)+q23*f_in%p(mesh%nx-3,j,k)
                end if
                if (mesh%is<=mesh%nx .and. mesh%ie>=mesh%nx) then
                    f_out%p(mesh%nx,j,k) = q11*f_in%p(mesh%nx-1,j,k)+q12*f_in%p(mesh%nx-2,j,k)+q13*f_in%p(mesh%nx-3,j,k)
                end if
            end do
        end do
    case("y")
        do k = mesh%ks, mesh%ke
            if (mesh%js<=1 .and. mesh%je>=1) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,1,k) = q11*f_in%p(i,1,k)+q12*f_in%p(i,2,k)+q13*f_in%p(i,3,k)
                end do
            end if
            if (mesh%js<=2 .and. mesh%je>=2) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,2,k) = q21*f_in%p(i,1,k)+q22*f_in%p(i,2,k)+q23*f_in%p(i,3,k)
                end do
            end if
            if (mesh%js<=3 .and. mesh%je>=3) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,3,k) = q31*f_in%p(i,1,k)+q32*f_in%p(i,2,k)+q33*f_in%p(i,3,k)+q34*f_in%p(i,4,k)
                end do
            end if
            if (mesh%js<=4 .and. mesh%je>=4) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,4,k) = q41*f_in%p(i,1,k)+q42*f_in%p(i,2,k)+q43*f_in%p(i,3,k)+q44*f_in%p(i,4,k)+q45*f_in%p(i,5,k)
                end do
            end if
            do j = max(mesh%js, 5), min(mesh%je, mesh%ny-4)
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,k) = (9.0_8*(f_in%p(i,j  ,k) + f_in%p(i,j-1,k)) &
                                           -(f_in%p(i,j+1,k) + f_in%p(i,j-2,k)))/16.0_8
                end do
            end do
            if (mesh%js<=mesh%ny-3 .and. mesh%je>=mesh%ny-3) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-3,k) = q41*f_in%p(i,mesh%ny-1,k)+q42*f_in%p(i,mesh%ny-2,k)+q43*f_in%p(i,mesh%ny-3,k)+q44*f_in%p(i,mesh%ny-4,k)+q45*f_in%p(i,mesh%ny-5,k)
                end do
            end if
            if (mesh%js<=mesh%ny-2 .and. mesh%je>=mesh%ny-2) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-2,k) = q31*f_in%p(i,mesh%ny-1,k)+q32*f_in%p(i,mesh%ny-2,k)+q33*f_in%p(i,mesh%ny-3,k)+q34*f_in%p(i,mesh%ny-4,k)
                end do
            end if
            if (mesh%js<=mesh%ny-1 .and. mesh%je>=mesh%ny-1) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny-1,k) = q21*f_in%p(i,mesh%ny-1,k)+q22*f_in%p(i,mesh%ny-2,k)+q23*f_in%p(i,mesh%ny-3,k)
                end do
            end if
            if (mesh%js<=mesh%ny .and. mesh%je>=mesh%ny) then
                do i = mesh%is, mesh%ie
                    f_out%p(i,mesh%ny,k) = q11*f_in%p(i,mesh%ny-1,k)+q12*f_in%p(i,mesh%ny-2,k)+q13*f_in%p(i,mesh%ny-3,k)
                end do
            end if
        end do
    case("z")
        if (mesh%ks<=1 .and. mesh%ke>=1) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,1) = q11*f_in%p(i,j,1)+q12*f_in%p(i,j,2)+q13*f_in%p(i,j,3)
                end do
            end do
        end if
        if (mesh%ks<=2 .and. mesh%ke>=2) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,2) = q21*f_in%p(i,j,1)+q22*f_in%p(i,j,2)+q23*f_in%p(i,j,3)
                end do
            end do
        end if
        if (mesh%ks<=3 .and. mesh%ke>=3) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,3) = q31*f_in%p(i,j,1)+q32*f_in%p(i,j,2)+q33*f_in%p(i,j,3)+q34*f_in%p(i,j,4)
                end do
            end do
        end if
        if (mesh%ks<=4 .and. mesh%ke>=4) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,4) = q41*f_in%p(i,j,1)+q42*f_in%p(i,j,2)+q43*f_in%p(i,j,3)+q44*f_in%p(i,j,4)+q45*f_in%p(i,j,5)
                end do
            end do
        end if
        do k = max(mesh%ks, 5), min(mesh%ke, mesh%nz-4)
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,k) = (9.0_8*(f_in%p(i,j,k  ) + f_in%p(i,j,k-1)) &
                                           -(f_in%p(i,j,k+1) + f_in%p(i,j,k-2)))/16.0_8
                end do
            end do
        end do
        if (mesh%ks<=mesh%nz-3 .and. mesh%ke>=mesh%nz-3) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-3) = q41*f_in%p(i,j,mesh%nz-1)+q42*f_in%p(i,j,mesh%nz-2)+q43*f_in%p(i,j,mesh%nz-3)+q44*f_in%p(i,j,mesh%nz-4)+q45*f_in%p(i,j,mesh%nz-5)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz-2 .and. mesh%ke>=mesh%nz-2) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-2) = q31*f_in%p(i,j,mesh%nz-1)+q32*f_in%p(i,j,mesh%nz-2)+q33*f_in%p(i,j,mesh%nz-3)+q34*f_in%p(i,j,mesh%nz-4)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz-1 .and. mesh%ke>=mesh%nz-1) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz-1) = q21*f_in%p(i,j,mesh%nz-1)+q22*f_in%p(i,j,mesh%nz-2)+q23*f_in%p(i,j,mesh%nz-3)
                end do
            end do
        end if
        if (mesh%ks<=mesh%nz .and. mesh%ke>=mesh%nz) then
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    f_out%p(i,j,mesh%nz) = q11*f_in%p(i,j,mesh%nz-1)+q12*f_in%p(i,j,mesh%nz-2)+q13*f_in%p(i,j,mesh%nz-3)
                end do
            end do
        end if
    case default
        call parcomm_global%abort("Wrong direction in sbp_interp_c2i_42_mod")
    end select

end subroutine apply_tile

end module sbp_interp_c2i_42_mod
