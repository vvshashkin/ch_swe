module baroclinic_instability_test_parameters_mod

    use const_mod,  only : pi

    implicit none

    real(kind=8), parameter :: p0 = 1e5_8
    real(kind=8), parameter :: sigma_0 = 0.252_8 !jet core level pressure 252hPa
    real(kind=8), parameter :: sigma_t = 0.2_8   !tropopause level 200hPa
    real(kind=8), parameter :: u0 = 35.0_8
    real(kind=8), parameter :: Gamma = 0.005_8 !Temperature lapse rate K/m
    real(kind=8), parameter :: T0 = 288.0_8, delta_T = 4.8e5
    !perturbation position:
    real(kind=8), parameter :: lam_c = pi/9.0_8, phi_c = 2.0_8*pi/9.0_8
    real(kind=8), parameter :: xc = cos(phi_c)*cos(lam_c), yc = cos(phi_c)*sin(lam_c), zc = sin(phi_c)

end module
