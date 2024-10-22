module metric_mod

use mesh_mod,      only : mesh_t
use orography_mod, only : orography_t
use parcomm_mod,   only : parcomm_global

implicit none

!Abstract type to aquire horizontal grid characteristics from panel coordinates
type, abstract :: metric_t
    real(kind=8) :: scale      ! projection scale factor (radius of the sphere in the case of equiangular cubed sphere)
    real(kind=8) :: vertical_scale ! vertical coordinate scaling factor
    real(kind=8) :: alpha0, beta0 ! lower bound of panel coordinates
    real(kind=8) :: alpha1, beta1 ! upper bound of panel coordinates

    real(kind=8) :: omega ! rotation speed of the coord system
    logical      :: coriolis_const = .false.
    real(kind=8), allocatable :: rotation_axis(:)!direction of the angular velocity vector

    real(kind=8), allocatable :: rotation_matrix(:,:)

contains

    !New interfaces
    procedure(calculate_r_orog),  deferred :: calculate_r_orog
    procedure(calculate_r_2d),    deferred :: calculate_r_2d !no orography case
    generic :: calculate_r => calculate_r_orog, calculate_r_2d

    procedure(calculate_h),       deferred :: calculate_h

    procedure(calculate_vec_cov_orog), deferred :: calculate_a1_orog
    procedure(calculate_vec_cov_2d),   deferred :: calculate_a1_2d
    generic :: calculate_a1 => calculate_a1_orog, calculate_a1_2d

    procedure(calculate_vec_cov_orog), deferred :: calculate_a2_orog
    procedure(calculate_vec_cov_2d),   deferred :: calculate_a2_2d
    generic :: calculate_a2 => calculate_a2_orog, calculate_a2_2d

    procedure(calculate_vec4_orog),    deferred :: calculate_a3_orog
    procedure(calculate_vec_cov_2d),   deferred :: calculate_a3_2d
    generic :: calculate_a3 => calculate_a3_orog, calculate_a3_2d

    procedure(calculate_vec_contra_orog), deferred :: calculate_b1_orog
    procedure(calculate_vec_contra_2d),   deferred :: calculate_b1_2d
    generic :: calculate_b1 => calculate_b1_orog, calculate_b1_2d

    procedure(calculate_vec_contra_orog), deferred :: calculate_b2_orog
    procedure(calculate_vec_contra_2d),   deferred :: calculate_b2_2d
    generic :: calculate_b2 => calculate_b2_orog, calculate_b2_2d

    procedure(calculate_vec_contra4_orog), deferred :: calculate_b3_orog
    procedure(calculate_vec_contra4_2d),   deferred :: calculate_b3_2d
    generic :: calculate_b3 => calculate_b3_orog, calculate_b3_2d

    procedure(calculate_metric_tensor_orog), deferred :: calculate_Q_orog
    procedure(calculate_metric_tensor_2d),   deferred :: calculate_Q_2d
    generic :: calculate_Q => calculate_Q_orog, calculate_Q_2d

    procedure(calculate_metric_tensor_orog), deferred :: calculate_Qi_orog
    procedure(calculate_metric_tensor_2d),   deferred :: calculate_Qi_2d
    generic :: calculate_Qi => calculate_Qi_orog, calculate_Qi_2d

    procedure(calculate_jacobian_orog), deferred :: calculate_J_orog
    procedure(calculate_jacobian_2d),   deferred :: calculate_J_2d
    generic :: calculate_J => calculate_J_orog, calculate_J_2d

    procedure(calculate_Christoffel_orog), deferred :: calculate_G_orog
    procedure(calculate_Christoffel_2d), deferred   :: calculate_G_2d
    generic :: calculate_G => calculate_G_orog, calculate_G_2d

    procedure(calc_Christoffel1_i), deferred :: calculate_S

    procedure(cart_to_native_transform), deferred :: transform_cartesian_to_native

    ! subroutine interfaces
    procedure :: set_curvilinear_mesh
    procedure :: transform_xyz_to_native

end type metric_t

abstract interface

    pure function calculate_r_orog(this, panel_ind, alpha, beta, eta, h_surf, h_top) result(r)
        import metric_t
        class(metric_t), intent(in) :: this
        integer(kind=4), intent(in) :: panel_ind
        real(kind=8),    intent(in) :: alpha, beta, eta
        real(kind=8),    intent(in) :: h_surf, h_top
        real(kind=8)                :: r(3)
    end function calculate_r_orog
    pure function calculate_r_2d(this, panel_ind, alpha, beta) result(r)
        import metric_t
        class(metric_t), intent(in) :: this
        integer(kind=4), intent(in) :: panel_ind
        real(kind=8),    intent(in) :: alpha, beta
        real(kind=8)                :: r(3)
    end function calculate_r_2d
    pure function calculate_h(this,panel_ind,alpha,beta,eta,h_surf,h_top) result(h)
        import metric_t
        class(metric_t), intent(in) :: this
        integer(kind=4), intent(in) :: panel_ind
        real(kind=8),    intent(in) :: alpha, beta, eta
        real(kind=8),    intent(in) :: h_surf, h_top
        real(kind=8)                :: h
    end function calculate_h
    pure function calculate_vec_cov_orog(this, panel_ind, alpha, beta, eta, &
                                         h_surf, dcov_h_surf, h_top) result(a)
        import metric_t
        class(metric_t), intent(in) :: this
        integer(kind=4), intent(in) :: panel_ind
        real(kind=8),    intent(in) :: alpha, beta, eta
        real(kind=8),    intent(in) :: h_surf, dcov_h_surf, h_top
        real(kind=8)                :: a(4)
    end function calculate_vec_cov_orog
    pure function calculate_vec4_orog(this, panel_ind, alpha, beta, eta, &
                                         h_surf,h_top) result(a)
        import metric_t
        class(metric_t), intent(in) :: this
        integer(kind=4), intent(in) :: panel_ind
        real(kind=8),    intent(in) :: alpha, beta, eta
        real(kind=8),    intent(in) :: h_surf, h_top
        real(kind=8)                :: a(4)
    end function calculate_vec4_orog
    pure function calculate_vec_cov_2d(this, panel_ind, alpha, beta) result(a)
        import metric_t
        class(metric_t), intent(in) :: this
        integer(kind=4), intent(in) :: panel_ind
        real(kind=8),    intent(in) :: alpha, beta
        real(kind=8)                :: a(4)
    end function calculate_vec_cov_2d

    pure function calculate_vec_contra_orog(this, panel_ind, alpha, beta, eta, &
                                            h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top) result(b)
        import metric_t
        class(metric_t), intent(in) :: this
        integer(kind=4), intent(in) :: panel_ind
        real(kind=8),    intent(in) :: alpha, beta, eta
        real(kind=8),    intent(in) :: h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top
        real(kind=8)                :: b(3)
    end function calculate_vec_contra_orog
    pure function calculate_vec_contra_2d(this, panel_ind, alpha, beta) result(b)
        import metric_t
        class(metric_t), intent(in) :: this
        integer(kind=4), intent(in) :: panel_ind
        real(kind=8),    intent(in) :: alpha, beta
        real(kind=8)                :: b(3)
    end function calculate_vec_contra_2d

    pure function calculate_vec_contra4_orog(this, panel_ind, alpha, beta, eta, &
                                            h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top) result(b)
        import metric_t
        class(metric_t), intent(in) :: this
        integer(kind=4), intent(in) :: panel_ind
        real(kind=8),    intent(in) :: alpha, beta, eta
        real(kind=8),    intent(in) :: h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top
        real(kind=8)                :: b(4)
    end function calculate_vec_contra4_orog
    pure function calculate_vec_contra4_2d(this, panel_ind, alpha, beta) result(b)
        import metric_t
        class(metric_t), intent(in) :: this
        integer(kind=4), intent(in) :: panel_ind
        real(kind=8),    intent(in) :: alpha, beta
        real(kind=8)                :: b(4)
    end function calculate_vec_contra4_2d

    pure function calculate_metric_tensor_orog(this, panel_ind, alpha, beta, eta, &
                                            h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top) result(Q)
        import metric_t
        class(metric_t), intent(in) :: this
        integer(kind=4), intent(in) :: panel_ind
        real(kind=8),    intent(in) :: alpha, beta, eta
        real(kind=8),    intent(in) :: h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top
        real(kind=8)                :: Q(6)
    end function calculate_metric_tensor_orog
    pure function calculate_metric_tensor_2d(this, panel_ind, alpha, beta) result(Q)
        import metric_t
        class(metric_t), intent(in) :: this
        integer(kind=4), intent(in) :: panel_ind
        real(kind=8),    intent(in) :: alpha, beta
        real(kind=8)                :: Q(6)
    end function calculate_metric_tensor_2d

    pure function calculate_jacobian_orog(this, panel_ind, alpha, beta, eta, h_surf, h_top) result(J)
        import metric_t
        class(metric_t), intent(in) :: this
        integer(kind=4), intent(in) :: panel_ind
        real(kind=8),    intent(in) :: alpha, beta, eta
        real(kind=8),    intent(in) :: h_surf, h_top
        real(kind=8)                :: J
    end function calculate_jacobian_orog
    pure function calculate_jacobian_2d(this, panel_ind, alpha, beta) result(J)
        import metric_t
        class(metric_t), intent(in) :: this
        integer(kind=4), intent(in) :: panel_ind
        real(kind=8),    intent(in) :: alpha, beta
        real(kind=8)                :: J
    end function calculate_jacobian_2d

    pure function calculate_Christoffel_orog(this, panel_ind, alpha, beta, eta, &
                                            h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top) result(G)
        import metric_t
        class(metric_t), intent(in) :: this
        integer(kind=4), intent(in) :: panel_ind
        real(kind=8),    intent(in) :: alpha, beta, eta
        real(kind=8),    intent(in) :: h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top
        real(kind=8)                :: G(3,3,3)
    end function calculate_Christoffel_orog
    pure function calculate_Christoffel_2d(this, panel_ind, alpha, beta) result(G)
        import metric_t
        class(metric_t), intent(in) :: this
        integer(kind=4), intent(in) :: panel_ind
        real(kind=8),    intent(in) :: alpha, beta
        real(kind=8)                :: G(3,3,3)
    end function calculate_Christoffel_2d

    pure function calc_Christoffel1_i(this, panel_ind, alpha, beta) result(S)
        import metric_t
        class(metric_t),    intent(in) :: this
        integer(kind=4),    intent(in) :: panel_ind
        real(kind=8),       intent(in) :: alpha, beta
        real(kind=8)                   :: S(3,2)
    end function calc_Christoffel1_i

    subroutine cart_to_native_transform(this,panel_ind, alpha, beta, r)
        import metric_t
        class(metric_t), intent(in)  :: this
        integer(kind=4), intent(out) :: panel_ind
        real(kind=8),    intent(out) :: alpha, beta
        real(kind=8),    intent(in)  :: r(3)
    end subroutine cart_to_native_transform

end interface

contains

subroutine set_curvilinear_mesh(this,mesh,orog)
    class(metric_t),        intent(in)    :: this
    type(mesh_t),           intent(inout) :: mesh
    type(orography_t),      intent(in), optional :: orog

    call parcomm_global%abort("set_tile_mesh is not implemented for this metric class")
end subroutine

subroutine transform_xyz_to_native(this,alpha,beta, panel_ind,x,y,z,mesh)
    use grid_field_mod, only : grid_field_t, tile_field_t
    use mesh_mod,       only : mesh_t

    class(metric_t),    intent(in)    :: this
    type(grid_field_t), intent(inout) :: alpha, beta, panel_ind
    type(grid_field_t), intent(in)    :: x, y, z 
    type(mesh_t),       intent(in)    :: mesh

    real(kind=8)    :: r(3)
    integer(kind=4) :: t, i, j, k, p_ind

    do t = mesh%ts, mesh%te
        do k = mesh%tile(t)%ks, mesh%tile(t)%ke
            do j = mesh%tile(t)%js, mesh%tile(t)%je
                do i = mesh%tile(t)%is, mesh%tile(t)%ie
                    r(1:3) = [x%tile(t)%p(i,j,k), y%tile(t)%p(i,j,k), z%tile(t)%p(i,j,k)]
                    call this%transform_cartesian_to_native(p_ind,alpha%tile(t)%p(i,j,k),beta%tile(t)%p(i,j,k), r)
                    panel_ind%tile(t)%p(i,j,k) = p_ind
                end do
            end do
        end do
    end do

end subroutine

end module metric_mod
