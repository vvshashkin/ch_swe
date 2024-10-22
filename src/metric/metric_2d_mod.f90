module metric_2d_mod

use mesh_mod,      only : mesh_t
use parcomm_mod,   only : parcomm_global

implicit none

! This abstract class implements 2d general curvilinear
! metric at the logicaly paneled domains

type, public, abstract :: metric_2d_t
    real(kind=8) :: scale
    real(kind=8) :: alpha0, beta0 ! lower bound of panel coordinates
    real(kind=8) :: alpha1, beta1 ! upper bound of panel coordinates
    real(kind=8) :: omega ! rotation speed of the coord system

    real(kind=8), allocatable :: rotation_axis(:) !direction of the angular velocity vector
    real(kind=8), allocatable :: rotation_matrix(:,:)
contains
    procedure(calc_metric_vec_quantity_i),    deferred :: calc_r
    procedure(calc_metric_vec_quantity_i),    deferred :: calc_a1, calc_a2
    procedure(calc_metric_vec_quantity_i),    deferred :: calc_b1, calc_b2
    procedure(calc_metric_tensor_i),          deferred :: calc_Q, calc_Qi 
    procedure(calc_J_i),                      deferred :: calc_J
    procedure(calc_Christoffel_i),            deferred :: calc_G
    procedure(calc_Christoffel1_i),           deferred :: calc_S
    procedure(transform_cart2native_i),       deferred :: transform_cart2native
    procedure(set_curvilinear_mesh_i),        deferred :: set_curvilinear_mesh

    procedure :: transform_xyz_to_native
    procedure :: calc_cartesian_hor_wind
end type metric_2d_t

abstract interface
    pure function calc_J_i(this, panel_ind, alpha, beta) result(J)
        import metric_2d_t
        class(metric_2d_t), intent(in) :: this
        integer(kind=4),    intent(in) :: panel_ind
        real(kind=8),       intent(in) :: alpha, beta
        real(kind=8)                   :: J
    end function calc_J_i
    pure function calc_metric_vec_quantity_i(this, panel_ind, alpha, beta) result(vec)
        import metric_2d_t
        class(metric_2d_t), intent(in) :: this
        integer(kind=4),    intent(in) :: panel_ind
        real(kind=8),       intent(in) :: alpha, beta
        real(kind=8)                   :: vec(3)
    end function calc_metric_vec_quantity_i
    pure function calc_metric_tensor_i(this, panel_ind, alpha, beta) result(tensor)
        import metric_2d_t
        class(metric_2d_t), intent(in) :: this
        integer(kind=4),    intent(in) :: panel_ind
        real(kind=8),       intent(in) :: alpha, beta
        real(kind=8)                   :: tensor(3)
    end function calc_metric_tensor_i
    pure function calc_Christoffel_i(this, panel_ind, alpha, beta) result(G)
        import metric_2d_t
        class(metric_2d_t), intent(in) :: this
        integer(kind=4),    intent(in) :: panel_ind
        real(kind=8),       intent(in) :: alpha, beta
        real(kind=8)                   :: G(2, 2, 2)
    end function calc_Christoffel_i
    pure function calc_Christoffel1_i(this, panel_ind, alpha, beta) result(S)
        import metric_2d_t
        class(metric_2d_t), intent(in) :: this
        integer(kind=4),    intent(in) :: panel_ind
        real(kind=8),       intent(in) :: alpha, beta
        real(kind=8)                   :: S(3,2)
    end function calc_Christoffel1_i
    subroutine transform_cart2native_i(this, panel_ind, alpha, beta, r)
        import metric_2d_t
        class(metric_2d_t), intent(in)  :: this
        integer(kind=4),    intent(out) :: panel_ind
        real(kind=8),       intent(out) :: alpha, beta
        real(kind=8),       intent(in)  :: r(3)
    end subroutine transform_cart2native_i
    subroutine set_curvilinear_mesh_i(this, mesh)
        import metric_2d_t, mesh_t
        class(metric_2d_t), intent(in)    :: this
        type(mesh_t),       intent(inout) :: mesh
    end subroutine
end interface 

contains

subroutine transform_xyz_to_native(this,alpha,beta, panel_ind,x,y,z,mesh)
    use grid_field_mod, only : grid_field_t, tile_field_t
    use mesh_mod,       only : mesh_t

    class(metric_2d_t), intent(in)    :: this
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
                    call this%transform_cart2native(p_ind,alpha%tile(t)%p(i,j,k),beta%tile(t)%p(i,j,k), r)
                    panel_ind%tile(t)%p(i,j,k) = p_ind
                end do
            end do
        end do
    end do

end subroutine

subroutine calc_cartesian_hor_wind(this, vx, vy, vz, u, v, &
                                   panel_ind, alpha, beta, &
                                   mesh, native_components_type)

    use grid_field_mod, only : grid_field_t, tile_field_t
    use mesh_mod,       only : mesh_t

    class(metric_2d_t),   intent(in)    :: this
    type(grid_field_t),   intent(inout) :: vx, vy, vz
    type(grid_field_t),   intent(in)    :: u, v, alpha, beta, panel_ind
    type(mesh_t),         intent(in)    :: mesh
    character(len=*),     intent(in)    :: native_components_type

    integer(kind=4) :: i, j, k, t
    real(kind=8)    :: a1(3), a2(3), r(3)

    select case (native_components_type)
    case ("contravariant")

    do t = mesh%ts, mesh%te
        do k = mesh%tile(t)%ks, mesh%tile(t)%ke
            do j = mesh%tile(t)%js, mesh%tile(t)%je
                do i = mesh%tile(t)%is, mesh%tile(t)%ie

                    a1(1:3) = this%calc_a1(int(panel_ind%tile(t)%p(i,j,k)), alpha%tile(t)%p(i,j,k), beta%tile(t)%p(i,j,k))
                    a2(1:3) = this%calc_a2(int(panel_ind%tile(t)%p(i,j,k)), alpha%tile(t)%p(i,j,k), beta%tile(t)%p(i,j,k))

                    vx%tile(t)%p(i,j,k) = u%tile(t)%p(i,j,k)*a1(1)+v%tile(t)%p(i,j,k)*a2(1)
                    vy%tile(t)%p(i,j,k) = u%tile(t)%p(i,j,k)*a1(2)+v%tile(t)%p(i,j,k)*a2(2)
                    vz%tile(t)%p(i,j,k) = u%tile(t)%p(i,j,k)*a1(3)+v%tile(t)%p(i,j,k)*a2(3)

                end do
            end do
        end do
    end do

    case ("covariant")
        call parcomm_global%abort("shallow atmosphere metric, calc cartesian hor wind, covariant native components not implemented")
    case default
        call parcomm_global%abort("shallow atmosphere metric, calc cartesian hor wind, unknown native components type: "//native_components_type)
    end select

end subroutine

end module metric_2d_mod