module coriolis_C_contra_lincons_mod

!Coriolis formulation for C-grid contravariant components
!conserving energy for linear SWE.
!
!The formulation is fuv = Qi*inv(J)*C*uv,
!ie. covariant components of Coriolis force 
!are first computed (and multiplied by J)
!then divided by J in velocity points and finaly 
!Qi is applied to transform to contravariant components

!C = (        0         | Ih2u Jh^2*f Iv2h)
!    (-Ih2v Jh^2*f Iu2h |      0          )

!Energy conservation proof relies on hypothetical inverse of Qi:
! uv^T*inv(Qi)^T*J*C*uv = 0

use grid_field_mod,                only : grid_field_t, tile_field_t
use domain_mod,                    only : domain_t
use abstract_coriolis_mod,         only : coriolis_operator_t
use mesh_mod,                      only : mesh_t, tile_mesh_t
use abstract_interpolators2d_mod,  only : interpolator2d_vec2vec_t
use abstract_co2contra_mod,        only : co2contra_operator_t
use vec_math_mod,                  only : divide_by_J_self

type, public, extends(coriolis_operator_t) :: coriolis_C_contra_lincons_t

    type(grid_field_t) :: fJ2  ! coriolis parameter multiplied by J**2
    type(grid_field_t) :: u_at_f, v_at_f ! input u and v interpolated to f-points and then u and v cov tendencies
    type(grid_field_t) :: fu, fv !covariant tends at u,v-points

    type(mesh_t), pointer :: f_mesh !mesh at f-points

    class(co2contra_operator_t), allocatable     :: co2contra
    class(interpolator2d_vec2vec_t), allocatable :: interp_v2f_op, interp_f2v_op

contains
    procedure, public :: calc_coriolis
    procedure, public :: calc_coriolis_vec_inv
    procedure, public :: calc_coriolis_contra

end type coriolis_C_contra_lincons_t

contains

subroutine calc_coriolis_contra(this, cor_ut, cor_vt, ut, vt, domain)
    class(coriolis_C_contra_lincons_t), intent(inout) :: this
    type(domain_t),                     intent(in)    :: domain
    type(grid_field_t),                 intent(inout) :: ut, vt !contravariant components
    type(grid_field_t),                 intent(inout) :: cor_ut, cor_vt !contravariant tentedcies

    integer(kind=4) :: t

    call this%interp_v2f_op%interp2d_vec2vec(this%u_at_f, this%v_at_f, ut, vt, domain)

    do t = this%f_mesh%ts, this%f_mesh%te
        call calc_coriolis_fpoints_tile(this%u_at_f%tile(t), this%v_at_f%tile(t), &
                                        this%fJ2%tile(t), this%f_mesh%tile(t))
    end do

    call this%interp_f2v_op%interp2d_vec2vec(this%fu, this%fv, this%u_at_f, this%v_at_f, domain)

    call divide_by_J_self(this%fu, domain%mesh_u)
    call divide_by_J_self(this%fv, domain%mesh_v)

    call this%co2contra%transform(cor_ut, cor_vt, this%fu, this%fv, domain)

end subroutine calc_coriolis_contra

subroutine calc_coriolis_fpoints_tile(u, v, f, mesh)

    use mesh_mod, only : tile_mesh_t

    type(tile_field_t),   intent(inout) :: u, v
    type(tile_field_t),   intent(in)    :: f
    type(tile_mesh_t),    intent(in)    :: mesh

    integer(kind=4) :: i, j, k
    real(kind=8)    :: ut_tend, vt_tend

    do k=mesh%ks, mesh%ke
        do j=mesh%js, mesh%je
            do i=mesh%is, mesh%ie

                ut_tend = v%p(i,j,k)*f%p(i,j,1)
                vt_tend =-u%p(i,j,k)*f%p(i,j,1)

                u%p(i,j,k) = ut_tend
                v%p(i,j,k) = vt_tend

            end do
        end do
    end do

end subroutine calc_coriolis_fpoints_tile

subroutine calc_coriolis(this, cor_u, cor_v, ut, vt, domain)
    class(coriolis_C_contra_lincons_t), intent(inout) :: this
    type(domain_t),                     intent(in)    :: domain
    type(grid_field_t),                 intent(inout) :: ut, vt!contravariant components
    type(grid_field_t),                 intent(inout) :: cor_u, cor_v

    integer(kind=4) :: t

    call domain%parcomm%abort("calc_coriolis is not implemented!")

end subroutine calc_coriolis

subroutine calc_coriolis_vec_inv(this, cor_u, cor_v, hu, hv, h, curl, domain)
    class(coriolis_C_contra_lincons_t), intent(inout) :: this
    type(domain_t),                     intent(in)    :: domain
    type(grid_field_t),                 intent(inout) :: hu, hv! massflux contravariant components
    type(grid_field_t),                 intent(inout) :: h, curl
    type(grid_field_t),                 intent(inout) :: cor_u, cor_v

    call domain%parcomm%abort("calc_coriolis_vec_inv is not implemented!")

end subroutine calc_coriolis_vec_inv

end module coriolis_C_contra_lincons_mod