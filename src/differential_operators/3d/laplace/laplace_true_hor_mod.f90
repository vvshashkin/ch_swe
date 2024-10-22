module laplace_true_hor_mod

use abstract_laplace_mod,      only : laplace_operator_t
use abstract_co2contra_3d_mod, only : co2contra_3d_operator_t
use abstract_div_3d_mod,       only : div_3d_operator_t
use abstract_grad_3d_mod,      only : grad_3d_operator_t
use grid_field_mod,            only : grid_field_t
use mesh_mod,                  only : mesh_t
use domain_mod,                only : domain_t

implicit none

type, extends(laplace_operator_t) :: laplace_true_hor_t
    class(grad_3d_operator_t),      allocatable :: grad
    class(div_3d_operator_t),       allocatable :: div
    class(co2contra_3d_operator_t), allocatable :: co2contra
    type(grid_field_t) :: fx, fy, fz, fxt, fyt, fzt, hor_factor

    contains
        procedure :: calc_laplace
end type

contains

subroutine calc_laplace(this, f1, f, domain)
    class(laplace_true_hor_t), intent(inout) :: this
    type(domain_t),            intent(in)    :: domain
    type(grid_field_t),        intent(inout) :: f
    !output:
    type(grid_field_t),        intent(inout) :: f1

    call this%grad%calc_grad(this%fx,this%fy,this%fz,f,domain)

    call this%co2contra%transform(this%fxt,this%fyt,this%fzt,this%fx,this%fy,this%fz,domain)

    !Project to the true horizontal plane:
    call this%fz%assign_prod(1.0_8,this%fz,this%hor_factor,domain%mesh_w)
    call this%fzt%update(-1.0_8,this%fz,domain%mesh_w)
    !apply boundary conditions for mass conservation
    call apply_zero_vert_boundary_conds(this%fzt,domain%mesh_w)

    call this%div%calc_div(f1,this%fxt, this%fyt, this%fzt, domain)

end subroutine calc_laplace

subroutine apply_zero_vert_boundary_conds(fzt,mesh)
    type(grid_field_t), intent(inout) :: fzt
    type(mesh_t),       intent(in)    :: mesh

    integer(kind=4) :: t, i, j, ks, ke

    do t = mesh%ts, mesh%te
        ks = mesh%tile(t)%ks; ke = mesh%tile(t)%ke
        do j = mesh%tile(t)%js, mesh%tile(t)%je
            do i = mesh%tile(t)%is, mesh%tile(t)%ie
                fzt%tile(t)%p(i,j,ks) = 0.0_8
                fzt%tile(t)%p(i,j,ke) = 0.0_8
            end do
        end do
    end do
end subroutine

end module
