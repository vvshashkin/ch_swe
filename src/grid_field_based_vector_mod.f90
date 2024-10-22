module grid_field_based_vector_mod

use parcomm_mod,             only : parcomm_global
use vector_mod,              only : abstract_vector_t
use grid_field_mod,          only : grid_field_t
use domain_mod,              only : domain_t
use mesh_mod,                only : mesh_t
use abstract_quadrature_mod, only : quadrature_t
use grid_field_factory_mod,  only : create_grid_field

implicit none


! Vector containing single grid-field, mesh and quadrature
type, public, extends(abstract_vector_t) :: grid_field_based_vector_t
    type(grid_field_t)               :: grid_field
    type(mesh_t),        pointer     :: mesh
    class(quadrature_t), allocatable :: quadrature
contains
    procedure :: init
    procedure :: create_similar
    procedure :: copy
    procedure :: set_scalar          ! y = alpha
    procedure :: update              ! y = y + alpha*x
    procedure :: assign_prod         ! z = x*y

    procedure :: assign_s1v1_s2v2    ! z = alpha*x + beta*y
    procedure :: assign_s1v1         ! z = alpha*x + beta*y

    procedure :: scale => scale_self ! y = alpha*y
    procedure :: norm
    procedure :: dot
    procedure :: output
end type grid_field_based_vector_t

contains

subroutine init(this, hw_xy, hw_z, mesh, quadrature)
    class(grid_field_based_vector_t), intent(inout) :: this
    integer(kind=4),                  intent(in)    :: hw_xy, hw_z
    type(mesh_t), pointer,            intent(in)    :: mesh
    type(quadrature_t),               intent(in)    :: quadrature

    this%mesh => mesh
    this%quadrature = quadrature
    call create_grid_field(this%grid_field, hw_xy, hw_z, mesh)

end subroutine init

subroutine create_similar(this, dest, domain)
    class(grid_field_based_vector_t),      intent(in)    :: this
    class(abstract_vector_t), allocatable, intent(inout) :: dest
    type(domain_t), intent(in)                           :: domain

    allocate(grid_field_based_vector_t :: dest)
    select type (dest)
    class is (grid_field_based_vector_t)
        dest%grid_field = this%grid_field%create_similar()
        dest%mesh => this%mesh
        dest%quadrature = this%quadrature
    class default
        call parcomm_global%abort("dest of wrong type")
    end select

end subroutine create_similar

subroutine copy(this, source, domain)
    class(grid_field_based_vector_t), intent(inout) :: this
    class(abstract_vector_t),         intent(in)    :: source
    type(domain_t),                   intent(in)    :: domain

    select type (source)
    class is (grid_field_based_vector_t)
        this%mesh => source%mesh
        this%grid_field = source%grid_field%copy()
    class default
        call parcomm_global%abort("source of wrong type")
    end select

end subroutine copy

subroutine set_scalar(this, scalar, domain)
    class(grid_field_based_vector_t), intent(inout) :: this
    real(kind=8),                     intent(in)    :: scalar
    type(domain_t),                   intent(in)    :: domain

    call this%grid_field%assign(scalar, this%mesh)

end subroutine set_scalar

subroutine scale_self(this, scalar, domain)
    class(grid_field_based_vector_t), intent(inout) :: this
    real(kind=8),                     intent(in)    :: scalar
    type(domain_t),                   intent(in)    :: domain

    call this%grid_field%assign(scalar, this%grid_field, this%mesh)

end subroutine scale_self

subroutine update(this, scalar, vec, domain)
    class(grid_field_based_vector_t), intent(inout) :: this
    real(kind=8),                     intent(in)    :: scalar
    class(abstract_vector_t),         intent(in)    :: vec
    type(domain_t),                   intent(in)    :: domain

    select type (vec)
    class is (grid_field_based_vector_t)
        call this%grid_field%update(scalar, vec%grid_field, this%mesh)
    class default
        call parcomm_global%abort("vec of wrong type")
    end select

end subroutine update

subroutine assign_s1v1_s2v2(this, scalar1, vec1, scalar2, vec2, domain)
    class(grid_field_based_vector_t), intent(inout) :: this
    real(kind=8),                     intent(in)    :: scalar1, scalar2
    class(abstract_vector_t),         intent(in)    :: vec1, vec2
    type(domain_t),                   intent(in)    :: domain

    select type (vec1)
    class is (grid_field_based_vector_t)
        select type (vec2)
        class is (grid_field_based_vector_t)
            call this%grid_field%assign(scalar1, vec1%grid_field, scalar2, vec2%grid_field, this%mesh)
        class default
            call parcomm_global%abort("vec2 of wrong type")
        end select
    class default
        call parcomm_global%abort("vec1 of wrong type")
    end select

end subroutine assign_s1v1_s2v2

subroutine assign_s1v1(this, scalar1, vec1, domain)
    class(grid_field_based_vector_t), intent(inout) :: this
    real(kind=8),                     intent(in)    :: scalar1
    class(abstract_vector_t),         intent(in)    :: vec1
    type(domain_t),                   intent(in)    :: domain

    select type (vec1)
    class is (grid_field_based_vector_t)
        call this%grid_field%assign(scalar1, vec1%grid_field, this%mesh)
    class default
        call parcomm_global%abort("vec1 of wrong type")
    end select

end subroutine assign_s1v1

subroutine assign_prod(this, vec1, vec2, domain)
    class(grid_field_based_vector_t), intent(inout) :: this
    class(abstract_vector_t),         intent(in)    :: vec1, vec2
    type(domain_t),                   intent(in)    :: domain

    select type (vec1)
    class is (grid_field_based_vector_t)
        select type (vec2)
        class is (grid_field_based_vector_t)
            call this%grid_field%assign_prod(1.0_8, vec1%grid_field, vec2%grid_field, this%mesh)
        class default
            call parcomm_global%abort("vec2 of wrong type")
        end select
    class default
        call parcomm_global%abort("vec1 of wrong type")
    end select

end subroutine assign_prod

function norm(this, domain)
    class(grid_field_based_vector_t), intent(in) :: this
    type(domain_t),                   intent(in) :: domain
    real(kind=8)                                 :: norm

    norm = this%quadrature%l2norm(this%grid_field, this%mesh, domain%parcomm)

end function norm

function dot(this, vec, domain)
    class(grid_field_based_vector_t), intent(in) :: this
    class(abstract_vector_t),         intent(in) :: vec
    type(domain_t),                   intent(in) :: domain
    real(kind=8)                                 :: dot

    select type (vec)
    class is (grid_field_based_vector_t)
        dot = this%quadrature%dot(this%grid_field, vec%grid_field, this%mesh, domain%parcomm)
    class default
        call parcomm_global%abort("vec of wrong type")
    end select

end function dot
subroutine output(this, filename, domain)

    use outputer_factory_mod,  only : create_master_paneled_outputer
    use outputer_abstract_mod, only : outputer_t

    class(grid_field_based_vector_t), intent(inout) :: this
    character(len=*),                 intent(in)    :: filename
    type(domain_t),                   intent(in)    :: domain

    class(outputer_t), allocatable :: outputer
    ! WORKAROUND for debuging purposes
    call create_master_paneled_outputer(outputer, "p", domain)

    call outputer%write(this%grid_field, domain, filename, 1)

end subroutine output
end module grid_field_based_vector_mod
