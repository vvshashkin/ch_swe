!Module defining abstract vector type for the iterative solvers

module vector_mod

use parcomm_mod, only : parcomm_global
use domain_mod,  only : domain_t

implicit none

private

type, public, abstract :: abstract_vector_t
private
contains
    procedure(create_similar_i), public, deferred :: create_similar
    procedure(copy_i),           public, deferred :: copy
    procedure(set_scalar_i),     public, deferred :: set_scalar   ! y = alpha
    procedure(update_i),         public, deferred :: update       ! y = y + alpha*x
    procedure(scale_i),          public, deferred :: scale        ! y = alpha*y
    procedure(assign_prod_i),    public, deferred :: assign_prod  ! z = x*y
    procedure(norm_i),           public, deferred :: norm
    procedure(dot_i),            public, deferred :: dot

    procedure(assign_s1v1_s2v2_i), public, deferred :: assign_s1v1_s2v2 ! z = alpha*x + beta*y
    procedure(assign_s1v1_i),      public, deferred :: assign_s1v1      ! z = alpha*x

    generic :: assign => assign_s1v1, assign_s1v1_s2v2

    procedure,                   public           :: output

end type abstract_vector_t

abstract interface
    subroutine create_similar_i(this, dest, domain)
        import abstract_vector_t, domain_t
        class(abstract_vector_t), intent(in)                 :: this
        class(abstract_vector_t), allocatable, intent(inout) :: dest
        type(domain_t), intent(in)                           :: domain
    end subroutine create_similar_i
end interface

abstract interface
    subroutine copy_i(this, source, domain)
        import abstract_vector_t, domain_t
        class(abstract_vector_t), intent(inout) :: this
        class(abstract_vector_t), intent(in)    :: source
        type(domain_t), intent(in)              :: domain
    end subroutine copy_i
end interface

abstract interface
    subroutine set_scalar_i(this, scalar, domain)
        import abstract_vector_t, domain_t
        class(abstract_vector_t), intent(inout) :: this
        real(kind=8),             intent(in)    :: scalar
        type(domain_t),           intent(in)    :: domain
    end subroutine set_scalar_i
end interface

abstract interface
    subroutine scale_i(this, scalar, domain)
        import abstract_vector_t, domain_t
        class(abstract_vector_t), intent(inout) :: this
        real(kind=8),             intent(in)    :: scalar
        type(domain_t),           intent(in)    :: domain
    end subroutine scale_i
end interface

abstract interface
    subroutine assign_prod_i(this, vec1, vec2, domain)
        import abstract_vector_t, domain_t
        class(abstract_vector_t), intent(inout) :: this
        class(abstract_vector_t), intent(in)    :: vec1, vec2
        type(domain_t),           intent(in)    :: domain
    end subroutine assign_prod_i
end interface

abstract interface
    subroutine update_i(this, scalar, vec, domain)
        import abstract_vector_t, domain_t
        class(abstract_vector_t), intent(inout) :: this
        real(kind=8),             intent(in)    :: scalar
        class(abstract_vector_t), intent(in)    :: vec
        type(domain_t),           intent(in)    :: domain
    end subroutine update_i
end interface

abstract interface
    subroutine assign_s1v1_s2v2_i(this, scalar1, vec1, scalar2, vec2, domain)
        import abstract_vector_t, domain_t
        class(abstract_vector_t), intent(inout) :: this
        real(kind=8),             intent(in)    :: scalar1, scalar2
        class(abstract_vector_t), intent(in)    :: vec1, vec2
        type(domain_t),           intent(in)    :: domain
    end subroutine assign_s1v1_s2v2_i
end interface

abstract interface
    subroutine assign_s1v1_i(this, scalar1, vec1, domain)
        import abstract_vector_t, domain_t
        class(abstract_vector_t), intent(inout) :: this
        real(kind=8),             intent(in)    :: scalar1
        class(abstract_vector_t), intent(in)    :: vec1
        type(domain_t),           intent(in)    :: domain
    end subroutine assign_s1v1_i
end interface

abstract interface
    function norm_i(this, domain) result(norm)
        import abstract_vector_t, domain_t
        class(abstract_vector_t), intent(in) :: this
        type(domain_t),           intent(in) :: domain
        real(kind=8)                         :: norm
    end function norm_i
end interface

abstract interface
    function dot_i(this, vec, domain) result(dot)
        import abstract_vector_t, domain_t
        class(abstract_vector_t), intent(in) :: this
        class(abstract_vector_t), intent(in) :: vec
        type(domain_t),           intent(in) :: domain
        real(kind=8)                         :: dot
    end function dot_i
end interface

contains

subroutine output(this, filename, domain)
    class(abstract_vector_t), intent(inout) :: this
    character(len=*),         intent(in)    :: filename
    type(domain_t),           intent(in)    :: domain

    call parcomm_global%abort("output routine is not implemented for the specific vector class")

end subroutine output

end module vector_mod
