module halo_default_exchange_mod

!halo-procedure type interface to exchanges

use halo_mod,          only : halo_t, halo_vec_t
use exchange_halo_mod, only : exchange_t

implicit none

type, extends(halo_t) :: halo_default_exchange_t

    class(exchange_t), allocatable  :: exch_halo

    contains

    procedure :: get_halo_scalar => get_default_exchange_scalar_halo

end type

type, extends(halo_vec_t) :: halo_vec_default_exchange_t

    class(exchange_t), allocatable  :: exch_halo

    contains

    procedure :: get_halo_vector => get_default_exchange_vector_halo

end type

contains

subroutine get_default_exchange_scalar_halo(this,f,domain,halo_width)
    use grid_field_mod, only : grid_field_t
    use domain_mod,     only : domain_t

    class(halo_default_exchange_t),  intent(inout) :: this
    class(grid_field_t),             intent(inout) :: f
    type(domain_t),                  intent(in)    :: domain
    integer(kind=4),                 intent(in)    :: halo_width

    call this%exch_halo%do(f, domain%parcomm)

end subroutine get_default_exchange_scalar_halo

subroutine get_default_exchange_vector_halo(this,u,v,domain,halo_width)
    use grid_field_mod, only : grid_field_t
    use domain_mod,     only : domain_t

    class(halo_vec_default_exchange_t),  intent(inout) :: this
    class(grid_field_t),                 intent(inout) :: u, v
    type(domain_t),                      intent(in)    :: domain
    integer(kind=4),                     intent(in)    :: halo_width

    call this%exch_halo%do_vec(u, v, domain%parcomm)

end subroutine get_default_exchange_vector_halo

end module
