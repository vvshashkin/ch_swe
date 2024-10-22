module abstract_particle_interp_mod

use domain_mod,          only : domain_t
use particles_mod,       only : particles_t
use particle_values_mod, only : particle_values_t
use grid_field_mod,      only : tile_field_t
use string_mod,          only : find_str_index

implicit none

type, abstract :: particle_interp_t

    contains

    procedure(calc_w_i), deferred :: calc_w

    procedure :: interp_fld_name
    procedure :: interp_fld_ind
    procedure :: interp_2d_levs

    generic :: interp => interp_fld_name, interp_fld_ind, interp_2d_levs

end type

type particle_interp_container_t
    class(particle_interp_t), allocatable :: interp
end type

abstract interface

    subroutine calc_w_i(this, particles, t, domain)

        import particle_interp_t, particles_t, domain_t

        class(particle_interp_t), intent(inout) :: this
        type(particles_t),        intent(in)    :: particles
        integer(kind=4),          intent(in)    :: t
        type(domain_t),           intent(in)    :: domain

    end subroutine

end interface

contains

subroutine interp_fld_name(this, f_interp, fin, field_name, domain)

    class(particle_interp_t), intent(in)    :: this
    type(particle_values_t),  intent(inout) :: f_interp
    type(tile_field_t),       intent(in)    :: fin
    character(len=*),         intent(in)    :: field_name
    type(domain_t),           intent(in)    :: domain

    integer(kind=4) :: fld_ind
    
    fld_ind = find_str_index(field_name, f_interp%field_names)

    if(fld_ind == 0) call domain%parcomm%abort("interp_fld_name, cannot find field: "//field_name)

    call this%interp_fld_ind(f_interp, fin, fld_ind, domain)

end subroutine

subroutine interp_fld_ind(this, f_interp, fin, field_ind, domain)

    class(particle_interp_t), intent(in)    :: this
    type(particle_values_t),  intent(inout) :: f_interp
    type(tile_field_t),       intent(in)    :: fin
    integer(kind=4),          intent(in)    :: field_ind
    type(domain_t),           intent(in)    :: domain


    call domain%parcomm%abort("interp_fld_ind not implemented")

end subroutine

subroutine interp_2d_levs(this, f_interp, fin, domain)

    class(particle_interp_t), intent(in)    :: this
    type(particle_values_t),  intent(inout) :: f_interp
    type(tile_field_t),       intent(in)    :: fin
    type(domain_t),           intent(in)    :: domain

    call domain%parcomm%abort("interp_fld_2d_levs not implemented")

end subroutine

end module