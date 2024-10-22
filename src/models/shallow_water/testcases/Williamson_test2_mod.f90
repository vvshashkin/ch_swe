module Williamson_test2_mod

use test_fields_mod,     only : scalar_field_generator_t, &
                                set_vector_test_field, set_scalar_test_field, &
                                solid_rotation_t, ts2_height_generator_t,     &
                                zero_scalar_field_generator_t
use domain_mod,          only : domain_t
use generic_config_mod,  only : generic_config_t
use stvec_flexible_mod,  only : stvec_flexible_t, stvec_t
use grid_field_mod,      only : grid_field_t
use abstract_scorer_mod, only : scorer_t
use vec_math_mod,        only : l2norm

use const_mod,           only : pi, Earth_omega, Earth_radii, &
                                Earth_grav, Earth_sidereal_T
implicit none

real(kind=8), private, parameter :: &
     h_mean_default = 29400/Earth_grav,                 &
     u0_default = pi*Earth_radii/6.0_8/Earth_sidereal_T

type, extends(scorer_t) :: WT2_scorer_t
     type(grid_field_t) :: h0, u0, v0
     type(grid_field_t) :: err_h
     real(kind=8)       :: l2_h, linf_h
     contains
         procedure print_scores
 end type

contains

subroutine setup_Williamson_test2(state, scorer, config, v_components_type, domain)

    class(generic_config_t), intent(inout) :: config
    character(len=*),        intent(in)    :: v_components_type
    type(domain_t),          intent(in)    :: domain

    class(stvec_t),          intent(inout) :: state
    class(scorer_t),         intent(out), allocatable :: scorer

    real(kind=8) :: u0, h_mean, omega, a, grav, rotation_axis(3)
    logical      :: linear_mode
    type(grid_field_t), pointer :: h, u, v
    type(ts2_height_generator_t) :: height_field
    type(solid_rotation_t)       :: velocity_field
    type(WT2_scorer_t), allocatable :: wt2_scorer

    !get parameters
    call config%get(u0,"u0",default = u0_default)
    call config%get(h_mean,"h_mean",default = h_mean_default)
    call config%get(linear_mode,"linear_mode",default=.false.)

    omega = domain%metric%omega
    rotation_axis = domain%metric%rotation_axis
    a = domain%metric%scale

    !ensure that correct default values are set:
    if(domain%parcomm%myid == 0) then

        if(u0 /= u0_default) &
            print *, "WARNING: Williamson test 2 using non-standard u0 = ", u0

        if(h_mean /= h_mean_default) &
            print *, "WARNING: Williamson test 2 using non-standard h_mean = ", h_mean

        if(omega /= Earth_omega) &
            print *, "WARNING: Williamson test 2 using non-standard omega = ", omega

        if(a /= Earth_radii) &
            print *, "WARNING: Williamson test 2 using non-standard planet radius = ", a

    end if

    height_field = ts2_height_generator_t(h_mean = h_mean, &
                                          u0 = u0, omega = omega, a = a, &
                                          grav = Earth_grav, axis = rotation_axis, &
                                          linear_mode = linear_mode)

    velocity_field = solid_rotation_t(u0 = u0, axis = rotation_axis)

    select type(state)
    class is (stvec_flexible_t)

        call state%get_field(h,"h")
        call set_scalar_test_field(h, height_field, domain%mesh_p, 0)

        call state%get_field(u,"u")
        call state%get_field(v,"v")
        call set_vector_test_field(u, v, velocity_field, domain%mesh_u, &
                                   domain%mesh_v, 0, v_components_type)

        ! initialize scoring instance:
        allocate(wt2_scorer)

        wt2_scorer%h0    = h%create_similar()
        wt2_scorer%err_h = h%create_similar()
        wt2_scorer%u0    = u%create_similar()
        wt2_scorer%v0    = v%create_similar()

        call wt2_scorer%h0%assign(h,domain%mesh_p)
        wt2_scorer%l2_h   = l2norm(wt2_scorer%h0, domain%mesh_p, domain%parcomm)
        wt2_scorer%linf_h = wt2_scorer%h0%maxabs(domain%mesh_p, domain%parcomm)

        call wt2_scorer%u0%assign(u,domain%mesh_p)
        call wt2_scorer%v0%assign(v,domain%mesh_p)

        call move_alloc(wt2_scorer, scorer)

    end select

end subroutine

subroutine print_scores(this, state, domain, time)

    class(wt2_scorer_t), intent(inout) :: this
    class(stvec_t),      intent(in)    :: state
    type(domain_t),      intent(in)    :: domain
    real(kind=8),        intent(in)    :: time

    type(grid_field_t), pointer :: h, u, v
    real(kind=8)                :: l2_h, linf_h

    character(len=256) :: str_buff, fmt_str

    select type(state)
    class is (stvec_flexible_t)

        call state%get_field(h, "h")
        call this%err_h%assign(1.0_8, h, -1.0_8, this%h0, domain%mesh_p)

        l2_h   = l2norm(this%err_h, domain%mesh_p, domain%parcomm) / this%l2_h
        linf_h = this%err_h%maxabs(domain%mesh_p, domain%parcomm) / this%linf_h

    end select

    write(str_buff,*) "Errors, t=", time/3600.0_8, " hours, l2_h = ", l2_h, ", linf_h =", linf_h

    if(domain%parcomm%myid == 0) then
        write(fmt_str,"(A,I8.8,A)") "(A",len(trim(str_buff)),")"
        print trim(fmt_str), str_buff
    end if

end subroutine

end module Williamson_test2_mod