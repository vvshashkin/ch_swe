module shallow_water_mod

use domain_mod,                    only : domain_t
use domain_factory_mod,            only : create_domain
use stvec_mod,                     only : stvec_t
use stvec_swm_mod,                 only : stvec_swm_t
use stvec_swm_factory_mod,         only : create_stvec_swm
use operator_mod,                  only : operator_t
use operator_swm_factory_mod,      only : create_swm_operator
use swm_output_diag_mod,           only : swm_output_diag_t
use swm_output_diag_factory_mod,   only : create_swm_output_diag
use timescheme_mod,                only : timescheme_t
use timescheme_factory_mod,        only : create_timescheme
use generic_config_mod,            only : generic_config_t
use test_fields_mod,               only : scalar_field_generator_t
use operator_swm_diff_mod,         only : operator_swm_diff_t
use operator_swm_diff_factory_mod, only : create_swm_diff_operator
use swm_forcing_mod,               only : swm_forcing_t
use vec_math_mod,                  only : l2norm
use const_mod,                     only : Earth_grav
use abstract_scorer_mod,           only : scorer_t

use timer_mod, only : init_timer_mod, stop_timer, start_timer, print_timer

implicit none

contains

subroutine run_shallow_water_model(config)

    class(generic_config_t), intent(inout) :: config

    class(generic_config_t), allocatable :: subconfig
    type(domain_t) :: domain

    class(stvec_t),           allocatable :: state, state_err, state_ex
    class(operator_t),        allocatable :: operator
    class(swm_forcing_t),     allocatable :: forcing
    class(timescheme_t),      allocatable :: timescheme, timescheme_diff
    class(scorer_t),          allocatable :: scorer

    character(len=:), allocatable :: scores_str

    type(operator_swm_diff_t),   allocatable :: operator_diff

    type(swm_output_diag_t) :: output_diag

    class(scalar_field_generator_t), allocatable :: orography_generator

    integer(kind=4) :: halo_width = 8

    character(:), allocatable :: timescheme_name, name, v_components_type

    real(kind=8)     :: dt, tau_write, simulation_time, tau_diagnostics
    integer(kind=4)  :: nstep_write, nstep_diagnostics, N

    real(kind=8)    :: time
    integer(kind=4) :: it

    !get time parameters
    call config%get(dt,"dt")
    call config%get(tau_write,"tau_write")
    call config%get(tau_diagnostics,"tau_diagnostics")
    nstep_write = nint(tau_write/dt)
    nstep_diagnostics = nint(tau_diagnostics/dt)
    call config%get(simulation_time,"simulation_time")

    !create domain
    call config%get(subconfig,"domain")
    call create_domain(domain, subconfig)

    !determine components type of prognostic wind
    call config%get(name,"swm_operator_section")
    call config%get(v_components_type,name//"%v_components_type")

    !create stvec
    call create_stvec_swm(state,     domain, halo_width, 0)
    call create_stvec_swm(state_ex,  domain, halo_width, 0)
    call create_stvec_swm(state_err, domain, halo_width, 0)

    !setup test case (initial conditions & orography)
    call config%get(name,"testcase_section")
    call config%get(subconfig,name)
    call setup_testcase(state, orography_generator, forcing, scorer, &
                        subconfig, v_components_type, domain)

    call config%get(name,"swm_operator_section")
    call config%get(subconfig,name)
    call create_swm_operator(operator, Earth_grav, subconfig, domain, &
                             orography_generator)    

    call config%get(timescheme_name,"timescheme_name")
    call create_timescheme(timescheme, state, timescheme_name)

    call config%get(name,"diffusion_section",default="None")
    if(name /= "None") then
        call config%get(subconfig,name)
        call config%get(timescheme_name,name//"%diff_time_scheme")
        call create_swm_diff_operator(operator_diff, subconfig, domain)
        call create_timescheme(timescheme_diff, state, timescheme_name)
    end if

    call config%get(name,"output_diag_section")
    call config%get(subconfig,name,default=config%get_empty_subconfig())
    call create_swm_output_diag(output_diag, subconfig, orography_generator, &
                                v_components_type, domain)

    call output_diag%print_cfl_diag(state, dt, domain)
    call output_diag%write_fields(state, 1, domain)
    if(allocated(scorer)) call scorer%print_scores(state, domain, 0.0_8)

    do it = 1, int(simulation_time/dt)

        if(allocated(forcing)) call forcing%apply(state, dt, domain)

        call timescheme%step(state, operator, domain, dt)

        if(allocated(operator_diff)) &
            call timescheme_diff%step(state, operator_diff, domain, dt)


        if(mod(it, nstep_diagnostics) == 0) then

            time = it*dt

            call output_diag%print_integrals_diag(state, domain)
            call output_diag%print_integral_tends(operator, state, domain, "dynamics")
            if(allocated(operator_diff)) &
                call output_diag%print_integral_tends(operator_diff, state, domain, "diffusion")

            if(allocated(scorer)) call scorer%print_scores(state, domain, time)

        end if

        if(mod(it,nstep_write) == 0) then

            call output_diag%write_fields(state, int(it/nstep_write)+1, domain)
            if(domain%parcomm%myid == 0) &
                print *, "rec", int(it/nstep_write)+1

        end if

    end do

end subroutine run_shallow_water_model

subroutine setup_testcase(state, orography_generator, forcing, scorer, &
                          config, v_components_type, domain)

    use Williamson_test2_mod,          only : setup_Williamson_test2
    use Williamson_test5_mod,          only : setup_Williamson_test5
    use RH4_testcase_mod,              only : setup_RH4_test
    use barotropic_instability_mod,    only : setup_barotropic_instability_test
    use Eldred_testcase_mod,           only : setup_Eldred_testcase
    use gaussian_hill_linear_test_mod, only : setup_Gaussian_hill_linear_test
    use test_fields_mod,               only : zero_scalar_field_generator_t
    use oscillating_gaussian_hill_test_mod, only : setup_oscillating_Gaussian_hill_test

    class(generic_config_t), intent(inout) :: config
    character(len=*),        intent(in)    :: v_components_type
    type(domain_t),          intent(in)    :: domain

    class(stvec_t),                               intent(inout) :: state
    class(scalar_field_generator_t), allocatable, intent(out)   :: orography_generator
    class(swm_forcing_t),            allocatable, intent(out)   :: forcing
    class(scorer_t),                 allocatable, intent(out)   :: scorer

    character(len=:), allocatable :: testcase_name

    call config%get(testcase_name,"testcase_name")

    select case(testcase_name)
    case("Williamson_test2", "TS2", "test2")
        call setup_Williamson_test2(state, scorer, config, v_components_type, domain)
    case("MIRW","Williamson_test5")
        call setup_Williamson_test5(state, orography_generator, config, v_components_type, domain)
    case("RH4","Williamson_test6")
        call setup_RH4_test(state, config,v_components_type,domain)
    case("Galewsky", "barotropic_instability")
        call setup_barotropic_instability_test(state, config, v_components_type, domain)
    case("Eldred_test", "Eldred_testcase")
        call setup_Eldred_testcase(state,orography_generator,forcing,config,v_components_type,domain)
    case("Gaussian_hill_linear")
        call setup_Gaussian_hill_linear_test(state, scorer, config, v_components_type, domain)
    case("Oscillating_gauss_linear")
        call setup_oscillating_gaussian_hill_test(state, config, v_components_type, domain)
    case default
        call domain%parcomm%abort("shallow water model setup_testcase, unknown case name: "//testcase_name)
    end select

    if(.not. allocated(orography_generator)) orography_generator = zero_scalar_field_generator_t()

end subroutine

end module shallow_water_mod