module barotropic_inst_mod

use domain_mod,                 only : domain_t
use domain_factory_mod,         only : create_domain
use stvec_mod,                  only : stvec_t
use stvec_swm_mod,              only : stvec_swm_t
use stvec_swm_factory_mod,      only : create_stvec_swm
use operator_mod,               only : operator_t
use operator_swm_factory_mod,   only : create_swm_operator
use timescheme_mod,             only : timescheme_t
use timescheme_factory_mod,     only : create_timescheme
use outputer_abstract_mod,      only : outputer_t, outputer_vector_t
use outputer_factory_mod,       only : create_master_paneled_outputer,&
                                     create_latlon_outputer, create_latlon_vec_outputer
use parcomm_mod,                only : parcomm_global
use generic_config_mod,         only : generic_config_t
use config_tools_mod,           only : get_config_type_by_extension

use operator_swm_mod,           only : operator_swm_t
use operator_swm_vec_inv_imex_mod, only : operator_swm_vec_inv_imex_t

use operator_swm_diff_mod,         only : operator_swm_diff_t
use operator_swm_diff_factory_mod, only : create_swm_diff_operator

use const_mod,  only : pi, Earth_grav, Earth_omega, Earth_radii, pi, Earth_sidereal_T

use test_fields_mod, only : barotropic_instability_height_generator_t, &
                            barotropic_instability_wind_generator_t
use key_value_mod, only : key_value_r8_t

use grid_field_mod,         only : grid_field_t
use grid_field_factory_mod, only : create_grid_field
use abstract_curl_mod,      only : curl_operator_t
use curl_factory_mod,       only : create_curl_operator
use abstract_co2contra_mod, only : co2contra_operator_t
use co2contra_factory_mod,  only : create_co2contra_operator

use mpi

implicit none

real(kind=8), parameter :: grav   = Earth_grav
real(kind=8), parameter :: H0     = 11000._8
real(kind=8), parameter :: u0     = 80.0_8
real(kind=8), parameter :: h_pert = 120.0_8

type(barotropic_instability_height_generator_t) :: height_field
type(barotropic_instability_wind_generator_t)   :: velocity_field

contains

subroutine run_barotropic_inst()

    use vec_math_mod, only : l2norm
    use namelist_read_mod, only : read_namelist_as_str

    class(generic_config_t), allocatable :: config, subconfig
    type(domain_t)     :: domain

    class(stvec_t),           allocatable :: state, state_err, state_ex
    class(operator_t),        allocatable :: operator
    class(timescheme_t),      allocatable :: timescheme
    class(outputer_t),        allocatable :: outputer, outputer_curl
    class(outputer_vector_t), allocatable :: outputer_vec

    type(grid_field_t) :: curl

    type(operator_swm_diff_t),   allocatable :: operator_diff
    class(timescheme_t),         allocatable :: timescheme_diff

    integer(kind=4) :: halo_width = 8

    character(:), allocatable :: config_file, config_string, timescheme_name, name, &
                                 v_components_type
    type(key_value_r8_t) :: diagnostics

    class(curl_operator_t), allocatable :: curl_op
    class(co2contra_operator_t), allocatable :: co2contra_op
    type(grid_field_t) :: u, v

    real(kind=8)      :: a, omega
    real(kind=8), allocatable :: rotation_axis(:)
    real(kind=8)      :: dt, tau_write, tau_diagnostics, simulation_time
    integer(kind=4)   :: nstep_write, nstep_diagnostics

    real(kind=8)    :: time, l2err, l2_ex, wall_time
    real(kind=8)    :: l2_err_h, l2_ex_h, l2_err_u, l2_ex_u
    real(kind=8)    :: linf_err_h, linf_ex_h, linf_err_u, linf_ex_u

    integer(kind=4) :: it, N, N_quad_points
    logical :: found

    !get and parse config from file
    config_file = "swm_model.cfg"
    call read_namelist_as_str(config_string, config_file, parcomm_global%myid)
    config = get_config_type_by_extension(config_file)
    call config%parse(config_string)

    call config%get(N,"domain%N")
    N_quad_points = 10*N

    !Ensure that correct default values are set:
    call config%get(omega, "domain%metric%omega",found=found,default=Earth_omega)
    if(.not. found) call config%set("domain%metric%omega",Earth_omega)
    if(found .and. omega /= Earth_omega) &
        print *, "WARNING: BAROTROPIC_INST using non-standard omega = ", omega

    call config%get(a, "domain%metric%planet_radius",found=found,default=Earth_radii)
    if(.not. found) call config%set("domain%metric%planet_radius",Earth_radii)
    if(found .and. a /= Earth_radii) &
        print *, "WARNING: BAROTROPIC_INST using non-standard planet radius = ", a

    call config%get(rotation_axis,"domain%metric%rotation_axis",found=found, &
                    default = [0.0_8, 0.0_8, 1.0_8])
    if(.not. found) call config%set("domain%metric%rotation_axis",rotation_axis)
    if(found .and. (rotation_axis(1) /= 0.0_8 .or. &
                    rotation_axis(2) /= 0.0_8 .or. rotation_axis(3) /= 1.0_8)) &
        print *, "WARNING: BAROTROPIC_INST using non-standard rotation axis = ", rotation_axis

    !get time parameters
    call config%get(dt,"dt")
    call config%get(tau_write,"tau_write")
    call config%get(tau_diagnostics,"tau_diagnostics")
    nstep_write = nint(tau_write/dt)
    nstep_diagnostics = nint(tau_diagnostics/dt)
    call config%get(simulation_time,"simulation_time")

    !Set field generators
    height_field = create_barotropic_instability_height_field_generator(H0 = H0, &
                              omega = omega, grav = Earth_grav, &
                              h_pert = h_pert, u0 = u0,   &
                              a = a, Nq = N_quad_points)

    velocity_field = barotropic_instability_wind_generator_t(u0=u0)

    if (parcomm_global%myid==0) print*, "Gravity Wave CFL = ", &
            real(dt*sqrt(H0* grav)    &
            /(2*pi/4/N)/a,4)

    call config%get(subconfig,"domain")
    call create_domain(domain, subconfig)

    call create_stvec_swm(state,     domain, halo_width, 0)
    call create_stvec_swm(state_ex,  domain, 0         , 0)
    call create_stvec_swm(state_err, domain, 0         , 0)

    call config%get(name,"swm_operator_section")
    call config%get(subconfig,name)
    call subconfig%get(v_components_type,"v_components_type")
    call create_swm_operator(operator, Earth_grav, subconfig, domain)

    call config%get(timescheme_name,"timescheme_name")
    call create_timescheme(timescheme, state, timescheme_name)

    call config%get(name,"diffusion_section")
    call config%get(subconfig,name)
    call config%get(timescheme_name,name//"%diff_time_scheme")
    call create_swm_diff_operator(operator_diff, subconfig, domain)
    call create_timescheme(timescheme_diff, state, timescheme_name)

    if(trim(domain%horizontal_staggering)=="C") then
        call create_curl_operator(curl_op,"curl_c_sbp_sat_42",domain)
        co2contra_op = create_co2contra_operator(domain,"co2contra_c_sbp42_new")
        call create_grid_field(u,7,0,domain%mesh_u)
        call create_grid_field(v,7,0,domain%mesh_v)
    end if

    if(domain%horizontal_staggering == "Ah") then
        call create_latlon_outputer(outputer,          2*domain%partition%Nh+1, 4*domain%partition%Nh, "xy", domain)
        call create_latlon_outputer(outputer_curl,     2*domain%partition%Nh+1, 4*domain%partition%Nh, "xy", domain)
        call create_latlon_vec_outputer(outputer_vec,  2*domain%partition%Nh+1, 4*domain%partition%Nh, "xy", "xy", &
                                   v_components_type, domain)
    else if(domain%horizontal_staggering == "Ch") then
        call create_latlon_outputer(outputer, 2*domain%partition%Nh+1, 4*domain%partition%Nh, "xy", domain)
        call create_latlon_vec_outputer(outputer_vec,  2*domain%partition%Nh+1, 4*domain%partition%Nh, "y", "x", &
                                   v_components_type, domain)
    else if(domain%horizontal_staggering == "C") then
        call create_latlon_outputer(outputer,      2*domain%partition%Nh+1, 4*domain%partition%Nh, "o", domain)
        call create_latlon_outputer(outputer_curl, 2*domain%partition%Nh+1, 4*domain%partition%Nh, "xy", domain)
        call create_latlon_vec_outputer(outputer_vec,  2*domain%partition%Nh+1, 4*domain%partition%Nh, "x", "y", &
                                       v_components_type, domain)
    else
        call parcomm_global%abort("This staggering is not implemented in"//&
                                  " barotropic instability swe test output:"//&
                                  domain%horizontal_staggering)
    end if

    call get_exact_solution(state, domain, v_components_type)
    call state_ex%assign(1.0_8,state,0.0_8,state,domain)

    call create_grid_field(curl, halo_width, 0, domain%mesh_q)

    select type(state)
    class is (stvec_swm_t)
        select type(operator)
        class is (operator_swm_t)
            call operator%curl_op%calc_curl(curl, state%u, state%v, domain)
            call outputer_curl%write(curl, domain, 'curl.dat', 1)
        class is (operator_swm_vec_inv_imex_t)
            call operator%curl_op%calc_curl(curl, state%u, state%v, domain)
            call outputer_curl%write(curl, domain, 'curl.dat', 1)
        class default
            if(trim(domain%horizontal_staggering)=="C") then
                call co2contra_op%transform2co(u,v,state%u,state%v,domain)
                call curl_op%calc_curl(curl,u,v,domain)
                call outputer_curl%write(curl, domain, 'curl.dat', 1)
            end if
        end select
        call outputer%write(state%h, domain, 'h.dat',1)
        call outputer_vec%write(state%u,state%v,domain,"u.dat","v.dat",1)
    end select

    select type(state_ex)
    class is (stvec_swm_t)
        l2_ex_h = l2norm(state_ex%h,  domain%mesh_p, domain%parcomm)
        l2_ex_u = sqrt(l2norm(state_ex%u,  domain%mesh_u, domain%parcomm)**2+&
                       l2norm(state_ex%v,  domain%mesh_v, domain%parcomm)**2)
        linf_ex_h = state_ex%h%maxabs(domain%mesh_p, domain%parcomm)
        linf_ex_u = max(state_ex%u%maxabs(domain%mesh_u, domain%parcomm), &
                        state_ex%v%maxabs(domain%mesh_v, domain%parcomm))
    end select

    wall_time = mpi_wtime()
    do it = 1, int(simulation_time/dt)

        call timescheme%step(state, operator, domain, dt)
        !call timescheme_diff%step(state, operator_diff, domain, dt)

        time = it*dt

        call state_err%assign(1.0_8,state,-1.0_8,state_ex,domain)

        if(mod(it, nstep_diagnostics) == 0) then
            diagnostics = operator%get_diagnostics(state, domain)
            if(parcomm_global%myid == 0) call diagnostics%print()

            select type(state_err)
            class is (stvec_swm_t)
                l2_err_h = l2norm(state_err%h,  domain%mesh_p, domain%parcomm) / l2_ex_h
                l2_err_u = sqrt(l2norm(state_err%u,  domain%mesh_u, domain%parcomm)**2+&
                                l2norm(state_err%v,  domain%mesh_v, domain%parcomm)**2) / l2_ex_u
                linf_err_h = state_err%h%maxabs(domain%mesh_p, domain%parcomm) / linf_ex_h
                linf_err_u = max(state_err%u%maxabs(domain%mesh_u, domain%parcomm), &
                                 state_err%v%maxabs(domain%mesh_v, domain%parcomm)) / linf_ex_u
                if (parcomm_global%myid==0) print '(A,F12.4,4(A,E15.7))', &
                                                  "Hours = ", real(time/3600 ,4), &
                                                  " l2_h =", real(l2_err_h,4), " linf_h = ", real(linf_err_h,4), &
                                                  " l2_u =", real(l2_err_u,4), " linf_u = ", real(linf_err_u,4)
                if(domain%parcomm%myid == 0) print *, "step wall-time: ", mpi_wtime()-wall_time
                wall_time = mpi_wtime()
            end select
        end if

        if(mod(it,nstep_write) == 0) then

            select type(state)
            class is (stvec_swm_t)

                select type(operator)
                class is (operator_swm_t)
                    call operator%curl_op%calc_curl(curl, state%u, state%v, domain)
                    call outputer_curl%write(curl, domain, 'curl.dat', int(it/nstep_write)+1)
                class is (operator_swm_vec_inv_imex_t)
                    call operator%curl_op%calc_curl(curl, state%u, state%v, domain)
                    call outputer_curl%write(curl, domain, 'curl.dat', int(it/nstep_write)+1)
                class default
                    if(trim(domain%horizontal_staggering)=="C") then
                        call co2contra_op%transform2co(u,v,state%u,state%v,domain)
                        call curl_op%calc_curl(curl,u,v,domain)
                        call outputer_curl%write(curl, domain, 'curl.dat', int(it/nstep_write)+1)
                    end if
                end select
                call outputer%write(state%h, domain, 'h.dat', int(it/nstep_write)+1)
                call outputer_vec%write(state%u,state%v,domain,"u.dat","v.dat",int(it/nstep_write)+1)
                l2err = l2norm(state%h, domain%mesh_p, domain%parcomm)
                if (parcomm_global%myid==0) print*, "Hours = ", real(time/3600 ,4), &
                                                    "rec ", int(it/nstep_write)+1
            end select
        end if
    end do

end subroutine run_barotropic_inst

subroutine get_exact_solution(state, domain, v_components_type)

    use test_fields_mod, only : set_vector_test_field, &
                                set_scalar_test_field

    class(stvec_t), intent(inout) :: state
    type(domain_t), intent(in)    :: domain
    character(*),   intent(in)    :: v_components_type

    select type(state)
    class is (stvec_swm_t)

        call set_scalar_test_field(state%h, height_field, domain%mesh_p, 0)
        call set_vector_test_field(state%u, state%v, velocity_field, &
              domain%mesh_u, domain%mesh_v, 0, v_components_type)

    end select


end subroutine get_exact_solution

function create_barotropic_instability_height_field_generator(H0, &
                              omega, grav, h_pert, u0, a, Nq) result (height_gen)

    use const_mod, only : pi
    use barotropic_instability_u_mod, only : barotropic_instability_u

    real(kind=8), intent(in) :: H0
    real(kind=8), intent(in) :: omega
    real(kind=8), intent(in) :: grav
    real(kind=8), intent(in) :: h_pert
    real(kind=8), intent(in) :: u0
    real(kind=8), intent(in) :: a
    integer(kind=4), intent(in) :: Nq

    type(barotropic_instability_height_generator_t) :: height_gen

    real(kind=8), parameter :: phi0 = pi/7d0, phi1 = .5d0*pi-phi0
    real(kind=8)            :: dphi
    real(kind=8)            :: phi, u, dh
    real(kind=8)            :: urhs
    integer(kind=4)         :: j, iq

    real(kind=8), parameter :: xq(3) = [-sqrt(0.6_8), 0.0_8, sqrt(0.6_8)]
    real(kind=8), parameter :: wq(3) = [5.0_8/9.0_8, 8.0_8/9.0_8, 5.0_8/9.0_8]

    height_gen%H0 = H0
    height_gen%Nq = Nq
    allocate(height_gen%H_zonal(-1:Nq+1))
    dphi = (phi1 - phi0) / real(Nq,8)
    height_gen%dphi = dphi
    height_gen%h_pert = h_pert

    height_gen%H_zonal(-1) = H0!0.0_8
    height_gen%H_zonal( 0) = H0!0.0_8
    do j=1, Nq
        dh = 0.0_8
        do iq=1, size(wq)
            phi   = phi0+(j-0.5_8)*dphi+dphi*xq(iq)*0.5_8
            u = barotropic_instability_u(u0, phi)
            urhs  =-(u**2*tan(phi)+2._8*omega*a*sin(phi)*u) / grav
            dh = dh+wq(iq)*urhs
        end do
        height_gen%H_zonal(j) = height_gen%H_zonal(j-1)+0.5_8*dh*dphi
    end do
    height_gen%H_zonal(Nq+1) = height_gen%H_zonal(Nq)

    height_gen%H_north = height_gen%H_zonal(Nq)

end function

end module barotropic_inst_mod
