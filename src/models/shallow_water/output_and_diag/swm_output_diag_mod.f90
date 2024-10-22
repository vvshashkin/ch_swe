module swm_output_diag_mod

use domain_mod,                   only : domain_t
use stvec_flexible_mod,           only : stvec_flexible_t, stvec_t
use outputer_abstract_mod,        only : outputer_t, outputer_vector_t

use abstract_grad_mod,            only : grad_operator_t
use abstract_div_mod,             only : div_operator_t
use abstract_curl_mod,            only : curl_operator_t
use abstract_co2contra_mod,       only : co2contra_operator_t
use abstract_interpolators2d_mod, only : interpolator2d_scalar2vec_t, &
                                         interpolator2d_vec2vec_t
use abstract_quadrature_mod,      only : quadrature_t
use grid_field_mod,               only : grid_field_t
use operator_mod,                 only : operator_t

use const_mod,                    only : Earth_grav

implicit none

type, public :: swm_output_diag_t

    class(outputer_t),        allocatable :: outputer_h, outputer_curl
    class(outputer_vector_t), allocatable :: outputer_v
    logical :: output_div, output_curl
    character(len=:), allocatable :: v_components_type

    class(curl_operator_t),             allocatable :: curl_op
    class(div_operator_t),              allocatable :: div_op
    class(co2contra_operator_t),        allocatable :: co2contra_op
    class(interpolator2d_scalar2vec_t), allocatable :: interp_h2uv
    class(interpolator2d_vec2vec_t),    allocatable :: interp_uv2q

    class(quadrature_t),         allocatable :: quadrature_h, quadrature_u, &
                                                quadrature_v, quadrature_q

    type(grid_field_t) :: h, u, v, curl, div, hu, hv, h_surf, hqu, hqv, fcori
    type(grid_field_t) :: hq, hq_tend, curl_tend

    contains

        procedure :: write_fields
        procedure :: print_cfl_diag
        procedure :: print_integrals_diag
        procedure :: print_integral_tends
        procedure :: calc_energies

end type swm_output_diag_t

contains

subroutine write_fields(this, state, irec, domain)

    class(swm_output_diag_t), intent(inout) :: this
    class(stvec_t),           intent(inout) :: state
    integer(kind=4),          intent(in)    :: irec
    type(domain_t),           intent(in)    :: domain

    type(grid_field_t), pointer :: h, u, v

    select type(state)
    class is (stvec_flexible_t)
        call state%get_field(h,"h")
        call state%get_field(u,"u")
        call state%get_field(v,"v")
    class default
        call domain%parcomm%abort("shallow water outputter works only with stvec_flexible_t and its extensions")
    end select

    if(allocated(this%outputer_h)) then
        call this%h%assign(1.0_8, h, 1.0_8, this%h_surf, domain%mesh_p)
        call this%outputer_h%write(this%h, domain, 'h.dat', irec)
    end if

    if(allocated(this%outputer_v)) &
        call this%outputer_v%write(u, v, domain, 'u.dat', 'v.dat', irec)

    if(this%output_curl) then

        if(this%v_components_type == "contravariant") then
            call this%co2contra_op%transform2co(this%u, this%v, u, v, domain)
        else
            call this%u%assign(u, domain%mesh_u)
            call this%v%assign(v, domain%mesh_v)
        end if

        call this%curl_op%calc_curl(this%curl, this%u, this%v, domain)

        call this%outputer_curl%write(this%curl, domain, "curl.dat", irec)

    end if

    if(this%output_div) then

        if(this%v_components_type == "covariant") then
            call this%co2contra_op%transform(this%u, this%v, u, v, domain)
        else
            call this%u%assign(u, domain%mesh_u)
            call this%v%assign(v, domain%mesh_v)
        end if

        call this%div_op%calc_div(this%div, this%u, this%v, domain)

        call this%outputer_h%write(this%div, domain, "div.dat", irec)

    end if
    
end subroutine write_fields

subroutine print_CFL_diag(this, state, dt, domain)

    class(swm_output_diag_t), intent(inout) :: this
    class(stvec_t),           intent(inout) :: state
    real(kind=8),             intent(in)    :: dt
    type(domain_t),           intent(in)    :: domain

    type(grid_field_t), pointer :: h, u, v
    real(kind=8) :: cfl_gv, cfl_vx, cfl_vy, dx, dy
    integer(kind=4) :: t

    select type(state)
    class is (stvec_flexible_t)
        call state%get_field(h,"h")
        call state%get_field(u,"u")
        call state%get_field(v,"v")
    class default
        call domain%parcomm%abort("shallow water diagnostics works only with stvec_flexible_t and its extensions")
    end select

    if(this%v_components_type == "covariant") then
        call this%co2contra_op%transform(this%u, this%v, u, v, domain)
    else
        call this%co2contra_op%transform2co(this%u, this%v, u, v, domain)
    end if

    do t = domain%mesh_p%ts, domain%mesh_p%te

        dx = domain%mesh_p%scale*domain%mesh_p%tile(t)%hx
        dy = domain%mesh_p%scale*domain%mesh_p%tile(t)%hy

        call this%h%tile(t)%assign(Earth_grav/min(dx,dy)**2, h%tile(t), domain%mesh_p%tile(t))
        call this%u%tile(t)%assign_prod(1._8/dx**2,this%u%tile(t),u%tile(t),domain%mesh_u%tile(t))
        call this%v%tile(t)%assign_prod(1._8/dy**2,this%v%tile(t),v%tile(t),domain%mesh_v%tile(t))
    
    end do
    
    cfl_gv = sqrt(this%h%maximum(domain%mesh_p, domain%parcomm))*dt
    cfl_vx = sqrt(this%u%maxabs(domain%mesh_u, domain%parcomm))*dt
    cfl_vy = sqrt(this%v%maxabs(domain%mesh_v, domain%parcomm))*dt

    if(domain%parcomm%myid == 0) &
        print '(A,F15.7,A,2F15.7)', "Gravity waves CFL-max = ", cfl_gv, &
                                    ", Advective CFL-max (x,y) = ", cfl_vx, cfl_vy
end subroutine

subroutine print_integrals_diag(this, state, domain)

    class(swm_output_diag_t), intent(inout) :: this
    class(stvec_t),           intent(inout) :: state
    type(domain_t),           intent(in)    :: domain

    type(grid_field_t), pointer :: h, u, v
    real(kind=8) :: h_min, h_max, mass, TE, PE, KE, Enstrophy

    select type(state)
    class is (stvec_flexible_t)
        call state%get_field(h,"h")
        call state%get_field(u,"u")
        call state%get_field(v,"v")
    class default
        call domain%parcomm%abort("shallow water diagnostics works only with stvec_flexible_t and its extensions")
    end select

    h_min = h%minimum(domain%mesh_p, domain%parcomm)
    h_max = h%maximum(domain%mesh_p, domain%parcomm)
    mass  = this%quadrature_h%mass(h, domain%mesh_p, domain%parcomm)

    call this%calc_energies(te, ke, pe, Enstrophy, h, u, v, domain)

    if(domain%parcomm%myid == 0) &
        print '(7(A,E24.16))', "hmin = ", h_min, " hmax = ", h_max, " mass = ", mass, &
                 " TE = ", te, " KE = ", ke, " PE = ", pe, " Enstrophy = ", Enstrophy

end subroutine print_integrals_diag

subroutine calc_energies(this, te, ke, pe, enstr, h, u, v, domain)

    class(swm_output_diag_t), intent(inout) :: this
    class(grid_field_t),      intent(inout) :: h, u, v
    type(domain_t),           intent(in)    :: domain

    real(kind=8), intent(out) :: te, ke, pe, enstr

    if(this%v_components_type == "covariant") then
        call this%co2contra_op%transform(this%u, this%v, u, v, domain)
        call this%curl_op%calc_curl(this%curl, u, v, domain)
    else
        call this%co2contra_op%transform2co(this%u, this%v, u, v, domain)
        call this%curl_op%calc_curl(this%curl, this%u, this%v, domain)
    end if

    call this%interp_h2uv%interp2d_scalar2vec(this%hu, this%hv, h, domain)
    !For enstrophy:
    if(allocated(this%interp_uv2q)) then
        call this%interp_uv2q%interp2d_vec2vec(this%hqu, this%hqv, this%hu, this%hv, domain)
        call this%hqu%assign(0.5_8, this%hqu, 0.5_8, this%hqv, domain%mesh_q)
    else
        !assume colocated mesh:
        call this%hqu%assign(this%hu,domain%mesh_u)
    end if
    !For energy:
    call this%hu%assign_prod(1.0_8,this%hu,this%u,domain%mesh_u)
    call this%hv%assign_prod(1.0_8,this%hv,this%v,domain%mesh_v)

    call this%hu%assign_prod(0.5_8,this%hu,u,domain%mesh_u)
    call this%hv%assign_prod(0.5_8,this%hv,v,domain%mesh_v)

    call this%h%assign(1.0_8,h,1.0_8,this%h_surf,domain%mesh_p)
    call this%h%assign_prod(0.5_8*Earth_grav,this%h,this%h,domain%mesh_p)

    ke = this%quadrature_u%mass(this%hu,domain%mesh_u,domain%parcomm)+&
         this%quadrature_v%mass(this%hv,domain%mesh_v,domain%parcomm)
    pe   = this%quadrature_h%mass(this%h,domain%mesh_p,domain%parcomm)

    te = ke+pe

    call this%curl%update(1.0_8, this%fcori, domain%mesh_q)

    call this%curl%assign_prod(0.5_8, this%curl, this%curl, domain%mesh_q)
    call this%curl%assign_ratio(1.0_8, this%curl, this%hqu, domain%mesh_q)

    enstr = this%quadrature_q%mass(this%curl, domain%mesh_q, domain%parcomm)

end subroutine

subroutine print_integral_tends(this, operator, state, domain, comment)

    class(swm_output_diag_t), intent(inout) :: this
    class(operator_t),        intent(inout) :: operator
    class(stvec_t),           intent(inout) :: state
    type(domain_t),           intent(in)    :: domain
    character(len=*),         intent(in)    :: comment

    class(stvec_t), allocatable :: tends
    type(grid_field_t), pointer :: h, u, v, h_tend, u_tend, v_tend
    real(kind=8) :: ke_tend, pe_tend, enstrophy_tend

    call state%create_similar(tends)
    call operator%apply(tends, state, domain)

    select type(state)
    class is (stvec_flexible_t)
        call state%get_field(h,"h")
        call state%get_field(u,"u")
        call state%get_field(v,"v")
    class default
        call domain%parcomm%abort("shallow water diagnostics works only with stvec_flexible_t and its extensions")
    end select

    select type(tends)
    class is (stvec_flexible_t)
        call tends%get_field(h_tend,"h")
        call tends%get_field(u_tend,"u")
        call tends%get_field(v_tend,"v")
    class default
        call domain%parcomm%abort("shallow water diagnostics works only with stvec_flexible_t and its extensions")
    end select

    if(this%v_components_type == "covariant") then
        call this%co2contra_op%transform(this%u, this%v, u, v, domain)
    else
        call this%co2contra_op%transform2co(this%u, this%v, u, v, domain)
    end if

    call this%interp_h2uv%interp2d_scalar2vec(this%hu, this%hv, h_tend, domain)
    call this%hu%assign_prod(1.0_8,this%hu,this%u,domain%mesh_u)
    call this%hv%assign_prod(1.0_8,this%hv,this%v,domain%mesh_v)

    call this%hu%assign_prod(0.5_8,this%hu,u,domain%mesh_u)
    call this%hv%assign_prod(0.5_8,this%hv,v,domain%mesh_v)

    ke_tend = this%quadrature_u%mass(this%hu,domain%mesh_u,domain%parcomm)+&
              this%quadrature_v%mass(this%hv,domain%mesh_v,domain%parcomm)

    !below we assume that co2contra is symmetric in terms of uu quadrature:
    call this%interp_h2uv%interp2d_scalar2vec(this%hu, this%hv, h, domain)
    call this%hu%assign_prod(1.0_8,this%hu,this%u,domain%mesh_u)
    call this%hv%assign_prod(1.0_8,this%hv,this%v,domain%mesh_v)

    call this%hu%assign_prod(1.0_8,this%hu,u_tend,domain%mesh_u)
    call this%hv%assign_prod(1.0_8,this%hv,v_tend,domain%mesh_v)

    ke_tend = this%quadrature_u%mass(this%hu,domain%mesh_u,domain%parcomm)+&
              this%quadrature_v%mass(this%hv,domain%mesh_v,domain%parcomm)+&
              ke_tend
    
    call this%h%assign(1.0_8,h,1.0_8,this%h_surf,domain%mesh_p)
    call this%h%assign_prod(Earth_grav,this%h,h_tend,domain%mesh_p)
    pe_tend = this%quadrature_h%mass(this%h,domain%mesh_p,domain%parcomm)

    !Enstrophy tend
    call this%interp_h2uv%interp2d_scalar2vec(this%hu, this%hv, h, domain)
    if(allocated(this%interp_uv2q)) then
        call this%interp_uv2q%interp2d_vec2vec(this%hq, this%hqu, this%hu, this%hv, domain)
    else
        call this%hq%assign(this%hu,domain%mesh_q)
    end if

    call this%interp_h2uv%interp2d_scalar2vec(this%hu, this%hv, h_tend, domain)
    if(allocated(this%interp_uv2q)) then
        call this%interp_uv2q%interp2d_vec2vec(this%hq_tend, this%hqu, this%hu, this%hv, domain)
    else
        call this%hq_tend%assign(this%hu,domain%mesh_q)
    end if

    if(this%v_components_type == "covariant") then
        call this%curl_op%calc_curl(this%curl,u,v,domain)
        call this%curl_op%calc_curl(this%curl_tend,u_tend,v_tend,domain)
    else
        call this%curl_op%calc_curl(this%curl,this%u,this%v,domain)
        call this%co2contra_op%transform2co(this%u,this%v,u_tend,v_tend,domain)
        call this%curl_op%calc_curl(this%curl_tend,this%u,this%v,domain)
    end if
    call this%curl%update(1.0_8,this%fcori,domain%mesh_q)

    call this%curl_tend%assign_prod(1.0_8,this%curl,this%curl_tend,domain%mesh_q)
    call this%curl_tend%assign_ratio(1.0_8,this%curl_tend,this%hq,domain%mesh_q)

    call this%curl%assign_prod(-0.5_8,this%curl,this%curl,domain%mesh_q)
    call this%curl%assign_prod(1.0_8,this%curl,this%hq_tend,domain%mesh_q)
    call this%curl%assign_ratio(1.0_8,this%curl,this%hq,domain%mesh_q)
    call this%curl%assign_ratio(1.0_8,this%curl,this%hq,domain%mesh_q)

    enstrophy_tend = &
              this%quadrature_q%mass(this%curl_tend,domain%mesh_q,domain%parcomm)+&
              this%quadrature_q%mass(this%curl,domain%mesh_q,domain%parcomm)

    if(domain%parcomm%myid == 0) &
        print '(2(A,E20.12))', comment // ": TE_tend = ", ke_tend+pe_tend, &
                 ", Enstrophy_tend = ", enstrophy_tend

end subroutine print_integral_tends

end module swm_output_diag_mod