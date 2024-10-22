module operator_skew_swm_mod

use stvec_mod,      only : stvec_t
use domain_mod,     only : domain_t
use operator_mod,   only : operator_t
use grid_field_mod, only : grid_field_t, tile_field_t
use mesh_mod,       only : mesh_t, tile_mesh_t
  
use abstract_div_mod,        only : div_operator_t
use abstract_flux_div_mod,   only : flux_div_operator_t
use abstract_grad_mod,       only : grad_operator_t
use abstract_coriolis_mod,   only : coriolis_operator_t
use abstract_massflux_mod,   only : massflux_operator_t
use abstract_co2contra_mod,  only : co2contra_operator_t
use abstract_quadrature_mod, only : quadrature_t

use abstract_interpolators2d_mod, only : interpolator2d_vec2vec_t
use shallow_atm_trajectory_tools_mod, only : calc_shallow_atm_xyz_wind, calc_hor_native_wind

use stvec_swm_mod, only : stvec_swm_t
use key_value_mod, only : key_value_r8_t
 
implicit none

type, public, extends(operator_t) :: operator_skew_swm_t
    
    character(:),                    allocatable :: v_components_type
    class (co2contra_operator_t),    allocatable :: co2contra_op
    class(flux_div_operator_t),      allocatable :: h_flux_div
    class(div_operator_t),           allocatable :: div_op
    class(grad_operator_t),          allocatable :: grad_op
    class(coriolis_operator_t),      allocatable :: coriolis_op
    class(massflux_operator_t),      allocatable :: massflux_op
    class(interpolator2d_vec2vec_t), allocatable :: v2p_interp, p2v_interp
    
    type(grid_field_t) :: hx, hy, up, vp
    type(grid_field_t) :: vx, vy, vz
    type(grid_field_t) :: vx_tend_adv, vy_tend_adv, vz_tend_adv
    type(grid_field_t) :: vx_tend_metric, vy_tend_metric, vz_tend_metric
    type(grid_field_t) :: div, p_buff
    ! type(grid_field_t) :: h_surf !orography
    ! type(grid_field_t) :: h_total !h+h_surf
    ! type(grid_field_t) :: div, grad_x, grad_y, curl
    ! type(grid_field_t) :: cor_u, cor_v
    ! type(grid_field_t) :: KE !kinetic energy
    ! type(grid_field_t) :: ut, vt !contravariant components
    ! type(grid_field_t) :: hu, hv !mass fluxes in continuty eq.
    
    real(kind=8) :: grav
    
    contains
    procedure, public :: apply
    procedure, public :: get_diagnostics
    ! procedure, public :: calc_energy

end type operator_skew_swm_t

contains

subroutine apply(this, vout, vin, domain)
    class(operator_skew_swm_t), intent(inout) :: this
    class(stvec_t),             intent(inout) :: vout !inout to enable preallocated vectors
    class(stvec_t),             intent(inout) :: vin
    type(domain_t),             intent(in)    :: domain

    call vout%assign(0.0_8, domain)

    select type (vin)
    type is (stvec_swm_t)
    select type (vout)
    type is (stvec_swm_t)

        call this%grad_op%calc_grad(this%hx, this%hy, vin%h, domain)
        call this%co2contra_op%transform(vout%u, vout%v, this%hx, this%hy, domain)
        call this%coriolis_op%calc_coriolis_contra(this%hx, this%hy, vin%u, vin%v, domain)

        call vout%u%assign(-this%grav,vout%u,1.0_8,this%hx,domain%mesh_u)
        call vout%v%assign(-this%grav,vout%v,1.0_8,this%hy,domain%mesh_v)

        call this%div_op%calc_div(this%div,vin%u, vin%v, domain)

        call this%v2p_interp%interp2d_vec2vec(this%up, this%vp, vin%u, vin%v, domain)
        call calc_shallow_atm_xyz_wind(this%vx, this%vy,this%vz,this%up, this%vp, domain%mesh_p)
 
        ! call this%grad_op%calc_grad(this%hx,this%hy,this%vx,domain)
        ! call this%hx%assign_prod(1.0_8, this%hx,vin%u,domain%mesh_u)
        ! call this%hy%assign_prod(1.0_8, this%hy,vin%v,domain%mesh_v)
        ! call this%v2p_interp%interp2d_vec2vec(this%up,this%vp,this%hx,this%hy,domain)
        ! call this%vx_tend_adv%assign(1.0_8,this%up,1.0_8,this%vp,domain%mesh_p)

        ! call this%grad_op%calc_grad(this%hx,this%hy,this%vy,domain)
        ! call this%hx%assign_prod(1.0_8, this%hx,vin%u,domain%mesh_u)
        ! call this%hy%assign_prod(1.0_8, this%hy,vin%v,domain%mesh_v)
        ! call this%v2p_interp%interp2d_vec2vec(this%up,this%vp,this%hx,this%hy,domain)
        ! call this%vy_tend_adv%assign(1.0_8,this%up,1.0_8,this%vp,domain%mesh_p)

        ! call this%grad_op%calc_grad(this%hx,this%hy,this%vz,domain)
        ! call this%hx%assign_prod(1.0_8, this%hx,vin%u,domain%mesh_u)
        ! call this%hy%assign_prod(1.0_8, this%hy,vin%v,domain%mesh_v)
        ! call this%v2p_interp%interp2d_vec2vec(this%up,this%vp,this%hx,this%hy,domain)
        ! call this%vz_tend_adv%assign(1.0_8,this%up,1.0_8,this%vp,domain%mesh_p)

        call this%h_flux_div%calc_flux_div(this%vx_tend_adv, this%vx, vin%u, vin%v, domain)
        call this%p_buff%assign_prod( 1.0_8, this%div, this%vx, domain%mesh_p)
        call this%vx_tend_adv%update(-1.0_8, this%p_buff, domain%mesh_p)

        call this%h_flux_div%calc_flux_div(this%vy_tend_adv, this%vy, vin%u, vin%v, domain)
        call this%p_buff%assign_prod( 1.0_8, this%div, this%vy, domain%mesh_p)
        call this%vy_tend_adv%update(-1.0_8, this%p_buff, domain%mesh_p)

        call this%h_flux_div%calc_flux_div(this%vz_tend_adv, this%vz, vin%u, vin%v, domain)
        call this%p_buff%assign_prod( 1.0_8, this%div, this%vz, domain%mesh_p)
        call this%vz_tend_adv%update(-1.0_8, this%p_buff, domain%mesh_p)

        call calc_shallow_atm_metric_force(this%vx_tend_metric, this%vy_tend_metric, this%vz_tend_metric, this%vx, this%vy, this%vz, domain%mesh_p)

        call this%vx_tend_adv%update(-1.0_8,this%vx_tend_metric, domain%mesh_p)
        call this%vy_tend_adv%update(-1.0_8,this%vy_tend_metric, domain%mesh_p)
        call this%vz_tend_adv%update(-1.0_8,this%vz_tend_metric, domain%mesh_p)

        call calc_hor_native_wind(this%up, this%vp, this%vx_tend_adv, this%vy_tend_adv, this%vz_tend_adv, domain%mesh_p)
        call this%p2v_interp%interp2d_vec2vec(this%hx,this%hy, this%up, this%vp, domain)
        call vout%u%update(-1.0_8,this%hx,domain%mesh_u)
        call vout%v%update(-1.0_8,this%hy,domain%mesh_v)

        call this%h_flux_div%calc_flux_div(vout%h, vin%h, vin%u, vin%v, domain)
        call vout%h%assign(-1.0_8,vout%h,domain%mesh_p)
    class default
        call domain%parcomm%abort("operator_skew_swm_t: vout of wrong type")
    end select
    class default
        call domain%parcomm%abort("operator_skew_swm_t: vin of wrong type")
    end select

end subroutine

function get_diagnostics(this, v, domain) result(diagnostics)

    use key_value_mod, only : key_value_r8_t

    class(operator_skew_swm_t), intent(inout) :: this
    class(stvec_t),             intent(inout) :: v
    type(domain_t),             intent(in)    :: domain

    type(key_value_r8_t)  :: diagnostics

    integer(kind=4), parameter :: ndiag = 3
    real(kind=8) :: te, ke, pe, enstrophy

    allocate(diagnostics%keys(ndiag))
    allocate(diagnostics%values(ndiag))

    select type(v)
    type is (stvec_swm_t)

        diagnostics%keys(1)%str = "hmin"
        diagnostics%values(1) = v%h%minimum(domain%mesh_p, domain%parcomm)
        diagnostics%keys(2)%str = "hmax"
        diagnostics%values(2) = v%h%maximum(domain%mesh_p, domain%parcomm)
        diagnostics%keys(3)%str = "mass"
        diagnostics%values(3) = 0.0_8!this%quadrature_h%mass(v%h, domain%mesh_p, domain%parcomm)

    class default
        call domain%parcomm%abort("wrong type of v in operator_swm_t get_diagnostics")
    end select

end function get_diagnostics

subroutine calc_shallow_atm_metric_force(fx, fy, fz, vx, vy, vz, mesh)

    type(grid_field_t), intent(inout) :: fx, fy, fz
    type(grid_field_t), intent(in)    :: vx, vy, vz
    type(mesh_t),       intent(in)    :: mesh

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        call calc_shallow_atm_metric_force_tile(fx%tile(t), fy%tile(t), fz%tile(t), &
                                                vx%tile(t), vy%tile(t), vz%tile(t), &
                                                mesh%tile(t), mesh%scale)
    end do

end subroutine

subroutine calc_shallow_atm_metric_force_tile(fx, fy, fz, vx, vy, vz,&
                                              mesh, scale)

    type(tile_field_t), intent(inout) :: fx, fy, fz
    type(tile_field_t), intent(in)    :: vx, vy, vz
    type(tile_mesh_t),  intent(in)    :: mesh
    real(kind=8),       intent(in)    :: scale

    integer(kind=4) :: i, j, k
    real(kind=8)    :: f

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                f = (vx%p(i,j,k)**2+vy%p(i,j,k)**2+vz%p(i,j,k)**2) / scale
                fx%p(i,j,k) = -f*mesh%rx(i,j,k)
                fy%p(i,j,k) = -f*mesh%ry(i,j,k)
                fz%p(i,j,k) = -f*mesh%rz(i,j,k)
            end do
        end do
    end do

end subroutine


end module operator_skew_swm_mod