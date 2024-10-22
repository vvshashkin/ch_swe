module coriolis_3d_factory_mod

use abstract_coriolis_3d_mod,           only : coriolis3d_operator_t
use shallow_atm_coriolis_colocated_mod, only : shallow_atm_coriolis_colocated_t
use shallow_atm_coriolis_C_mod,         only : shallow_atm_coriolis_C_t
use domain_mod,                         only : domain_t
use grid_field_factory_mod,             only : create_grid_field, create_grid_field_2d
use interpolator2d_factory_mod,         only : create_vec2vec_interpolator2d
use parcomm_mod,                        only : parcomm_global

implicit none

contains

subroutine create_coriolis_3d_operator(coriolis_op,coriolis_op_name,domain)
    class(coriolis3d_operator_t), allocatable, intent(out) :: coriolis_op
    character(len=*),                          intent(in)  :: coriolis_op_name
    type(domain_t),                            intent(in)  :: domain

    select case (coriolis_op_name)
    case("shallow_atm_colocated_coriolis")
        call create_shallow_atm_coriolis_colocated(coriolis_op,coriolis_op_name,domain)
    case("shallow_atm_C_sbp21_corriolis", "shallow_atm_C_sbp42_corriolis")
        call create_shallow_atm_coriolis_C_sbp(coriolis_op,coriolis_op_name,domain)
    case default
        call parcomm_global%abort("create_coriolis_3d_operator, unknown coriolis op name"// coriolis_op_name)
    end select
end subroutine

subroutine create_shallow_atm_coriolis_colocated(coriolis_op,coriolis_op_name,domain)
    class(coriolis3d_operator_t), allocatable, intent(out) :: coriolis_op
    character(len=*),                          intent(in)  :: coriolis_op_name
    type(domain_t),                            intent(in)  :: domain

    type(shallow_atm_coriolis_colocated_t), allocatable :: coriolis
    integer(kind=4) :: i, j, k, t, pind
    real(kind=8)    :: alpha, beta, cori_f, Q(6), Jac, omega, axis(3), r(3)

    allocate(coriolis)

    call create_grid_field(coriolis%u2u, 0, 0, domain%mesh_u)
    call create_grid_field(coriolis%u2v, 0, 0, domain%mesh_u)
    call create_grid_field(coriolis%v2u, 0, 0, domain%mesh_u)
    call create_grid_field(coriolis%v2v, 0, 0, domain%mesh_u)

    omega = domain%mesh_u%omega
    axis(1:3) = domain%mesh_u%rotation_axis(1:3)

    do t = domain%mesh_u%ts, domain%mesh_u%te
        pind = domain%mesh_u%tile(t)%panel_ind
        do k = domain%mesh_u%tile(t)%ks, domain%mesh_u%tile(t)%ke
            do j = domain%mesh_u%tile(t)%js, domain%mesh_u%tile(t)%je
                beta = domain%mesh_u%tile(t)%get_beta(j)
                do i = domain%mesh_u%tile(t)%is, domain%mesh_u%tile(t)%ie
                    alpha = domain%mesh_u%tile(t)%get_alpha(i)
                    Q(1:6) = domain%metric%calculate_Q_2d(pind,alpha,beta)
                    Jac = domain%metric%calculate_J_2d(pind,alpha,beta)
                    r(1:3) = [domain%mesh_u%tile(t)%rx(i,j,k),domain%mesh_u%tile(t)%ry(i,j,k),&
                              domain%mesh_u%tile(t)%rz(i,j,k)]
                    cori_f = 2.0_8*omega*sum(axis(1:3)*r(1:3))

                    coriolis%u2u%tile(t)%p(i,j,k) = cori_f*Q(2)/Jac
                    coriolis%u2v%tile(t)%p(i,j,k) =-cori_f*Q(1)/Jac
                    coriolis%v2u%tile(t)%p(i,j,k) = cori_f*Q(3)/Jac
                    coriolis%v2v%tile(t)%p(i,j,k) =-cori_f*Q(2)/Jac
                end do
            end do
        end do
    end do

    call move_alloc(coriolis,coriolis_op)

end subroutine

subroutine create_shallow_atm_coriolis_C_sbp(coriolis_op,coriolis_op_name,domain)
    class(coriolis3d_operator_t), allocatable, intent(out) :: coriolis_op
    character(len=*),                          intent(in)  :: coriolis_op_name
    type(domain_t),                            intent(in)  :: domain

    type(shallow_atm_coriolis_C_t), allocatable :: coriolis
    integer(kind=4) :: i, j, k, t, pind, halo_width
    real(kind=8)    :: alpha, beta, cori_f, Q(6), Jac, omega, axis(3), r(3)
    character(len=:), allocatable :: pvec2uv_interp_name, uv2pvec_interp_name

    allocate(coriolis)

    select case(coriolis_op_name)
    case("shallow_atm_C_sbp21_corriolis")
        uv2pvec_interp_name = "interp2d_uv2pvec_C_sbp21"
        pvec2uv_interp_name = "interp2d_pvec2uv_C_sbp21"
        halo_width = 2
    case("shallow_atm_C_sbp42_corriolis")
        uv2pvec_interp_name = "interp2d_uv2pvec_C_sbp42"
        pvec2uv_interp_name = "interp2d_pvec2uv_C_sbp42"
        halo_width = 4
    case default
        call parcomm_global%abort("unknown shallow atmosphere C-staggere coriolis operator "//coriolis_op_name)
    end select

    call create_vec2vec_interpolator2d(coriolis%interp2d_uv2p, uv2pvec_interp_name, domain)
    call create_vec2vec_interpolator2d(coriolis%interp2d_p2uv, pvec2uv_interp_name, domain)

    call create_grid_field_2d(coriolis%Q2d_uu, 0, domain%mesh_u)
    call create_grid_field_2d(coriolis%Ju2d  , 0, domain%mesh_u)
    call create_grid_field_2d(coriolis%Q2d_vv, 0, domain%mesh_v)
    call create_grid_field_2d(coriolis%Jv2d  , 0, domain%mesh_v)
    call create_grid_field_2d(coriolis%Q2d_uv, 0, domain%mesh_p)
    call create_grid_field_2d(coriolis%Jp2d  , 0, domain%mesh_p)
    call create_grid_field_2d(coriolis%fcori_p,0, domain%mesh_p)

    call create_grid_field(coriolis%up,   halo_width, 0, domain%mesh_p)
    call create_grid_field(coriolis%vp,   halo_width, 0, domain%mesh_p)
    call create_grid_field(coriolis%u_cov,halo_width, 0, domain%mesh_u)
    call create_grid_field(coriolis%v_cov,halo_width, 0, domain%mesh_v)

    omega = domain%mesh_u%omega
    axis(1:3) = domain%mesh_u%rotation_axis(1:3)

    do t = domain%mesh_u%ts, domain%mesh_u%te
        pind = domain%mesh_u%tile(t)%panel_ind
        do j = domain%mesh_u%tile(t)%js, domain%mesh_u%tile(t)%je
            beta = domain%mesh_u%tile(t)%get_beta(j)
            do i = domain%mesh_u%tile(t)%is, domain%mesh_u%tile(t)%ie
                alpha = domain%mesh_u%tile(t)%get_alpha(i)
                Q(1:6) = domain%metric%calculate_Q_2d(pind,alpha,beta)
                coriolis%Q2d_uu%tile(t)%p(i,j,1) = Q(1)
                coriolis%Ju2d%tile(t)%p(i,j,1) = &
                               domain%metric%calculate_J_2d(pind,alpha,beta)
            end do
        end do
    end do

    do t = domain%mesh_v%ts, domain%mesh_v%te
        pind = domain%mesh_v%tile(t)%panel_ind
        do j = domain%mesh_v%tile(t)%js, domain%mesh_v%tile(t)%je
            beta = domain%mesh_v%tile(t)%get_beta(j)
            do i = domain%mesh_v%tile(t)%is, domain%mesh_v%tile(t)%ie
                alpha = domain%mesh_v%tile(t)%get_alpha(i)
                Q(1:6) = domain%metric%calculate_Q_2d(pind,alpha,beta)
                coriolis%Q2d_vv%tile(t)%p(i,j,1) = Q(3)
                coriolis%Jv2d%tile(t)%p(i,j,1) = &
                               domain%metric%calculate_J_2d(pind,alpha,beta)
            end do
        end do
    end do

    do t = domain%mesh_p%ts, domain%mesh_p%te
        pind = domain%mesh_p%tile(t)%panel_ind
        do j = domain%mesh_p%tile(t)%js, domain%mesh_p%tile(t)%je
            beta = domain%mesh_p%tile(t)%get_beta(j)
            do i = domain%mesh_p%tile(t)%is, domain%mesh_p%tile(t)%ie
                alpha = domain%mesh_p%tile(t)%get_alpha(i)
                Q(1:6) = domain%metric%calculate_Q_2d(pind,alpha,beta)
                coriolis%Q2d_uv%tile(t)%p(i,j,1) = Q(2)
                coriolis%Jp2d%tile(t)%p(i,j,1) = &
                                   domain%metric%calculate_J_2d(pind,alpha,beta)
                r(1:3) = [domain%mesh_p%tile(t)%rx(i,j,1),domain%mesh_p%tile(t)%ry(i,j,1),&
                          domain%mesh_p%tile(t)%rz(i,j,1)]
                coriolis%fcori_p%tile(t)%p(i,j,1) = 2.0_8*omega*sum(axis(1:3)*r(1:3))
            end do
        end do
    end do

    call move_alloc(coriolis,coriolis_op)

end subroutine

end module
