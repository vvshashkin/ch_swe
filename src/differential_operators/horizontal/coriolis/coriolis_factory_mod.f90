module coriolis_factory_mod

use domain_mod,             only : domain_t
use abstract_coriolis_mod,  only : coriolis_operator_t
use parcomm_mod,            only : parcomm_global
use grid_field_mod,         only : grid_field_t

implicit none

private

public :: create_coriolis, calc_coriolis_parameter

contains

subroutine create_coriolis(coriolis_op, coriolis_op_name, domain)

    class(coriolis_operator_t), allocatable, intent(out) :: coriolis_op
    character(len=*),                        intent(in)  :: coriolis_op_name
    type(domain_t),                          intent(in)  :: domain

    select case(coriolis_op_name)

    case("coriolis_colocated")
        call create_coriolis_colocated(coriolis_op, domain)
    case("coriolis_Cgrid_sbp21")
        call create_coriolis_Cgrid_sbp(coriolis_op, "W21_stagered_interp_i2c", &
                                          "interp2d_uv2pvec_C_sbp21", &
                                          "interp2d_pvec2uv_C_sbp21", domain)
    case("coriolis_Cgrid_sbp42")
        call create_coriolis_Cgrid_sbp(coriolis_op, "W42_stagered_interp_i2c", &
                                       "interp2d_uv2pvec_C_sbp42", &
                                       "interp2d_pvec2uv_C_sbp42", domain)
    case("coriolis_Cgrid_sbp63")
        call create_coriolis_Cgrid_sbp(coriolis_op, "W63_stagered_interp_i2c", &
                                       "interp2d_uv2pvec_C_sbp63", &
                                       "interp2d_pvec2uv_C_sbp63", domain)
    case("coriolis_Cgrid_noncons_sbp21")
        call create_coriolis_Cgrid_noncons_sbp(coriolis_op, "interp2d_uv2pvec_C_sbp21", &
                                          "interp2d_pvec2uv_C_sbp21", "co2contra_c_sbp21_new", domain)
    case("coriolis_Cgrid_noncons_sbp42")
        call create_coriolis_Cgrid_noncons_sbp(coriolis_op, "interp2d_uv2pvec_C_sbp42", &
                                          "interp2d_pvec2uv_C_sbp42", "co2contra_c_sbp42_new", domain)
    case("coriolis_Cgrid_noncons_sbp63")
        call create_coriolis_Cgrid_noncons_sbp(coriolis_op, "interp2d_uv2pvec_C_sbp63", &
                                          "interp2d_pvec2uv_C_sbp63", "co2contra_c_sbp63_new", domain)
    case("coriolis_Ch_sbp21", "coriolis_Ch_sbp21_opt2")
        call create_coriolis_Ch_sbp(coriolis_op, coriolis_op_name, &
                                          "interp2d_uv2pvec_Ch_sbp21", &
                                          "interp2d_pvec2uv_Ch_sbp21", &
                                          "interp2d_qvec2uv_Ch_sbp21", domain)
    case("coriolis_Ch_sbp42", "coriolis_Ch_sbp42_opt2")
        call create_coriolis_Ch_sbp(coriolis_op, coriolis_op_name, &
                                          "interp2d_uv2pvec_Ch_sbp42", &
                                          "interp2d_pvec2uv_Ch_sbp42", &
                                          "interp2d_qvec2uv_Ch_sbp42", domain)
    case("coriolis_Ch_sbp63", "coriolis_Ch_sbp63_opt2")
        call create_coriolis_Ch_sbp(coriolis_op, coriolis_op_name, &
                                          "interp2d_uv2pvec_Ch_sbp63", &
                                          "interp2d_pvec2uv_Ch_sbp63", &
                                          "interp2d_qvec2uv_Ch_sbp63", domain)
    case default
        call parcomm_global%abort("Unknown coriolis operator: "//coriolis_op_name)
    end select

end subroutine create_coriolis

subroutine create_coriolis_Cgrid_sbp(coriolis_op,sbp_i2c_interp_name, v2h_interp_name, &
                                     h2v_interp_name, domain)

    use coriolis_Cgrid_mod,     only : coriolis_Cgrid_t
    use grid_field_factory_mod, only : create_grid_field

    use interpolator2d_factory_mod,   only : create_vec2vec_interpolator2d
    use interpolator_w2h_factory_mod, only : create_w2h_interpolator

    class(coriolis_operator_t), allocatable, intent(out) :: coriolis_op
    character(len=*),                        intent(in)  :: sbp_i2c_interp_name, &
                                                            h2v_interp_name, v2h_interp_name
    type(domain_t),                          intent(in)  :: domain

    type(coriolis_Cgrid_t), allocatable :: cor_Cgrid

    integer(kind=4) :: halo_width

    !WORKAROUND
    halo_width = 5

    allocate(cor_Cgrid)

    call create_grid_field(cor_Cgrid%f, 0, 0, domain%mesh_p)
    call calc_coriolis_parameter(cor_Cgrid%f, domain%mesh_p, domain%metric)

    call create_grid_field(cor_Cgrid%Ghu, halo_width+1, 0, domain%mesh_u)
    call create_grid_field(cor_Cgrid%Ghv, halo_width+1, 0, domain%mesh_v)

    call create_grid_field(cor_Cgrid%Ghu_p, 0, 0, domain%mesh_p)
    call create_grid_field(cor_Cgrid%Ghv_p, 0, 0, domain%mesh_p)

    call create_grid_field(cor_Cgrid%curl_p,  0, 0, domain%mesh_p)

    call create_grid_field(cor_Cgrid%cor_u_p, halo_width+1, 0, domain%mesh_p)
    call create_grid_field(cor_Cgrid%cor_v_p, halo_width+1, 0, domain%mesh_p)


    call create_vec2vec_interpolator2d(cor_Cgrid%interp_v2h_op, v2h_interp_name, domain)
    call create_vec2vec_interpolator2d(cor_Cgrid%interp_h2v_op, h2v_interp_name, domain)
    call create_w2h_interpolator(cor_Cgrid%interp_w2h_op, sbp_i2c_interp_name, domain)

    call move_alloc(cor_Cgrid, coriolis_op)

end subroutine create_coriolis_Cgrid_sbp

subroutine create_coriolis_Cgrid_noncons_sbp(coriolis_op, sbp_i2c_interp_name, &
                                             sbp_c2i_interp_name, co2contra_op_name, domain)

    use coriolis_Cgrid_noncons_mod,   only : coriolis_Cgrid_noncons_t
    use sbp_factory_mod,              only : create_sbp_operator
    use grid_field_factory_mod,       only : create_grid_field
    use co2contra_factory_mod,        only : create_co2contra_operator
    use interpolator2d_factory_mod,   only : create_vec2vec_interpolator2d

    class(coriolis_operator_t), allocatable, intent(out) :: coriolis_op
    character(len=*),                        intent(in)  :: sbp_i2c_interp_name, &
                                                            sbp_c2i_interp_name, &
                                                            co2contra_op_name
    type(domain_t),                          intent(in)  :: domain

    type(coriolis_Cgrid_noncons_t), allocatable :: cor_Cgrid

    integer(kind=4) :: halo_width

    !WORKAROUND
    halo_width = 5

    allocate(cor_Cgrid)

    call create_grid_field(cor_Cgrid%f, 0, 0, domain%mesh_p)
    call calc_coriolis_parameter(cor_Cgrid%f, domain%mesh_p, domain%metric)
    call create_grid_field(cor_Cgrid%u,halo_width+1, 0, domain%mesh_u)
    call create_grid_field(cor_Cgrid%v,halo_width+1, 0, domain%mesh_v)
    call create_grid_field(cor_Cgrid%uh,halo_width+1, 0, domain%mesh_p)
    call create_grid_field(cor_Cgrid%vh,halo_width+1, 0, domain%mesh_p)
    cor_Cgrid%co2contra = create_co2contra_operator(domain, co2contra_op_name)

    ! call create_v2h_interpolator(cor_Cgrid%interp_v2h_op, sbp_i2c_interp_name, domain)
    call create_vec2vec_interpolator2d(cor_Cgrid%interp_v2h_op, sbp_i2c_interp_name, domain)
    call create_vec2vec_interpolator2d(cor_Cgrid%interp_h2v_op, sbp_c2i_interp_name, domain)

    call move_alloc(cor_Cgrid, coriolis_op)

end subroutine create_coriolis_Cgrid_noncons_sbp

subroutine create_coriolis_Ch_sbp(coriolis_op, coriolis_op_name, &
                                  v2h_interp_name, h2v_interp_name, &
                                  q2v_interp_name, domain)

use coriolis_Ch_mod,              only : coriolis_Ch_t, coriolis_CH_opt2_t
use grid_field_factory_mod,       only : create_grid_field
use vec_math_mod,                 only : multiply_by_J_self
use interpolator2d_factory_mod,   only : create_vec2vec_interpolator2d
use halo_factory_mod,             only : create_vector_halo_procedure
use exchange_factory_mod,         only : create_symmetric_halo_vec_exchange_Ch

class(coriolis_operator_t), allocatable, intent(out) :: coriolis_op
character(len=*),                        intent(in)  :: coriolis_op_name
character(len=*),                        intent(in)  :: h2v_interp_name
character(len=*),                        intent(in)  :: v2h_interp_name
character(len=*),                        intent(in)  :: q2v_interp_name
type(domain_t),                          intent(in)  :: domain

class(coriolis_Ch_t), allocatable :: cor_Ch
integer(kind=4), parameter :: halo_width = 5
logical :: is_opt2

is_opt2 = .false.
if(coriolis_op_name(len(coriolis_op_name)-3:len(coriolis_op_name))=="opt2") then
    is_opt2 = .true.
end if

!WORKAROUND
if(is_opt2) then
    allocate(coriolis_Ch_opt2_t ::cor_Ch)
else
    allocate(coriolis_Ch_t ::cor_Ch)
end if

call create_grid_field(cor_Ch%fJ2, 0, 0, domain%mesh_p)
call cor_Ch%fJ2%assign(0.0_8, domain%mesh_p)
call calc_coriolis_parameter(cor_Ch%fJ2, domain%mesh_p, domain%metric)
if(.not. is_opt2) then
    call multiply_by_J_self(cor_Ch%fJ2,domain%mesh_p)
    call multiply_by_J_self(cor_Ch%fJ2,domain%mesh_p)
end if

call create_grid_field(cor_Ch%up, halo_width, 0, domain%mesh_p)
call create_grid_field(cor_Ch%vp, halo_width, 0, domain%mesh_p)
call create_grid_field(cor_Ch%curl_u, halo_width, 0, domain%mesh_u)
call create_grid_field(cor_Ch%curl_v, halo_width, 0, domain%mesh_v)
call create_grid_field(cor_Ch%curl_pu, halo_width, 0, domain%mesh_p)
call create_grid_field(cor_Ch%curl_pv, halo_width, 0, domain%mesh_p)

call create_vec2vec_interpolator2d(cor_Ch%interp_v2h_op, v2h_interp_name, domain)
call create_vec2vec_interpolator2d(cor_Ch%interp_h2v_op, h2v_interp_name, domain)
call create_vec2vec_interpolator2d(cor_Ch%interp_q2v_op, q2v_interp_name, domain)

if(.not. is_opt2) then
    call create_vector_halo_procedure(cor_Ch%uv_sync, domain, 1, "ecs_Ah_vec_sync_covariant")
else
    call create_vector_halo_procedure(cor_Ch%uv_sync, domain, 1, "ecs_Ah_vec_sync_contra")
end if

cor_Ch%CH_uv_exchange = create_symmetric_halo_vec_exchange_Ch(domain%partition, domain%parcomm, domain%topology, 1, "cross")

call move_alloc(cor_Ch, coriolis_op)

end subroutine create_coriolis_Ch_sbp

subroutine create_coriolis_colocated(coriolis_op, domain)

    use coriolis_colocated_mod, only : coriolis_colocated_t
    use grid_field_factory_mod, only : create_grid_field

    class(coriolis_operator_t), allocatable, intent(out) :: coriolis_op
    type(domain_t),                          intent(in)  :: domain

    type(coriolis_colocated_t), allocatable :: coriolis_colocated

    allocate(coriolis_colocated)

    call create_grid_field(coriolis_colocated%f, 0, 0, domain%mesh_u)

    call calc_coriolis_parameter(coriolis_colocated%f, domain%mesh_p, domain%metric)

    call move_alloc(coriolis_colocated, coriolis_op)

end subroutine create_coriolis_colocated

subroutine calc_coriolis_parameter(f, mesh, metric)

    use metric_mod,     only : metric_t
    use mesh_mod,       only : mesh_t
    use sph_coords_mod, only : cart2sph

    type(grid_field_t), intent(inout) :: f
    type(mesh_t),       intent(in)  :: mesh
    class(metric_t),    intent(in)  :: metric

    integer(kind=4) :: t, i, j, k
    real(kind=8) :: pr, n(3), lam, phi

    if(metric%coriolis_const) then

        do t = mesh%ts, mesh%te
            do k = mesh%tile(t)%ks, mesh%tile(t)%ke
                do j = mesh%tile(t)%js, mesh%tile(t)%je
                    do i = mesh%tile(t)%is, mesh%tile(t)%ie
                        f%tile(t)%p(i,j,k) = mesh%omega
                    end do
                end do
            end do
        end do

    else

        do t = mesh%ts, mesh%te
            do k = mesh%tile(t)%ks, mesh%tile(t)%ke
                do j = mesh%tile(t)%js, mesh%tile(t)%je
                    do i = mesh%tile(t)%is, mesh%tile(t)%ie
                        n = [mesh%tile(t)%rx(i,j,1),mesh%tile(t)%ry(i,j,1),mesh%tile(t)%rz(i,j,1)]
                        pr = dot_product(mesh%rotation_axis, n) / sqrt(dot_product(n,n))
                        f%tile(t)%p(i,j,k) = 2*mesh%omega*pr
                    end do
                end do
            end do
        end do

    end if

end subroutine calc_coriolis_parameter

end module coriolis_factory_mod
