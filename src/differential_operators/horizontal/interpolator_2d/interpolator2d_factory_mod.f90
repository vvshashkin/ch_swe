module interpolator2d_factory_mod

use abstract_interpolators2d_mod,    only : interpolator2d_scalar2vec_t, &
                                            interpolator2d_vec2vec_t
use interpolator_p2uv_sbp_C_mod,     only : interpolator_p2uv_sbp_C_t,   &
                                            interpolator_pvec2uv_sbp_C_t,&
                                            interpolator_pvec2uv_sbp_Ch_t
use p2uv_colocated_mod,              only : p2uv_colocated_t
use interpolator_uv2p_mod,           only : interpolator2d_uv2p_sbp_C_t
use interpolator_q2uv_sbp_mod,       only : interpolator_q2uv_sbp_Ch_t
use interpolator_uv2q_sbp_mod,       only : interpolator_uv2q_sbp_t
use domain_mod,                      only : domain_t
use sbp_factory_mod,                 only : create_sbp_operator
use exchange_factory_mod,            only : create_o_points_halo_exchange,  &
                                            create_z_points_halo_exchange,  &
                                            create_xy_points_halo_exchange, &
                                            create_symmetric_halo_vec_exchange_C, &
                                            create_symmetric_halo_vec_exchange_Ch, &
                                            create_symmetric_halo_vec_exchange_Cz
use parcomm_mod,                     only : parcomm_global

use sbp_interp_mod,                  only : sbp_interp_t
use sbp_interp_c2i_21_mod,           only : sbp_interp_c2i_21_t
use sbp_interp_i2c_21_mod,           only : sbp_interp_i2c_21_t
use sbp_interp_i2c_42_mod,           only : sbp_interp_i2c_42_t
use sbp_interp_c2i_42_mod,           only : sbp_interp_c2i_42_t
use sbp_interp_i2c_63_mod,           only : sbp_interp_i2c_63_t
use sbp_interp_c2i_63_mod,           only : sbp_interp_c2i_63_t

implicit none

contains

subroutine create_scalar2vec_interpolator2d(interpolator2d, interp2d_name, domain)
    class(interpolator2d_scalar2vec_t), allocatable, intent(out) :: interpolator2d
    character(len=*),                                intent(in)  :: interp2d_name
    type(domain_t),                                  intent(in)  :: domain

    select case(interp2d_name)
    case("interp2d_p2uv_C_sbp21")
        call create_p2uv_sbp_interpolator(interpolator2d, sbp_interp_c2i_21_t(), &
                                         domain,halo_width=1,staggering="C")
    case("interp2d_p2uv_C_sbp42")
        call create_p2uv_sbp_interpolator(interpolator2d, sbp_interp_c2i_42_t(), &
                                         domain,halo_width=3,staggering="C")
    case("interp2d_p2uv_C_sbp63")
        call create_p2uv_sbp_interpolator(interpolator2d, sbp_interp_c2i_63_t(), &
                                         domain,halo_width=5,staggering="C")
    case("interp2d_p2uv_Ch_sbp21")
        call create_p2uv_sbp_interpolator(interpolator2d, sbp_interp_i2c_21_t(), &
                                         domain,halo_width=1,staggering="Ch")
    case("interp2d_p2uv_Ch_sbp42")
        call create_p2uv_sbp_interpolator(interpolator2d, sbp_interp_i2c_42_t(), &
                                         domain,halo_width=3,staggering="Ch")
    case("interp2d_p2uv_Ch_sbp63")
        call create_p2uv_sbp_interpolator(interpolator2d, sbp_interp_i2c_63_t(), &
                                         domain,halo_width=5,staggering="Ch")
    case("interp2d_p2uv_identity")
        interpolator2d = p2uv_colocated_t()
    case default
        call parcomm_global%abort("create_scalar2vec_interpolator2d, unknown interpolator name: "//&
                                  interp2d_name)
    end select

end subroutine create_scalar2vec_interpolator2d

subroutine create_vec2vec_interpolator2d(interpolator2d, interp2d_name, domain)
    class(interpolator2d_vec2vec_t), allocatable,    intent(out) :: interpolator2d
    character(len=*),                                intent(in)  :: interp2d_name
    type(domain_t),                                  intent(in)  :: domain

    select case(interp2d_name)
    case("interp2d_pvec2uv_C_sbp21")
        call create_pvec2uv_sbp_interpolator(interpolator2d, sbp_interp_c2i_21_t(),&
                                         domain, halo_width=1, is_z = .false.)
    case("interp2d_pvec2uv_C_sbp42")
        call create_pvec2uv_sbp_interpolator(interpolator2d, sbp_interp_c2i_42_t(),&
                                         domain, halo_width=3, is_z = .false.)
    case("interp2d_pvec2uv_C_sbp63")
        call create_pvec2uv_sbp_interpolator(interpolator2d, sbp_interp_c2i_63_t(),&
                                         domain, halo_width=5, is_z = .false.)
    case("interp2d_pvec2uv_C_sbp21_z")
        call create_pvec2uv_sbp_interpolator(interpolator2d, sbp_interp_c2i_21_t(),&
                                         domain, halo_width=1, is_z = .true.)
    case("interp2d_pvec2uv_C_sbp42_z")
        call create_pvec2uv_sbp_interpolator(interpolator2d, sbp_interp_c2i_42_t(),&
                                         domain, halo_width=3, is_z = .true.)
    case("interp2d_pvec2uv_C_sbp63_z")
        call create_pvec2uv_sbp_interpolator(interpolator2d, sbp_interp_c2i_63_t(),&
                                         domain, halo_width=5, is_z = .true.)
    case("interp2d_uv2pvec_C_sbp21")
        call create_uv2pvec_sbp_interpolator(interpolator2d, sbp_interp_i2c_21_t(),&
                                             domain, halo_width=1,is_z=.false.)
    case("interp2d_uv2pvec_C_sbp42")
        call create_uv2pvec_sbp_interpolator(interpolator2d, sbp_interp_i2c_42_t(),&
                                             domain, halo_width=3,is_z=.false.)
    case("interp2d_uv2pvec_C_sbp63")
        call create_uv2pvec_sbp_interpolator(interpolator2d, sbp_interp_i2c_63_t(),&
                                             domain, halo_width=5,is_z=.false.)
    case("interp2d_uv2pvec_C_sbp21_z")
        call create_uv2pvec_sbp_interpolator(interpolator2d, sbp_interp_i2c_21_t(),&
                                             domain, halo_width=1,is_z=.true.)
    case("interp2d_uv2pvec_C_sbp42_z")
        call create_uv2pvec_sbp_interpolator(interpolator2d, sbp_interp_i2c_42_t(),&
                                             domain, halo_width=3,is_z=.true.)
    case("interp2d_uv2pvec_C_sbp63_z")
        call create_uv2pvec_sbp_interpolator(interpolator2d, sbp_interp_i2c_63_t(),&
                                             domain, halo_width=5,is_z=.true.)
    case("interp2d_pvec2uv_Ch_sbp21")
        call create_pvec2uv_Ch_sbp_interpolator(interpolator2d, sbp_interp_i2c_21_t(),&
                                             domain, halo_width=1)
    case("interp2d_pvec2uv_Ch_sbp42")
        call create_pvec2uv_Ch_sbp_interpolator(interpolator2d, sbp_interp_i2c_42_t(),&
                                             domain, halo_width=3)
    case("interp2d_pvec2uv_Ch_sbp63")
        call create_pvec2uv_Ch_sbp_interpolator(interpolator2d, sbp_interp_i2c_63_t(),&
                                             domain, halo_width=5)
    case("interp2d_uv2pvec_Ch_sbp21")
        call create_uv2pvec_sbp_Ch_interpolator(interpolator2d, sbp_interp_c2i_21_t(), &
                                                        domain, halo_width=1)
    case("interp2d_uv2pvec_Ch_sbp42")
        call create_uv2pvec_sbp_Ch_interpolator(interpolator2d, sbp_interp_c2i_42_t(), &
                                                        domain, halo_width=3)
    case("interp2d_uv2pvec_Ch_sbp63")
        call create_uv2pvec_sbp_Ch_interpolator(interpolator2d, sbp_interp_c2i_63_t(), &
                                                        domain, halo_width=6)
    case("interp2d_qvec2uv_Ch_sbp21")
        call create_qvec2uv_Ch_sbp_interpolator(interpolator2d, sbp_interp_c2i_21_t(),&
                                             domain, halo_width=1)
    case("interp2d_qvec2uv_Ch_sbp42")
        call create_qvec2uv_Ch_sbp_interpolator(interpolator2d, sbp_interp_c2i_42_t(),&
                                             domain, halo_width=3)
    case("interp2d_qvec2uv_Ch_sbp63")
        call create_qvec2uv_Ch_sbp_interpolator(interpolator2d, sbp_interp_c2i_63_t(),&
                                             domain, halo_width=5)
    case("interp2d_uv2qvec_C_sbp21")
        call create_uv2qvec_sbp_interpolator(interpolator2d, sbp_interp_c2i_21_t(),&
                                             domain, halo_width=1, is_Ch=.false.)
    case("interp2d_uv2qvec_C_sbp42")
        call create_uv2qvec_sbp_interpolator(interpolator2d, sbp_interp_c2i_42_t(),&
                                             domain, halo_width=3, is_Ch=.false.)
    case("interp2d_uv2qvec_C_sbp63")
        call create_uv2qvec_sbp_interpolator(interpolator2d, sbp_interp_c2i_63_t(),&
                                             domain, halo_width=5, is_Ch=.false.)
    case("interp2d_uv2qvec_Ch_sbp21")
        call create_uv2qvec_sbp_interpolator(interpolator2d, sbp_interp_i2c_21_t(),&
                                             domain, halo_width=1, is_Ch=.true.)
    case("interp2d_uv2qvec_Ch_sbp42")
        call create_uv2qvec_sbp_interpolator(interpolator2d, sbp_interp_i2c_42_t(),&
                                             domain, halo_width=3, is_Ch=.true.)
    case("interp2d_uv2qvec_Ch_sbp63")
        call create_uv2qvec_sbp_interpolator(interpolator2d, sbp_interp_i2c_63_t(),&
                                             domain, halo_width=5, is_Ch=.true.)
    case default
        call parcomm_global%abort("create_vec2vec_interpolator2d, unknown interpolator name: "//&
                                  interp2d_name)
    end select

end subroutine create_vec2vec_interpolator2d

subroutine create_p2uv_sbp_interpolator(interpolator_p2uv, sbp_interp_op, &
                                        domain, halo_width, staggering)

    class(interpolator2d_scalar2vec_t), allocatable, intent(out) :: interpolator_p2uv
    class(sbp_interp_t),                             intent(in)  :: sbp_interp_op
    type(domain_t),                                  intent(in)  :: domain
    integer(kind=4),                                 intent(in)  :: halo_width
    character(len=*),                                intent(in)  :: staggering

    type(interpolator_p2uv_sbp_C_t), allocatable :: interpolator_p2uv_sbp

    allocate(interpolator_p2uv_sbp)

    interpolator_p2uv_sbp%sbp_interp_p2uv = sbp_interp_op

    select case (staggering)
    case("C")
        interpolator_p2uv_sbp%exchange = &
            create_o_points_halo_exchange(domain%partition, domain%parcomm, &
                                      domain%topology, halo_width, 'cross')
    case("Ch")
        interpolator_p2uv_sbp%exchange = &
        create_xy_points_halo_exchange(domain%partition, domain%parcomm, &
                                  domain%topology, halo_width, 'cross')
    case default
        call domain%parcomm%abort("create_p2uv_sbp_interpolator, unsupported staggering: "//staggering//", currently C and Ch are vailable")
    end select

    call move_alloc(interpolator_p2uv_sbp, interpolator_p2uv)
end subroutine create_p2uv_sbp_interpolator

subroutine create_pvec2uv_sbp_interpolator(interpolator_p2uv, sbp_interp_op, &
                                                        domain, halo_width, is_z)

    class(interpolator2d_vec2vec_t), allocatable, intent(out) :: interpolator_p2uv
    class(sbp_interp_t),                          intent(in)  :: sbp_interp_op
    type(domain_t),                               intent(in)  :: domain
    integer(kind=4),                              intent(in)  :: halo_width
    logical,                                      intent(in)  :: is_z

    type(interpolator_pvec2uv_sbp_C_t), allocatable :: interpolator_p2uv_sbp

    allocate(interpolator_p2uv_sbp)

    interpolator_p2uv_sbp%sbp_interp_p2uv = sbp_interp_op

    if(is_z) then
        interpolator_p2uv_sbp%exchange = &
          create_z_points_halo_exchange(domain%partition, domain%parcomm, &
                                      domain%topology, halo_width, 'full')
        call domain%get_mesh(interpolator_p2uv_sbp%mesh_u,"xz")
        call domain%get_mesh(interpolator_p2uv_sbp%mesh_v,"yz")
    else
        interpolator_p2uv_sbp%exchange = &
          create_o_points_halo_exchange(domain%partition, domain%parcomm, &
                                      domain%topology, halo_width, 'full')
        call domain%get_mesh(interpolator_p2uv_sbp%mesh_u,"x")
        call domain%get_mesh(interpolator_p2uv_sbp%mesh_v,"y")
    end if

    call move_alloc(interpolator_p2uv_sbp, interpolator_p2uv)
end subroutine create_pvec2uv_sbp_interpolator

subroutine create_uv2pvec_sbp_interpolator(interpolator_uv2p, sbp_interp_op, &
                                                        domain, halo_width,is_z)

    class(interpolator2d_vec2vec_t), allocatable, intent(out) :: interpolator_uv2p
    class(sbp_interp_t),                          intent(in)  :: sbp_interp_op
    type(domain_t),                               intent(in)  :: domain
    integer(kind=4),                              intent(in)  :: halo_width
    logical,                                      intent(in)  :: is_z

    type(interpolator2d_uv2p_sbp_C_t), allocatable :: interpolator_uv2p_sbp

    allocate(interpolator_uv2p_sbp)

    interpolator_uv2p_sbp%sbp_interp_uv2p = sbp_interp_op

    if(is_z) then
        interpolator_uv2p_sbp%exchange = &
            create_symmetric_halo_vec_exchange_Cz(domain%partition, domain%parcomm, &
                                              domain%topology, halo_width, 'full')
        call domain%get_mesh(interpolator_uv2p_sbp%mesh_p,"z")
    else
        interpolator_uv2p_sbp%exchange = &
            create_symmetric_halo_vec_exchange_C(domain%partition, domain%parcomm, &
                                              domain%topology, halo_width, 'full')
        call domain%get_mesh(interpolator_uv2p_sbp%mesh_p,"o")
    end if

    call move_alloc(interpolator_uv2p_sbp, interpolator_uv2p)
end subroutine create_uv2pvec_sbp_interpolator

subroutine create_pvec2uv_Ch_sbp_interpolator(interpolator_p2uv, sbp_interp_op, &
                                                             domain, halo_width)

    class(interpolator2d_vec2vec_t), allocatable, intent(out) :: interpolator_p2uv
    class(sbp_interp_t),                          intent(in)  :: sbp_interp_op
    type(domain_t),                               intent(in)  :: domain
    integer(kind=4),                              intent(in)  :: halo_width

    type(interpolator_pvec2uv_sbp_Ch_t), allocatable :: interpolator_p2uv_sbp

    allocate(interpolator_p2uv_sbp)

    interpolator_p2uv_sbp%sbp_interp_p2uv = sbp_interp_op

    interpolator_p2uv_sbp%exchange = &
          create_xy_points_halo_exchange(domain%partition, domain%parcomm, &
                                       domain%topology, halo_width, 'full')

    call move_alloc(interpolator_p2uv_sbp, interpolator_p2uv)
end subroutine create_pvec2uv_Ch_sbp_interpolator

subroutine create_uv2pvec_sbp_Ch_interpolator(interpolator_uv2p, sbp_interp_op, &
                                                        domain, halo_width)

    class(interpolator2d_vec2vec_t), allocatable, intent(out) :: interpolator_uv2p
    class(sbp_interp_t),                          intent(in)  :: sbp_interp_op
    type(domain_t),                               intent(in)  :: domain
    integer(kind=4),                              intent(in)  :: halo_width

    type(interpolator2d_uv2p_sbp_C_t), allocatable :: interpolator_uv2p_sbp

    allocate(interpolator_uv2p_sbp)

    interpolator_uv2p_sbp%sbp_interp_uv2p = sbp_interp_op

    interpolator_uv2p_sbp%exchange = &
        create_symmetric_halo_vec_exchange_Ch(domain%partition, domain%parcomm, &
                                          domain%topology, halo_width, 'full')
    call domain%get_mesh(interpolator_uv2p_sbp%mesh_p,"xy")

    call move_alloc(interpolator_uv2p_sbp, interpolator_uv2p)
end subroutine create_uv2pvec_sbp_Ch_interpolator

subroutine create_qvec2uv_Ch_sbp_interpolator(interpolator_q2uv, sbp_interp_op, &
                                                             domain, halo_width)

    class(interpolator2d_vec2vec_t), allocatable, intent(out) :: interpolator_q2uv
    class(sbp_interp_t),                          intent(in)  :: sbp_interp_op
    type(domain_t),                               intent(in)  :: domain
    integer(kind=4),                              intent(in)  :: halo_width

    type(interpolator_q2uv_sbp_Ch_t), allocatable :: interpolator_q2uv_sbp

    allocate(interpolator_q2uv_sbp)

    interpolator_q2uv_sbp%sbp_interp_q2uv = sbp_interp_op

    interpolator_q2uv_sbp%exchange = &
          create_o_points_halo_exchange(domain%partition, domain%parcomm, &
                                       domain%topology, halo_width, 'full')

    call move_alloc(interpolator_q2uv_sbp, interpolator_q2uv)
end subroutine create_qvec2uv_Ch_sbp_interpolator

subroutine create_uv2qvec_sbp_interpolator(interpolator_uv2q, sbp_interp_op, &
                                                  domain, halo_width, is_Ch)

    class(interpolator2d_vec2vec_t), allocatable, intent(out) :: interpolator_uv2q
    class(sbp_interp_t),                          intent(in)  :: sbp_interp_op
    type(domain_t),                               intent(in)  :: domain
    integer(kind=4),                              intent(in)  :: halo_width
    logical,                                      intent(in)  :: is_Ch

    type(interpolator_uv2q_sbp_t), allocatable :: interpolator_uv2q_sbp

    allocate(interpolator_uv2q_sbp)

    interpolator_uv2q_sbp%is_Ch = is_Ch

    interpolator_uv2q_sbp%sbp_interp_uv2q = sbp_interp_op

    interpolator_uv2q_sbp%exchange = &
          create_symmetric_halo_vec_exchange_C(domain%partition, domain%parcomm, &
                                              domain%topology, halo_width, 'full')

    call move_alloc(interpolator_uv2q_sbp, interpolator_uv2q)
end subroutine create_uv2qvec_sbp_interpolator

end module interpolator2d_factory_mod
