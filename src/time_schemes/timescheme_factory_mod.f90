module timescheme_factory_mod

use stvec_mod,          only : stvec_t
use timescheme_mod,     only : timescheme_t
use parcomm_mod,        only : parcomm_global
use explicit_Eul1_mod,  only : explicit_Eul1_t
use implicit_Eul_mod,   only : implicit_Eul_t
use rk4_mod,            only : rk4_t
use ars343_mod,         only : ars343_t
use imex_M2_mod,        only : imex_M2_t, imex_M2_c
use crank_nikolson_mod, only : crank_nikolson_t
use lsrk_mod,           only : lsrk_t
use SISL_SETTLS_mod,    only : SISL_SETTLS_t
use SL_UKMO_mod,        only : SL_UKMO_t
use generic_config_mod, only : generic_config_t

implicit none

contains

subroutine create_timescheme(timescheme, v, timescheme_name, timescheme_config)
    class(timescheme_t), allocatable, intent(out) :: timescheme
    class(stvec_t),   intent(in) :: v !example of model state-vector
    character(len=*), intent(in) :: timescheme_name

    class(generic_config_t), intent(inout), optional :: timescheme_config

    select case(timescheme_name)
    case("explicit_Eul1")
        call create_explicit_Eul1_timescheme(timescheme, v)
    case("implicit_Eul")
        call create_implicit_Eul_timescheme(timescheme, v)
    case("crank_nikolson")
        call create_crank_nikolson_timescheme(timescheme, v)
    case("rk4")
        call create_rk4_timescheme(timescheme, v)
    case("ars343")
        call create_ars343_timescheme(timescheme,v)
    case("imex_M2a","imex_M2b","imex_M2c","imex_M2cn")
        call create_imex_M2_timescheme(timescheme,v,timescheme_name)
    case("lsrk")
        call create_lsrk_timescheme(timescheme,v)
    case("SISL_SETTLS")
        call create_SISL_SETTLS_timescheme(timescheme,v,timescheme_config)
    case("SL_UKMO")
        call create_SL_UKMO_timescheme(timescheme,v,timescheme_config)
    case default
        call parcomm_global%abort("unknown timescheme_name in create_timescheme: "// &
                                  timescheme_name)
    end select

end subroutine create_timescheme

subroutine create_explicit_Eul1_timescheme(timescheme, v)
    class(timescheme_t), allocatable, intent(out) :: timescheme
    class(stvec_t),                   intent(in)  :: v! model state-vector example

    type(explicit_Eul1_t), allocatable :: Eul1_timescheme

    allocate(Eul1_timescheme)
    call v%create_similar(Eul1_timescheme%tendency)
    call move_alloc(Eul1_timescheme, timescheme)

end subroutine create_explicit_Eul1_timescheme

subroutine create_implicit_Eul_timescheme(timescheme, v)
    class(timescheme_t), allocatable, intent(out) :: timescheme
    class(stvec_t),                   intent(in)  :: v! model state-vector example

    type(implicit_Eul_t), allocatable :: back_Eul_timescheme

    allocate(back_Eul_timescheme)
    call v%create_similar(back_Eul_timescheme%v_new)
    call move_alloc(back_Eul_timescheme, timescheme)

end subroutine create_implicit_Eul_timescheme

subroutine create_crank_nikolson_timescheme(timescheme, v)
    class(timescheme_t), allocatable, intent(out) :: timescheme
    class(stvec_t),                   intent(in)  :: v! model state-vector example

    type(crank_nikolson_t), allocatable :: crank_nikolson_timescheme

    allocate(crank_nikolson_timescheme)
    call v%create_similar(crank_nikolson_timescheme%v_new)
    call move_alloc(crank_nikolson_timescheme, timescheme)

end subroutine create_crank_nikolson_timescheme

subroutine create_rk4_timescheme(timescheme, v)

    class(timescheme_t), allocatable, intent(out) :: timescheme
    class(stvec_t),                   intent(in)  :: v !example of model state vector

    type(rk4_t), allocatable :: rk4

    allocate(rk4)

    !preallocate additional state vectors
    call v%create_similar(rk4%k1)
    call v%create_similar(rk4%k2)
    call v%create_similar(rk4%k3)
    call v%create_similar(rk4%k4)
    call v%create_similar(rk4%y)

    call move_alloc(rk4, timescheme)
end subroutine create_rk4_timescheme

subroutine create_ars343_timescheme(timescheme, v)
    class(timescheme_t), allocatable, intent(out) :: timescheme
    class(stvec_t),                   intent(in)  :: v! model state-vector example

    type(ars343_t), allocatable :: ars343

    allocate(ars343)
    call v%create_similar(ars343%y1)
    call v%create_similar(ars343%y2)
    call v%create_similar(ars343%y3)
    call v%create_similar(ars343%y4)
    call v%create_similar(ars343%q2)
    call v%create_similar(ars343%q3)
    call v%create_similar(ars343%q4)
    call v%create_similar(ars343%r)
    call v%create_similar(ars343%s)

    call move_alloc(ars343, timescheme)

end subroutine create_ars343_timescheme

subroutine create_imex_M2_timescheme(timescheme, v, timescheme_name)

    class(timescheme_t), allocatable, intent(out) :: timescheme
    class(stvec_t),                   intent(in)  :: v !example of model state vector
    character(len=*),                 intent(in)  :: timescheme_name

    type(imex_M2_t), allocatable :: imex_M2

    allocate(imex_M2)

    !preallocate additional state vectors
    call v%create_similar(imex_M2%ye)
    call v%create_similar(imex_M2%yi)
    call v%create_similar(imex_M2%s)
    call v%create_similar(imex_M2%r)
    call v%create_similar(imex_M2%res)

    imex_M2%c(1:5) = imex_M2_c(1:5)

    select case(timescheme_name)
    case("imex_M2a")
        imex_M2%d(1:6) = [3._8/11._8, 0._8, 3._8/11._8, 0._8, 0._8, 5._8/11._8]
    case("imex_M2b")
        imex_M2%d(1:6) = [0._8, 0._8, 3._8/5._8, 0._8, 0._8, 2._8/5._8]
    case("imex_M2c")
        imex_M2%d(1:6) = [2._8/7._8, 2._8/7._8, 0._8, 0._8, 0._8, 4._8/11._8]
    case("imex_M2cn")
        imex_M2%d(1:6) = [0.5_8, 0._8, 0._8, 0._8, 0._8, 0.5_8]
    case default
        call parcomm_global%abort("unknown imex_M2 scheme: "//timescheme_name//", stopping")
    end select

    call move_alloc(imex_M2, timescheme)
end subroutine create_imex_M2_timescheme

subroutine create_lsrk_timescheme(timescheme, v)
    class(timescheme_t), allocatable, intent(out) :: timescheme
    class(stvec_t),                   intent(in) :: v

    type(lsrk_t), allocatable :: lsrk_timescheme

    allocate(lsrk_timescheme)
    !preallocate additional state vectors
    call v%create_similar(lsrk_timescheme%k1)
    call v%create_similar(lsrk_timescheme%k2)
    call v%create_similar(lsrk_timescheme%temp)

    call move_alloc(lsrk_timescheme, timescheme)

end subroutine create_lsrk_timescheme

subroutine create_SISL_SETTLS_timescheme(timescheme, v, timescheme_config)
    class(timescheme_t), allocatable, intent(out) :: timescheme
    class(stvec_t),                   intent(in) :: v
    class(generic_config_t),          intent(inout), optional :: timescheme_config

    type(SISL_SETTLS_t), allocatable :: sisl_timescheme

    allocate(sisl_timescheme)
    !preallocate additional state vectors
    call v%create_similar(sisl_timescheme%v_dp)
    call v%create_similar(sisl_timescheme%rhs)
    call v%create_similar(sisl_timescheme%explicit_tend)
    call v%create_similar(sisl_timescheme%explicit_tend_old)
    call v%create_similar(sisl_timescheme%wind_arr)
    call v%create_similar(sisl_timescheme%wind_dp)

    call move_alloc(sisl_timescheme, timescheme)

end subroutine create_SISL_SETTLS_timescheme

subroutine create_SL_UKMO_timescheme(timescheme, v, timescheme_config)
    class(timescheme_t), allocatable, intent(out) :: timescheme
    class(stvec_t),                   intent(in) :: v
    class(generic_config_t),          intent(inout), optional :: timescheme_config

    type(SL_UKMO_t), allocatable :: sisl_timescheme

    allocate(sisl_timescheme)

    if(present(timescheme_config)) then
        call timescheme_config%get(sisl_timescheme%num_iter,"num_iter", default=sisl_timescheme%num_iter)
        call timescheme_config%get(sisl_timescheme%epsilon, "epsilon",  default=sisl_timescheme%epsilon)
    !else !if config is not present SL_UKMO_t defaults are used for num_iter and epsilon
    end if

    !preallocate additional state vectors
    call v%create_similar(sisl_timescheme%v_dp)
    call v%create_similar(sisl_timescheme%rhs)
    call v%create_similar(sisl_timescheme%wind_arr)
    call v%create_similar(sisl_timescheme%wind_dp)

    call move_alloc(sisl_timescheme, timescheme)

end subroutine create_SL_UKMO_timescheme

end module timescheme_factory_mod
