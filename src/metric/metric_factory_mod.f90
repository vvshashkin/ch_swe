module metric_factory_mod

use metric_mod,         only : metric_t
use topology_mod,       only : topology_t
use parcomm_mod,        only : parcomm_global
use config_metric_mod,  only : config_metric_t
use generic_config_mod, only : generic_config_t
use string_mod,         only : string_t

implicit none

character(len=*), private, parameter :: METRIC_2d_DEFAULT = "equiangular_cubed_sphere"
character(len=*), private, parameter :: VERTICAL_TRANSFORM_DEFAULT = "vertical_transform_default"
real(kind=8), private, parameter :: SCALE_DEFAULT = 1.0_8, VERTICAL_SCALE_DEFAULT = 1.0_8, OMEGA_DEFAULT = 1.0_8
real(kind=8), private, parameter :: ROTATION_AXIS_DEFAULT(3) = [0.0_8, 0.0_8, 1.0_8]
real(kind=8), private, parameter :: ROTATION_MATRIX_DEFAULT(9) = &
                                             [1.0_8,0.0_8,0.0_8, &
                                              0.0_8,1.0_8,0.0_8, &
                                              0.0_8,0.0_8,1.0_8 ]

interface create_metric
    module procedure :: create_metric_by_param
    module procedure :: create_metric_by_config
    module procedure :: create_metric_by_generic_config
end interface

contains

subroutine create_metric_by_generic_config(metric, topology, config)

    class(metric_t), allocatable, intent(out)    :: metric
    class(topology_t),            intent(in)     :: topology
    class(generic_config_t),      intent(inout)  :: config

    character(len=:), allocatable :: metric_type
    type(string_t), allocatable :: unexpected_vars(:), expected_vars(:)
    logical         :: check
    integer(kind=4) :: i

    expected_vars = [string_t("metric_type"),&
                     string_t("vertical_transform_name"), &
                     string_t("h_top"), &
                     string_t("planet_radius"), &
                     string_t("omega"), &
                     string_t("metric_2d_type"), &
                     string_t("rotation_matrix"),&
                     string_t("rotation_axis"), &
                     string_t("coriolis_constant")]

    check = config%check_no_unexpected(expected_vars,&
                        unexpected_varnames=unexpected_vars)
    if(.not. check) then
        do i = 1, size(unexpected_vars,1)
            call parcomm_global%print("metric factory WARNING: config variable "//unexpected_vars(i)%str//" presumably is not expected")
        end do
    end if

    call config%get(metric_type,"metric_type","metric type not set in config")

    select case(metric_type)
    case("shallow_atmosphere_metric")
        call create_shallow_atmosphere_metric(metric, topology, config)
    case default
        call parcomm_global%abort("Unknown metric_type var " // metric_type // &
                                                    " in metric_factory_mod")
    end select
end subroutine create_metric_by_generic_config

subroutine create_metric_by_config(metric, topology, metric_type, config)

    class(metric_t), allocatable, intent(out) :: metric
    class(topology_t),            intent(in)  :: topology
    character(len=*),             intent(in)  :: metric_type
    type(config_metric_t),        intent(in)  :: config

    select case(metric_type)
    case("shallow_atmosphere_metric")
        call create_shallow_atmosphere_metric_old(metric, topology, metric_type, config)
    case default
        call parcomm_global%abort("Unknown metric_type var " // metric_type // &
                                                    " in metric_factory_mod")
    end select
end subroutine create_metric_by_config

subroutine create_metric_by_param(metric, topology, metric_type)

    class(topology_t),            intent(in)  :: topology
    character(len=*),             intent(in)  :: metric_type
    class(metric_t), allocatable, intent(out) :: metric

    type(config_metric_t) :: config

    call config%set_defaults()
    call create_metric_by_config(metric,topology,metric_type,config)

end subroutine create_metric_by_param

subroutine create_shallow_atmosphere_metric(metric, topology, config)

    use shallow_atm_metric_mod,         only : shallow_atm_metric_t
    use vertical_transform_factory_mod, only : create_vertical_transform

    class(metric_t), allocatable, intent(out)    :: metric
    class(topology_t),            intent(in)     :: topology
    class(generic_config_t),      intent(inout)  :: config

    type(shallow_atm_metric_t), allocatable :: metric_shallow_atm
    character(len=:), allocatable :: metric_2d_type, vertical_transform_name
    real(kind=8), allocatable :: rotation_matrix(:)

    allocate(metric_shallow_atm)

    call config%get(metric_2d_type,"metric_2d_type",default=METRIC_2D_DEFAULT)
    call create_metric_2d(metric_shallow_atm%metric_2d, topology, metric_2d_type, config)

    metric_shallow_atm%scale = metric_shallow_atm%metric_2d%scale
    call config%get(metric_shallow_atm%vertical_scale,"h_top",default=VERTICAL_SCALE_DEFAULT)
    metric_shallow_atm%omega = metric_shallow_atm%metric_2d%omega
    metric_shallow_atm%rotation_axis = &
            metric_shallow_atm%metric_2d%rotation_axis
    metric_shallow_atm%rotation_matrix = &
            metric_shallow_atm%metric_2d%rotation_matrix
    metric_shallow_atm%alpha0 = metric_shallow_atm%metric_2d%alpha0
    metric_shallow_atm%alpha1 = metric_shallow_atm%metric_2d%alpha1
    metric_shallow_atm%beta0  = metric_shallow_atm%metric_2d%beta0
    metric_shallow_atm%beta1  = metric_shallow_atm%metric_2d%beta1

    call config%get(vertical_transform_name,"vertical_transform_name", &
                    default=VERTICAL_TRANSFORM_DEFAULT)
    metric_shallow_atm%vertical_transform = &
            create_vertical_transform(vertical_transform_name)

    call config%get(metric_shallow_atm%coriolis_const,"coriolis_constant",default=.false.)

    call move_alloc(metric_shallow_atm, metric)

end subroutine create_shallow_atmosphere_metric

subroutine create_shallow_atmosphere_metric_old(metric, topology, metric_type, config)

    use shallow_atm_metric_mod,         only : shallow_atm_metric_t
    use vertical_transform_factory_mod, only : create_vertical_transform

    class(metric_t), allocatable, intent(out) :: metric
    class(topology_t),            intent(in)  :: topology
    character(len=*),             intent(in)  :: metric_type
    type(config_metric_t),        intent(in)  :: config

    type(shallow_atm_metric_t), allocatable :: metric_shallow_atm

    allocate(metric_shallow_atm)

    call create_metric_2d_old(metric_shallow_atm%metric_2d, topology, config%metric_2d_type, config)

    metric_shallow_atm%scale           = config%scale
    metric_shallow_atm%vertical_scale  = config%vertical_scale
    metric_shallow_atm%omega           = config%omega
    metric_shallow_atm%rotation_axis   = config%rotation_axis
    metric_shallow_atm%rotation_matrix = config%rotation_matrix
    metric_shallow_atm%alpha0          = metric_shallow_atm%metric_2d%alpha0
    metric_shallow_atm%alpha1          = metric_shallow_atm%metric_2d%alpha1
    metric_shallow_atm%beta0           = metric_shallow_atm%metric_2d%beta0
    metric_shallow_atm%beta1           = metric_shallow_atm%metric_2d%beta1

    metric_shallow_atm%vertical_transform = &
                    create_vertical_transform(config%vertical_transform_name)

    call move_alloc(metric_shallow_atm, metric)

end subroutine create_shallow_atmosphere_metric_old

subroutine create_metric_2d(metric, topology, metric_type, config)

    use metric_2d_mod,          only : metric_2d_t

    class(metric_2d_t),  allocatable, intent(out)   :: metric
    class(topology_t),                intent(in)    :: topology
    character(len=*),                 intent(in)    :: metric_type
    class(generic_config_t),          intent(inout) :: config

    real(kind=8) :: scale, omega
    real(kind=8), allocatable :: rotation_axis(:), buff(:)
    real(kind=8) :: rotation_matrix(3,3)

    call config%get(scale,"planet_radius",default=SCALE_DEFAULT)
    call config%get(omega,"omega",default=OMEGA_DEFAULT)
    call config%get(rotation_axis,"rotation_axis", &
                    default=ROTATION_AXIS_DEFAULT)
    call config%get(buff,"rotation_matrix", &
                    default = ROTATION_MATRIX_DEFAULT)
    rotation_matrix = reshape(buff,[3,3])

    select case(metric_type)
    case("ecs","equiangular_cubed_sphere")
        call create_ecs_2d_metric(metric, topology, &
         scale, omega, rotation_matrix, rotation_axis)
    case default
        call parcomm_global%abort("Unknown metric_type var " // metric_type // &
                                                    " in metric_factory_mod")
    end select
    
end subroutine create_metric_2d

subroutine create_metric_2d_old(metric, topology, metric_type, config)

    use metric_2d_mod,          only : metric_2d_t

    class(metric_2d_t),  allocatable, intent(out) :: metric
    class(topology_t),                intent(in)  :: topology
    character(len=*),                 intent(in)  :: metric_type
    type(config_metric_t),            intent(in)  :: config

    select case(metric_type)
    case("ecs")
        call create_ecs_2d_metric(metric, topology, &
         config%scale, config%omega, config%rotation_matrix, config%rotation_axis)
    case default
        call parcomm_global%abort("Unknown metric_type var " // metric_type // &
                                                    " in metric_factory_mod")
    end select
    
end subroutine create_metric_2d_old

subroutine create_ecs_2d_metric(metric, topology, sphere_r, sphere_omega, rotation_matrix, rotation_axis)

    use cubed_sphere_topology_mod, only : cubed_sphere_topology_t
    use topology_mod,              only : topology_t
    use metric_2d_mod,             only : metric_2d_t
    use metric_2d_ecs_mod,         only : metric_2d_ecs_t
    use const_mod,                 only : pi

    class(topology_t),               intent(in)  :: topology
    class(metric_2d_t), allocatable, intent(out) :: metric
    real(kind=8),       optional,    intent(in)  :: sphere_r, sphere_omega
    real(kind=8),       optional,    intent(in)  :: rotation_matrix(3,3), rotation_axis(3)

    type(metric_2d_ecs_t), allocatable :: ecs_metric
    real(kind=8) :: local_rotation_matrix(3, 3), local_rotation_axis(3)
    real(kind=8) :: local_sphere_r, local_sphere_omega

    integer(kind=4) :: i_panel

    local_rotation_matrix = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])
    if (present(rotation_matrix)) local_rotation_matrix = rotation_matrix

    local_rotation_axis = [0, 0, 1]
    if (present(rotation_axis)) local_rotation_axis = rotation_axis

    local_sphere_omega = 1.0_8
    if (present(sphere_omega)) local_sphere_omega = sphere_omega

    local_sphere_r = 1.0_8
    if (present(sphere_r)) local_sphere_r = sphere_r

    allocate(ecs_metric)

    ecs_metric%alpha0 =-0.25*pi
    ecs_metric%beta0  =-0.25*pi
    ecs_metric%alpha1 = 0.25*pi
    ecs_metric%beta1  = 0.25*pi

    ecs_metric%scale = local_sphere_r
    ecs_metric%omega = local_sphere_omega

    ecs_metric%rotation_matrix = local_rotation_matrix
    ecs_metric%rotation_axis   = local_rotation_axis

    do i_panel = 1, 6
        ecs_metric%rotation_matrix_panel(1,1,i_panel) = &
            sum(topology%ex(1:3,i_panel)*local_rotation_matrix(1:3,1)) 
        ecs_metric%rotation_matrix_panel(2,1,i_panel) = &
            sum(topology%ey(1:3,i_panel)*local_rotation_matrix(1:3,1))
        ecs_metric%rotation_matrix_panel(3,1,i_panel) = &
            -sum(topology%n(1:3,i_panel)*local_rotation_matrix(1:3,1))

        ecs_metric%rotation_matrix_panel(1,2,i_panel) = &
            sum(topology%ex(1:3,i_panel)*local_rotation_matrix(1:3,2)) 
        ecs_metric%rotation_matrix_panel(2,2,i_panel) = &
            sum(topology%ey(1:3,i_panel)*local_rotation_matrix(1:3,2))
        ecs_metric%rotation_matrix_panel(3,2,i_panel) = &
            -sum(topology%n(1:3,i_panel)*local_rotation_matrix(1:3,2))

        ecs_metric%rotation_matrix_panel(1,3,i_panel) = &
            sum(topology%ex(1:3,i_panel)*local_rotation_matrix(1:3,3)) 
        ecs_metric%rotation_matrix_panel(2,3,i_panel) = &
            sum(topology%ey(1:3,i_panel)*local_rotation_matrix(1:3,3))
        ecs_metric%rotation_matrix_panel(3,3,i_panel) = &
            -sum(topology%n(1:3,i_panel)*local_rotation_matrix(1:3,3))
    end do

    select type(topology)
    class is (cubed_sphere_topology_t)
        ecs_metric%topology = topology
    class default
        call parcomm_global%abort("Wrong topology class in create_ecs_metric!")
    end select

    call move_alloc(ecs_metric, metric)

end subroutine create_ecs_2d_metric


end module metric_factory_mod
