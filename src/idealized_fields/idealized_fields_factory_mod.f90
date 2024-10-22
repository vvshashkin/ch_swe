module idealized_fields_factory_mod

use idealized_field_mod, only : idealized_field_t
use generic_config_mod,  only : generic_config_t
use parcomm_mod,         only : parcomm_global
use domain_mod,          only : domain_t
use string_mod,          only : string_t

implicit none

interface create_idealized_field
    module procedure :: create_testcase_idealized_field
end interface

contains

function create_testcase_idealized_field(field_name, component_type, &
                                         testcase_config, domain) result(field)
    character(len=*),        intent(in)    :: field_name
    character(len=*),        intent(in)    :: component_type
    class(generic_config_t), intent(inout) :: testcase_config
    type(domain_t),          intent(in)    :: domain

    class(idealized_field_t), allocatable  :: field

    type(string_t) :: testcase_name

    call testcase_config%get(testcase_name,"name")

    select case(testcase_name%str)
    case("solid_rotation")
        field = create_solid_rotation_field(field_name,component_type,&
                                                        testcase_config,domain)
    case("baroclinic_instability")
        field = create_baroclinic_instability_field(field_name,component_type,&
                                                        testcase_config,domain)
    case("GW")
        field = create_GW_field(field_name,component_type,testcase_config,domain)
    case default
        call parcomm_global%abort(__FILE__//": unknown testcase name "//testcase_name%str)
    end select
    return
end function

function create_solid_rotation_field(field_name,component_type,config,domain) result(field)

    use solid_rotation_fields_mod, only : solid_rotation_vector_field_t,     &
                                          solid_rotation_isotermic_theta_t,  &
                                          solid_rotation_isotermic_PExner_t, &
                                          solid_rotation_isotermic_rho_t,    &
                                          solid_rotation_Nbw_PExner_t,       &
                                          solid_rotation_Nbw_theta_t
    use const_mod,                 only : Earth_radii, Earth_grav

    character(len=*),        intent(in)    :: field_name
    character(len=*),        intent(in)    :: component_type
    class(generic_config_t), intent(inout) :: config
    type(domain_t),          intent(in)    :: domain

    class(idealized_field_t), allocatable  :: field

    type(string_t) :: mode
    real(kind=8)   :: omega, alpha, u0, T0, Nbw, sphere_rad, grav, peq, p0
    logical        :: found1, found2

    call config%get(mode,"mode",default=string_t("isotermic"))
    call config%get(T0,"T0",default=300.0_8)
    call config%get(omega,"omega",default=domain%mesh_p%omega)
    call config%get(alpha,"alpha",default=0.0_8)
    call config%get(u0,"u0",default=20.0_8)
    call config%get(sphere_rad,"R",default=domain%mesh_p%scale)
    call config%get(grav,"g",Earth_grav)
    call config%get(Nbw,"Nbw",default=0.01_8,found=found1)
    call config%get(peq,"peq",default=1e5_8,found=found2)
    if(mode%str=="isotermic" .and. (found1 .or. found2)) &
        call parcomm_global%abort(__FILE__//": solid_rotation_tescase flow config error: mode is isotermic, but Brunt-Vaisala frequency (Nbw) and/or pressure at equator (peq) are set")
    call config%get(p0,"p0",default=93000.0_8,found=found1)
    if(mode%str=="Nbw_const" .and. found1) &
        call parcomm_global%abort(__FILE__//": solid_rotation_tescase flow config error: mode is Nbw_const, but pole-point pressure (p0) is set")


    select case(field_name)
    case("u")
        field = solid_rotation_vector_field_t(component_direction="x",component_type=component_type,u0=u0,alpha=alpha)
        return
    case("v")
        field = solid_rotation_vector_field_t(component_direction="y",component_type=component_type,u0=u0,alpha=alpha)
        return
    case("w")
        field = solid_rotation_vector_field_t(component_direction="z",component_type=component_type,u0=u0,alpha=alpha)
        return
    end select

    select case(mode%str)
    case ("isotermic")
        select case(field_name)
        case("theta")
            field = solid_rotation_isotermic_theta_t(u0=u0,alpha=alpha,T0=T0,&
                                                     omega=omega,sphere_rad=sphere_rad,grav=grav,p0=p0)
        case("P","PExner")
            field = solid_rotation_isotermic_PExner_t(u0=u0,alpha=alpha,T0=T0,&
                                                      omega=omega,sphere_rad=sphere_rad,grav=grav,p0=p0)
        case("rho")
            field = solid_rotation_isotermic_rho_t(u0=u0,alpha=alpha,T0=T0,&
                                                   omega=omega,sphere_rad=sphere_rad,grav=grav,p0=p0)
        case default
            call parcomm_global%abort(__FILE__//": solid rotation testcase factory unknown field: "//field_name)
        end select

    case("Nbw_const")
        select case(field_name)
        case("theta")
            field = solid_rotation_Nbw_theta_t(u0=u0,alpha=alpha,T0=T0,omega=omega,Nbw=Nbw, &
                                               sphere_rad=sphere_rad,grav=grav,peq=peq)
        case("P","PExner")
            field = solid_rotation_Nbw_PExner_t(u0=u0,alpha=alpha,T0=T0,omega=omega,Nbw=Nbw, &
                                                sphere_rad=sphere_rad,grav=grav,peq=peq)
        case default
            call parcomm_global%abort(__FILE__//": solid rotation testcase factory unknown field: "//field_name)
        end select
    case default
        call parcomm_global%abort(__FILE__//": solid rotation testcase factory unknown mode: "//mode%str)
    end select
end function

function create_baroclinic_instability_field(field_name,component_type,config,domain) result(field)

    use baroclinic_instability_JW2006_mod, only : baroclinic_instability_wind_t,   &
                                                  baroclinic_instability_PExner_t, &
                                                  baroclinic_instability_theta_t,  &
                                                  baroclinic_instability_rho_t

    character(len=*),        intent(in)    :: field_name
    character(len=*),        intent(in)    :: component_type
    class(generic_config_t), intent(inout) :: config
    type(domain_t),          intent(in)    :: domain

    class(idealized_field_t), allocatable  :: field


    real(kind=8) :: u_pert

    call config%get(u_pert,"u_pert",default=1.0_8)

    select case(field_name)
    case("u")
        field = baroclinic_instability_wind_t(component_direction="x",component_type=component_type,u_pert=u_pert)
        return
    case("v")
        field = baroclinic_instability_wind_t(component_direction="y",component_type=component_type,u_pert=u_pert)
        return
    case("w")
        field = baroclinic_instability_wind_t(component_direction="z",component_type=component_type,u_pert=u_pert)
        return
    case("P","PExner")
        field = baroclinic_instability_PExner_t()
    case("theta")
        field = baroclinic_instability_theta_t()
    case("rho")
        field = baroclinic_instability_rho_t()
    case default
        call parcomm_global%abort(__FILE__//": baroclinic instability tescase factory error: unknown field name")
    end select

end function

function create_GW_field(field_name,component_type,config,domain) result(field)

    use solid_rotation_fields_mod,  only : solid_rotation_vector_field_t
    use GW_testcase_fields_mod,     only : GW_theta_t, GW_theta_linear_t, &
                                           GW_PExner_t, GW_rho_t
    use basic_idealized_fields_mod, only : zero_idealized_field_t
    use const_mod,                  only : Earth_radii, Earth_grav

    character(len=*),        intent(in)    :: field_name
    character(len=*),        intent(in)    :: component_type
    class(generic_config_t), intent(inout) :: config
    type(domain_t),          intent(in)    :: domain

    class(idealized_field_t), allocatable  :: field

    type(string_t) :: mode
    real(kind=8)   :: u0, Nbw, sphere_rad, d, theta_pert, Lz

    call config%get(u0,"u0",default=0.0_8)
    call config%get(Nbw,"Nbw",default=0.01_8)
    call config%get(sphere_rad,"R",default=domain%mesh_p%scale)
    call config%get(d,"d",default=5e3_8)
    call config%get(theta_pert,"theta_pert",1.0_8)
    call config%get(Lz,"Lz",1e4_8)
    call config%get(mode,"mode",default=string_t("non-linear"))
    if(mode%str == "linear" .and. u0 /= 0.0_8) &
        call parcomm_global%abort(__FILE__//" GW testcase config error: mode is linear, but u0 is set non zero")

    select case(field_name)
    case("u")
        field = solid_rotation_vector_field_t(component_direction="x",component_type=component_type,u0=u0,alpha=0.0_8)
        return
    case("v")
        field = solid_rotation_vector_field_t(component_direction="y",component_type=component_type,u0=u0,alpha=0.0_8)
        return
    case("w")
        field = solid_rotation_vector_field_t(component_direction="z",component_type=component_type,u0=u0,alpha=0.0_8)
        return
    end select

    select case(mode%str)
    case("linear")
        select case(field_name)
        case("theta")
            field = GW_theta_linear_t(u0=0.0_8,d=d,sphere_rad=sphere_rad,Nbw=Nbw,Lz=Lz,amp=theta_pert)
        case("PExner","P")
            field = zero_idealized_field_t()
        case default
            call parcomm_global%abort(__FILE__//" GW testcase factory, unknown field: "//field_name)
        end select
    case("non-linear","non_linear","non linear")
        select case(field_name)
        case("theta")
            field = GW_theta_t(u0=u0,d=d,sphere_rad=sphere_rad,Nbw=Nbw,Lz=Lz,amp=theta_pert)
        case("PExner","P")
            field = GW_PExner_t(u0=u0,sphere_rad=sphere_rad,Nbw=Nbw)
        case("rho")
            field = GW_rho_t(u0=u0,d=d,sphere_rad=sphere_rad,Nbw=Nbw,Lz=Lz,amp=theta_pert)
        case default
            call parcomm_global%abort(__FILE__//" GW testcase factory, unknown field: "//field_name)
        end select
    case default
        call parcomm_global%abort(__FILE__//" GW_testcase, wrong mode in config:"//mode%str//", use linear or non-linear")
    end select

end function

function create_vertical_profile_idealized_field(field_name,config,domain, vertical_derivative) result(field)

    use basic_idealized_fields_mod, only : zero_idealized_field_t
    use vertical_profiles_mod,      only : const_Nbw_profile_theta_t, &
                                           const_Nbw_profile_dtheta_dz_t, &
                                           const_Nbw_profile_P_t, &
                                           const_Nbw_profile_dP_dz_t, &
                                           isotermal_profile_P_t, &
                                           isotermal_profile_dP_dz_t, &
                                           isotermal_profile_theta_t, &
                                           isotermal_profile_dtheta_dz_t, &
                                           BW_profile_theta_t, &
                                           BW_profile_dtheta_dz_t, &
                                           BW_profile_P_t, BW_profile_dP_dz_t

    use const_mod,                  only : Earth_radii, Earth_grav

    character(len=*),        intent(in)    :: field_name
    class(generic_config_t), intent(inout) :: config
    type(domain_t),          intent(in)    :: domain
    logical,                 intent(in)    :: vertical_derivative

    class(idealized_field_t), allocatable  :: field

    character(len=:), allocatable :: profile_name
    real(kind=8) :: Nbw, T0, P0, grav

    call config%get(profile_name,"name",default="zero")

    select case(profile_name)
    case("zero")
        field = zero_idealized_field_t()
    case("const_Nbw")

        call config%get(Nbw,"Nbw",default=0.01_8)
        call config%get(T0,"T0",default=300.0_8)
        call config%get(grav,"grav",default=Earth_grav)

        select case (field_name)
        case("theta")
            if(vertical_derivative) then
                field = const_Nbw_profile_dtheta_dz_t(Nbw=Nbw,T0=T0,grav=grav)
            else
                field = const_Nbw_profile_theta_t(Nbw=Nbw,T0=T0,grav=grav)
            end if
        case("P")
            if(vertical_derivative) then
                field = const_Nbw_profile_dP_dz_t(Nbw=Nbw,T0=T0,grav=grav)
            else
                field = const_Nbw_profile_P_t(Nbw=Nbw,T0=T0,grav=grav)
            end if
        case default
        call parcomm_global%abort("create_vertical_profile_idealized_field, unknown const_Nbw field: "// field_name)
        end select

    case("isotermal")

        call config%get(T0,"T0",default=300.0_8)
        call config%get(P0,"p_surf",default=1.0e5_8)
        call config%get(grav,"grav",default=Earth_grav)

        select case (field_name)
        case("theta")
            if(vertical_derivative) then
                field = isotermal_profile_dtheta_dz_t(T0=T0, p_surf=P0, grav=grav)
            else
                field = isotermal_profile_theta_t(T0=T0, p_surf=P0, grav=grav)
            end if
        case("P")
            if(vertical_derivative) then
                field = isotermal_profile_dP_dz_t(T0=T0, p_surf=P0, grav=grav)
            else
                field = isotermal_profile_P_t(T0=T0, p_surf=P0, grav=grav)
            end if
        case default
            call parcomm_global%abort("create_vertical_profile_idealized_field, unknown const_Nbw field: "// field_name)
        end select
    
    case("BW_profile")

        select case (field_name)
        case("theta")
            if(vertical_derivative) then
                field = BW_profile_dtheta_dz_t()
            else
                field = BW_profile_theta_t()
            end if
        case("P")
            if(vertical_derivative) then
                field = BW_profile_dP_dz_t()
            else
                field = BW_profile_P_t()
            end if
        case default
            call parcomm_global%abort("create_vertical_profile_idealized_field, unknown const_Nbw field: "// field_name)
        end select


    case default
        call parcomm_global%abort("create_vertical_profile_idealized_field, unknown profile name: "// profile_name)
    end select

end function

end module
