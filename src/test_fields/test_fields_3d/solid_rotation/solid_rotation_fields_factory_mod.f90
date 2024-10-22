module solid_rotation_fields_factory_mod

use test_fields_3d_mod,            only : scalar_field3d_t, vector_field3d_t
use solid_rotation3d_therm_mod,    only : solid_rotation_theta_Nb_t, &
                                          solid_rotation_PExner_Nb_t,&
                                          solid_rotation_theta_isoT_t, &
                                          solid_rotation_PExner_isoT_t
use solid_rotation_wind_field_mod, only : solid_rotation_wind_field_t
use parcomm_mod,                   only : parcomm_global

implicit none

contains

subroutine create_solid_rotation_field_generators(background_type,&
                                                  u0,omega,sphere_rad,Nb,T0,grav,alpha,&
                                                  pref, peq, p0, &
                                                  theta_gen, Pexner_gen, wind_gen)
    character(len=*) :: background_type
    real(kind=8), intent(in) :: u0,omega,sphere_rad,grav,alpha
    real(kind=8), intent(in), optional :: Nb, T0, peq, pref, p0
    class(scalar_field3d_t), optional, allocatable, intent(out) :: theta_gen, Pexner_gen
    class(vector_field3d_t), optional, allocatable, intent(out) :: wind_gen

    real(kind=8) :: Nb_loc, T0_loc, p0_loc, pref_loc, peq_loc

    p0_loc = 930.0e2_8
    if(present(p0)) p0_loc = p0
    peq_loc = 1e5
    if(present(peq)) peq_loc = peq
    pref_loc = 1e5
    if(present(pref)) pref_loc = pref

    select case(background_type)
    case("Nb_const")
        Nb_loc = 0.01_8
        if(present(Nb)) Nb_loc = Nb

        if(present(theta_gen)) &
            theta_gen = solid_rotation_theta_Nb_t(u0=u0,omega=omega,sphere_rad=sphere_rad,&
                                                  Nb=Nb_loc,grav=grav,alpha=alpha,pref=pref_loc,peq=peq_loc)
        if(present(PExner_gen)) &
            PExner_gen = solid_rotation_PExner_Nb_t(u0=u0,omega=omega,sphere_rad=sphere_rad,&
                                                   Nb=Nb_loc,grav=grav,alpha=alpha,pref=pref_loc,peq=peq_loc)
        if(present(wind_gen)) &
            wind_gen =  solid_rotation_wind_field_t(u0=u0,alpha=alpha)
    case("isoterm")
        T0_loc = 300.0_8
        if(present(T0)) T0_loc = T0

        if(present(theta_gen)) &
            theta_gen = solid_rotation_theta_isoT_t(u0=u0,omega=omega,sphere_rad=sphere_rad,&
                                                  T0=T0_loc,grav=grav,alpha=alpha,pref=pref_loc,p0=p0_loc)
        if(present(PExner_gen)) &
            PExner_gen = solid_rotation_PExner_isoT_t(u0=u0,omega=omega,sphere_rad=sphere_rad,&
                                                   T0=T0_loc,grav=grav,alpha=alpha,pref=pref_loc,p0=p0_loc)
        if(present(wind_gen)) &
            wind_gen =  solid_rotation_wind_field_t(u0=u0,alpha=alpha)
    case default
        call parcomm_global%abort("create_solid_rotation_field_generators error, "//&
                                  "unknown background_type: "// background_type)
    end select
end subroutine create_solid_rotation_field_generators

end module solid_rotation_fields_factory_mod
