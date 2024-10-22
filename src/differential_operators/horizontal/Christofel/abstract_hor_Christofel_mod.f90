module abstract_hor_Christofel_mod
    use grid_field_mod, only : grid_field_t
    use mesh_mod,       only : mesh_t
    use parcomm_mod,    only : parcomm_global

    implicit none

    type, abstract :: hor_Christofel_t
    contains
        procedure :: add_contra_wind_advection
        procedure :: add_vector_advection
        procedure :: add_cov_wind_advection
        generic   :: add => add_contra_wind_advection, add_vector_advection
    end type

contains

    subroutine add_vector_advection(this,tend,s1,s2,ut,vt,mesh,component)
        class(hor_Christofel_t), intent(inout) :: this
        type(grid_field_t),      intent(inout) :: tend
        type(grid_field_t),      intent(in)    :: s1, s2, ut, vt
        type(mesh_t),            intent(in)    :: mesh
        character(len=*),        intent(in)    :: component

        call parcomm_global%abort("Christofel, add_vector_advection not implemented")
    end subroutine add_vector_advection

    subroutine add_contra_wind_advection(this,tend,ut,vt,mesh,component)
        class(hor_Christofel_t), intent(inout) :: this
        type(grid_field_t),      intent(inout) :: tend
        type(grid_field_t),      intent(in)    :: ut, vt
        type(mesh_t),            intent(in)    :: mesh
        character(len=*),        intent(in)    :: component

        call parcomm_global%abort("Christofel, add_contra_wind_advection not implemented")
    end subroutine add_contra_wind_advection

    subroutine add_cov_wind_advection(this,tend,ut,vt,mesh,component)
        class(hor_Christofel_t), intent(inout) :: this
        type(grid_field_t),      intent(inout) :: tend
        type(grid_field_t),      intent(in)    :: ut, vt
        type(mesh_t),            intent(in)    :: mesh
        character(len=*),        intent(in)    :: component

        call parcomm_global%abort("Christofel, add_cov_wind_advection not implemented")
    end subroutine add_cov_wind_advection

end module abstract_hor_Christofel_mod
