module halo_factory_mod

implicit none

contains

subroutine create_halo_procedure(halo,domain,halo_width,halo_type)
    use halo_mod,   only : halo_t
    use domain_mod, only : domain_t
    use ecs_halo_factory_mod, only : create_ecs_o_scalar_halo, create_ecs_xy_scalar_halo

    class(halo_t), allocatable, intent(out) :: halo
    class(domain_t),            intent(in)  :: domain
    integer(kind=4),            intent(in)  :: halo_width
    character(len=*),           intent(in)  :: halo_type

    if(halo_type=="A_default") then
        call create_default_exchange_halo_procedure(halo,domain,halo_width,points_type="o")
    else if(halo_type == "xy_default") then
        call create_default_exchange_halo_procedure(halo,domain,halo_width,points_type="xy")
    else if(halo_type == "xyz_default") then
        call create_default_exchange_halo_procedure(halo,domain,halo_width,points_type="xyz")
    else if(halo_type == "ECS_O") then
        call create_ecs_o_scalar_halo(halo,domain,halo_width,is_z_interfaces=.false.)
    else if(halo_type == "ECS_Oz") then
        call create_ecs_o_scalar_halo(halo,domain,halo_width,is_z_interfaces=.true.)
    else if(halo_type == "ECS_xy") then
        call create_ecs_xy_scalar_halo(halo,domain,halo_width,is_z_interfaces=.false.)
    else if(halo_type == "ECS_xyz") then
        call create_ecs_xy_scalar_halo(halo,domain,halo_width,is_z_interfaces=.true.)
    else if(halo_type == "Ah_scalar_sync") then
        call create_Ah_scalar_sync_halo_procedure(halo,.false.,domain,halo_width)
    else if(halo_type == "Ah_scalar_sync_z") then
        call create_Ah_scalar_sync_halo_procedure(halo,.true.,domain,halo_width)
    else
        call domain%parcomm%abort("unknown halo_type in create_halo_procedure: "// &
                                   halo_type)
    end if
end

subroutine create_vector_halo_procedure(halo,domain,halo_width,halo_type)
    use halo_mod,   only : halo_vec_t
    use domain_mod, only : domain_t
    use ecs_halo_vec_a_factory_mod,       only : create_ecs_A_vec_halo_procedure
    use ecs_halo_vec_xy_factory_mod,      only : create_ecs_Ah_vec_halo_procedure
    use ecs_halo_vec_c_factory_mod,       only : create_ecs_C_vec_halo_procedure, &
                                                 create_ecs_C_vec_covariant_halo_procedure
    use ecs_halo_Ah_vec_sync_factory_mod, only : create_ecs_Ah_vec_sync
    use ecs_halo_vec_Ch_factory_mod,      only : create_ecs_halo_vec_Ch

    class(halo_vec_t), allocatable, intent(out) :: halo
    class(domain_t),                intent(in)  :: domain
    integer(kind=4),                intent(in)  :: halo_width
    character(len=*),               intent(in)  :: halo_type

    if(halo_type=="A_vec_default") then
        call create_vec_default_exchange_halo_procedure(halo,domain,halo_width,points_type="o")
    else if(halo_type=="C_vec_default") then
        call create_C_vec_default_halo_procedure(halo,domain,halo_width)
    else if(halo_type=="ecs_A_vec") then
        call create_ecs_A_vec_halo_procedure(halo,domain,halo_width,is_covariant=.false.)
    else if(halo_type=="ecs_A_vec_covariant") then
        call create_ecs_A_vec_halo_procedure(halo,domain,halo_width,is_covariant=.true.)
    else if(halo_type=="ecs_C_vec") then
        call create_ecs_C_vec_halo_procedure(halo,domain,halo_width)
    else if(halo_type=="ecs_C_vec_covariant") then
        call create_ecs_C_vec_covariant_halo_procedure(halo,domain,halo_width)
    else if(halo_type=="ecs_Ah_vec" .or. halo_type == "ecs_xy_vec") then
            call create_ecs_Ah_vec_halo_procedure(halo,domain,halo_width,is_covariant=.false.)
    else if(halo_type=="ecs_Ah_vec_covariant" .or. halo_type == "ecs_xy_vec_covariant") then
            call create_ecs_Ah_vec_halo_procedure(halo,domain,halo_width,is_covariant=.true.)
    else if(halo_type=="ecs_Ah_vec_sync_contra") then
        call create_ecs_Ah_vec_sync(halo,domain,halo_width,"contravariant")
    else if(halo_type=="ecs_Ah_vec_sync_covariant") then
        call create_ecs_Ah_vec_sync(halo,domain,halo_width,"covariant")
    else if(halo_type=="ecs_Ch_vec") then
        call create_ecs_halo_vec_Ch(halo,domain,halo_width)
    else
        call domain%parcomm%abort("unknown halo_type in create_vector_halo_procedure: "// &
                                   halo_type)
    end if
end

subroutine create_default_exchange_halo_procedure(halo,domain,halo_width,points_type)
    use halo_mod,                  only : halo_t
    use domain_mod,                only : domain_t
    use halo_default_exchange_mod, only : halo_default_exchange_t
    use exchange_factory_mod,      only : create_o_points_halo_exchange, create_xy_points_halo_exchange, create_xyz_points_halo_exchange

    class(halo_t), allocatable, intent(out) :: halo
    class(domain_t),            intent(in)  :: domain
    integer(kind=4),            intent(in)  :: halo_width
    character(len=*),           intent(in)  :: points_type

    type(halo_default_exchange_t), allocatable :: halo_loc

    allocate(halo_loc)
    if(points_type == "o") then
        halo_loc%exch_halo = create_o_points_halo_exchange( &
                   domain%partition, domain%parcomm, domain%topology, halo_width, 'full')
    else if(points_type == "xy") then
        halo_loc%exch_halo = create_xy_points_halo_exchange( &
                    domain%partition, domain%parcomm, domain%topology, halo_width, 'full')
    else if(points_type == "xyz") then
        halo_loc%exch_halo = create_xyz_points_halo_exchange( &
                    domain%partition, domain%parcomm, domain%topology, halo_width, 'full')
    else
        call domain%parcomm%abort("create_default_exchange_halo_procedure error - unknown points type: "//points_type)
    end if

    call move_alloc(halo_loc, halo)

end subroutine 

subroutine create_Ah_scalar_sync_halo_procedure(halo,is_z_interfaces,domain,halo_width)
    use halo_mod,                only : halo_t
    use domain_mod,              only : domain_t
    use halo_Ah_scalar_sync_mod, only : halo_Ah_scalar_sync_t
    use exchange_factory_mod,    only : create_xy_points_halo_exchange, &
                                        create_xyz_points_halo_exchange

    class(halo_t), allocatable, intent(out) :: halo
    logical,                    intent(in)  :: is_z_interfaces
    class(domain_t),            intent(in)  :: domain
    integer(kind=4),            intent(in)  :: halo_width

    allocate(halo_Ah_scalar_sync_t :: halo)
    select type(halo)
    type is (halo_Ah_scalar_sync_t)
        halo%is_z_interfaces = is_z_interfaces
        if(is_z_interfaces) then
            halo%exch_halo = create_xyz_points_halo_exchange( &
                   domain%partition, domain%parcomm, domain%topology, 1, 'full')
        else
            halo%exch_halo = create_xy_points_halo_exchange( &
                   domain%partition, domain%parcomm, domain%topology, 1, 'full')
        end if
    end select
end

subroutine create_vec_default_exchange_halo_procedure(halo,domain,halo_width,points_type)
    use halo_mod,                  only : halo_vec_t
    use domain_mod,                only : domain_t
    use halo_default_exchange_mod, only : halo_vec_default_exchange_t
    use exchange_factory_mod,      only : create_o_points_halo_exchange

    class(halo_vec_t), allocatable, intent(out) :: halo
    class(domain_t),            intent(in)      :: domain
    integer(kind=4),            intent(in)      :: halo_width
    character(len=*),           intent(in)      :: points_type

    type(halo_vec_default_exchange_t), allocatable :: halo_loc

    allocate(halo_loc)

    if(points_type == "o") then

        halo_loc%exch_halo = create_o_points_halo_exchange(domain%partition, domain%parcomm, domain%topology, &
                                                     halo_width, 'full')
    else
        call domain%parcomm%abort("create_A_vec_default_halo_procedure error - unknown points_type: "// points_type)
    end if

    call move_alloc(halo_loc, halo)

end subroutine

subroutine create_C_vec_default_halo_procedure(halo_out,domain,halo_width)
    use halo_mod,               only : halo_vec_t
    use domain_mod,             only : domain_t
    use halo_C_default_mod,     only : halo_C_vec_default_t
    use exchange_factory_mod,   only : create_symmetric_halo_vec_exchange_C

    class(halo_vec_t), allocatable, intent(out) :: halo_out
    class(domain_t),            intent(in)      :: domain
    integer(kind=4),            intent(in)      :: halo_width

    type(halo_C_vec_default_t), allocatable :: halo
    integer(kind=4) :: t

    allocate(halo)
    !select type(halo)
    !type is (halo_C_vec_default_t)
        halo%exch_halo = create_symmetric_halo_vec_exchange_C(domain%partition, &
                                       domain%parcomm, domain%topology, halo_width, 'full')

        halo%ts_u = domain%mesh_u%ts; halo%te_u = domain%mesh_u%te;
        allocate(halo%is_left_edge(halo%ts_u:halo%te_u))
        allocate(halo%is_right_edge(halo%ts_u:halo%te_u))
        halo%ts_v = domain%mesh_v%ts; halo%te_v = domain%mesh_v%te;
        allocate(halo%is_top_edge(halo%ts_v:halo%te_v))
        allocate(halo%is_bottom_edge(halo%ts_v:halo%te_v))

        do t = halo%ts_u,halo%te_u
            halo%is_left_edge(t)  = (domain%mesh_u%tile(t)%is == 1)
            halo%is_right_edge(t) = (domain%mesh_u%tile(t)%ie == domain%partition%nh+1)
        end do
        do t = halo%ts_v,halo%te_v
            halo%is_bottom_edge(t) = (domain%mesh_v%tile(t)%js == 1)
            halo%is_top_edge(t)    = (domain%mesh_v%tile(t)%je == domain%partition%nh+1)
        end do

    !end select
    call move_alloc(halo,halo_out)
end

end module halo_factory_mod
