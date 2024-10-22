module regrid_factory_mod

use parcomm_mod,            only : parcomm_global
use abstract_regrid_mod,    only : regrid_t, regrid_vec_t
use latlon_regrid_mod,      only : latlon_regrid_t, latlon_regrid_vec_t, &
                                   indices_container_t, interp_t
use domain_mod,             only : domain_t
use halo_mod,               only : halo_t, halo_vec_t
use halo_factory_mod,       only : create_halo_procedure, create_vector_halo_procedure
use grid_field_factory_mod, only : create_grid_field
use mesh_mod,               only : mesh_t
use const_mod,              only : pi
use tiles_mod,              only : tiles_t

use interpolator2d_factory_mod, only : create_vec2vec_interpolator2d
use latlon_grid_generator_mod,  only : latlon_grid_generator_t

implicit none

type, private :: int_2d_array_t
    integer(kind=4), allocatable :: a(:,:)
end type

contains

subroutine create_latlon_regrid(regrid_out, domain, latlon_grid, interp_type, &
                                mesh_name)

    class(regrid_t), allocatable,  intent(out) :: regrid_out
    type(domain_t),                intent(in)  :: domain
    type(latlon_grid_generator_t), intent(in)  :: latlon_grid
    character(len=*),              intent(in)  :: interp_type
    character(len=*),              intent(in)  :: mesh_name

    type(latlon_regrid_t), allocatable :: regrid

    real(kind=8), allocatable :: my_alpha(:,:), my_beta(:,:) !alpha & beta of latlon points belonging to this process

    integer(kind=4) :: t

    type(mesh_t), pointer :: mesh
                
    allocate(regrid)


    call distribute_latlon_points_to_cubsphere(my_alpha, my_beta, &
                     regrid%msg_to_latlon, latlon_grid, domain, mesh_name)

    call init_halo_for_interp(regrid%halo, regrid%halo_width, mesh_name, &
                              interp_type, domain)

    !WORKAROUND: again halo_width == 8
    call domain%get_mesh(mesh, mesh_name)
    call create_grid_field(regrid%buff, 8, 0, mesh)

    allocate(regrid%latlon_interp(mesh%ts:mesh%te))
    do t = mesh%ts, mesh%te
        call regrid%latlon_interp(t)%init(my_alpha(:,t), my_beta(:,t), &
                            regrid%msg_to_latlon(t)%n, mesh%tile(t), interp_type)
    end do

    regrid%mesh_name  = mesh_name

    allocate(regrid%values(mesh%ts:mesh%te))
    do t = mesh%ts, mesh%te
        allocate(regrid%values(t)%a(regrid%msg_to_latlon(t)%n,mesh%tile(t)%nz))
    end do

    regrid%nz = mesh%tile(mesh%ts)%nz

    call move_alloc(regrid, regrid_out)

end subroutine create_latlon_regrid

subroutine create_latlon_vector_regrid(regrid_out, domain, latlon_grid, interp_type, &
                                       mesh_u_name, mesh_v_name, components_type)

    class(regrid_vec_t), allocatable, intent(out) :: regrid_out
    type(domain_t),                   intent(in)  :: domain
    type(latlon_grid_generator_t),    intent(in)  :: latlon_grid
    character(len=*),                 intent(in)  :: interp_type
    character(len=*),                 intent(in)  :: mesh_u_name, mesh_v_name
    character(len=*),                 intent(in)  :: components_type

    type(latlon_regrid_vec_t), allocatable :: regrid

    real(kind=8), allocatable :: my_alpha(:,:), my_beta(:,:) !alpha & beta of latlon points belonging to this process
    real(kind=8), allocatable :: my_lon(:,:), my_lat(:,:) !lon and lat of latlon points belonging to this process

    integer(kind=4) :: t

    type(mesh_t), pointer :: mesh

    character(len=:), allocatable :: mesh_interp_name
    
    allocate(regrid)

    if(mesh_u_name == "xy" .and. mesh_v_name == "xy") then
        mesh_interp_name = "xy"
    else if(mesh_u_name == "o" .and. mesh_v_name == "o") then
        mesh_interp_name = "o"
    else if(mesh_u_name == "x" .and. mesh_v_name == "y") then
        mesh_interp_name = "o"
        call create_vec2vec_interpolator2d(regrid%staggered2colocated, "interp2d_uv2pvec_C_sbp42", domain)
        call domain%get_mesh(regrid%mesh_u, mesh_u_name)
        call domain%get_mesh(regrid%mesh_v, mesh_v_name)
        call create_grid_field(regrid%buff_u_stag, 4, 0, regrid%mesh_u)
        call create_grid_field(regrid%buff_v_stag, 4, 0, regrid%mesh_v)
    else if(mesh_u_name == "y" .and. mesh_v_name == "x") then
        mesh_interp_name = "xy"
        call create_vec2vec_interpolator2d(regrid%staggered2colocated, "interp2d_uv2pvec_Ch_sbp42", domain)
        call domain%get_mesh(regrid%mesh_u, mesh_u_name)
        call domain%get_mesh(regrid%mesh_v, mesh_v_name)
        call create_grid_field(regrid%buff_u_stag, 4, 0, regrid%mesh_u)
        call create_grid_field(regrid%buff_v_stag, 4, 0, regrid%mesh_v)
    else
        call parcomm_global%abort("create_latlon_vector_regrid error - unsupported combination of mesh_u and mesh_v: "// mesh_u_name //", "//mesh_v_name)
    end if

    call domain%get_mesh(regrid%mesh, mesh_interp_name)

    call distribute_latlon_points_to_cubsphere(my_alpha, my_beta, &
                                               regrid%msg_to_latlon, latlon_grid,   &
                                               domain, mesh_interp_name)

    call init_halo_for_interp_vec(regrid%halo, regrid%halo_width, &
                                  mesh_interp_name, interp_type,  &
                                  components_type, domain)

    allocate(regrid%latlon_interp(regrid%mesh%ts:regrid%mesh%te))
    do t = regrid%mesh%ts, regrid%mesh%te
        call regrid%latlon_interp(t)%init(my_alpha(:,t), my_beta(:,t), &
                                          regrid%msg_to_latlon(t)%n, regrid%mesh%tile(t), &
                                          interp_type)
    end do

    !WORKAROUND: again halo_width == 8
    call create_grid_field(regrid%buff_u, 8, 0, regrid%mesh)
    call create_grid_field(regrid%buff_v, 8, 0, regrid%mesh)

    allocate(regrid%values_uv(regrid%mesh%ts:regrid%mesh%te))
    allocate(regrid%q(regrid%mesh%ts:regrid%mesh%te))
    do t = regrid%mesh%ts, regrid%mesh%te
        allocate(regrid%values_uv(t)%a(regrid%msg_to_latlon(t)%n,regrid%mesh%tile(t)%nz,2))
        call init_transform_to_latlon_vec_tile(regrid%q(t)%a, my_alpha(:,t), my_beta(:,t), &
                                               latlon_grid, regrid%msg_to_latlon(t),                    &
                                               components_type, domain,                    &
                                               domain%partition%panel_map(t))
    end do

    regrid%nz = regrid%mesh%tile(regrid%mesh%ts)%nz

    call move_alloc(regrid, regrid_out)

end subroutine create_latlon_vector_regrid

subroutine distribute_latlon_points_to_cubsphere(my_alpha, my_beta, &
                              msg_to_latlon, latlon_grid, domain, mesh_name)

    real(kind=8),              allocatable, intent(out) :: my_alpha(:,:), &
                                                           my_beta(:,:) !alpha & beta of latlon points belonging to this process
    type(indices_container_t), allocatable, intent(out) :: msg_to_latlon(:)

    type(latlon_grid_generator_t), intent(in) :: latlon_grid
    type(domain_t),   intent(in)              :: domain
    character(len=*), intent(in)              :: mesh_name

    type(mesh_t), pointer :: mesh     
    type(tiles_t) :: tiles

    type(int_2d_array_t) :: tile_map(domain%partition%num_panels)

    integer(kind=4) :: n_points_at_tile(domain%partition%num_panels*domain%partition%num_tiles)

    integer(kind=4) :: i, j, t, panel_ind, size, indx, indy
    real(kind=8) :: alpha, beta, alpha0, beta0, shift_alpha, shift_beta
    real(kind=8) :: r(3), zdx, zdy, hx, hy
    
    !get mesh and tiles description by name
    call domain%get_mesh(mesh, mesh_name)
    call domain%partition%get_tiles(mesh_name, tiles)

    !prepare map of tiles: panel_index, i, j -> tile_number
    do panel_ind = 1, domain%partition%num_panels
        allocate(tile_map(panel_ind)%a(tiles%Nx,tiles%Ny))
    end do
    do t = 1, tiles%Nt
        panel_ind = domain%partition%panel_map(t)
        tile_map(panel_ind)%a(tiles%tile(t)%is:tiles%tile(t)%ie,tiles%tile(t)%js:tiles%tile(t)%je) = t
    end do

    !prepare interpolated values message map:
    if(domain%parcomm%myid == 0) then
        allocate(msg_to_latlon(tiles%Nt))
    else
        allocate(msg_to_latlon(mesh%ts:mesh%te))      
    end if

    do t = lbound(msg_to_latlon,1), ubound(msg_to_latlon,1)
        size = latlon_grid%Nlon*latlon_grid%Nlat / 2 !upper estimate of latlon points number at each panel
        allocate(msg_to_latlon(t)%i(size))
        allocate(msg_to_latlon(t)%j(size))
        msg_to_latlon(t)%n = 0
    end do

    !distribute points for tiles
    allocate(my_alpha(latlon_grid%Nlon*latlon_grid%Nlat/2,mesh%ts:mesh%te))
    allocate(my_beta (latlon_grid%Nlon*latlon_grid%Nlat/2,mesh%ts:mesh%te))

    alpha0 = mesh%tile(mesh%ts)%alpha_0
    beta0  = mesh%tile(mesh%ts)%beta_0
    hx = mesh%tile(mesh%ts)%hx
    hy = mesh%tile(mesh%ts)%hy
    shift_alpha = mesh%tile(mesh%ts)%shift_i
    shift_beta  = mesh%tile(mesh%ts)%shift_j

    do j = 1, latlon_grid%Nlat
        do i = 1, latlon_grid%Nlon

            call latlon_grid%get_cartesian_coords(r,i,j)
            call domain%metric%transform_cartesian_to_native(panel_ind, alpha, beta, r)

            zdx = (alpha-alpha0)/hx+1-shift_alpha
            zdy = (beta -beta0) /hy+1-shift_beta
            indx = max(1,min(mesh%tile(mesh%ts)%nx,int(zdx)))
            indy = max(1,min(mesh%tile(mesh%ts)%ny,int(zdy)))

            t = tile_map(panel_ind)%a(indx, indy)

            if(domain%parcomm%myid == 0 .or. t >= mesh%ts .and. t <= mesh%te) then

                msg_to_latlon(t)%n = msg_to_latlon(t)%n + 1

                msg_to_latlon(t)%i(msg_to_latlon(t)%n) = i
                msg_to_latlon(t)%j(msg_to_latlon(t)%n) = j
                
            end if

            if (t >= mesh%ts .and. t <= mesh%te) then

                my_alpha(msg_to_latlon(t)%n,t) = alpha
                my_beta(msg_to_latlon(t)%n,t) = beta

            end if

        end do
    end do

    shrink : block

        use array_tools_mod, only : shrink_array

        do t = lbound(msg_to_latlon,1), ubound(msg_to_latlon,1)

            call shrink_array(msg_to_latlon(t)%i, msg_to_latlon(t)%n)
            call shrink_array(msg_to_latlon(t)%j, msg_to_latlon(t)%n)

        end do

    end block shrink

end subroutine

subroutine init_halo_for_interp(halo, halo_width, mesh_name, &
                                interp_type, domain)

    class(halo_t),  allocatable, intent(out) :: halo
    integer(kind=4),             intent(out) :: halo_width
    character(len=*),            intent(in)  :: mesh_name, interp_type
    type(domain_t),              intent(in)  :: domain

    integer(kind=4) :: t, t1
    type(mesh_t), pointer :: mesh

    if(mesh_name == "xy" .and. interp_type == "linear") then
        halo_width = 1
        call create_halo_procedure(halo,domain, halo_width, "xy_default")
    else if(mesh_name == "xy" .and. interp_type == "cubic") then
        halo_width = 2
        call create_halo_procedure(halo,domain, halo_width, "ECS_xy")
    else if(mesh_name == "xyz" .and. interp_type == "linear") then
        halo_width = 1
        call create_halo_procedure(halo,domain, halo_width, "xyz_default")
    else if(mesh_name == "xyz" .and. interp_type == "cubic") then
        halo_width = 2
        call create_halo_procedure(halo,domain, halo_width, "ECS_xyz")
    else if(mesh_name == "o" .and. (interp_type == "linear" .or. interp_type == "cubic")) then
        halo_width = 2
        call create_halo_procedure(halo,domain, halo_width, "ECS_O")
    else if(mesh_name == "z" .and. (interp_type == "linear" .or. interp_type == "cubic")) then
        halo_width = 2
        call create_halo_procedure(halo,domain, halo_width, "ECS_Oz")
    else
        call parcomm_global%abort("create_latlon_regrid_vec error - unsupported combination of mesh_name and interp_type: "// mesh_name//", "//interp_type)
    end if

end subroutine

subroutine init_halo_for_interp_vec(halo, halo_width, mesh_name, &
                                    interp_type, components_type, domain)

    class(halo_vec_t),  allocatable, intent(out) :: halo
    integer(kind=4),                 intent(out) :: halo_width

    character(len=*),            intent(in)  :: mesh_name, interp_type, components_type
    type(domain_t),              intent(in)  :: domain

    integer(kind=4) :: t, t1
    type(mesh_t), pointer :: mesh

    if(mesh_name == "xy" .and. interp_type == "linear" .and. &
       components_type == "contravariant") then
        halo_width = 1
        call create_vector_halo_procedure(halo,domain, halo_width, "ecs_xy_vec")
       else if(mesh_name == "xy" .and. interp_type == "linear" .and. &
               components_type == "covariant") then
               halo_width = 1
               call create_vector_halo_procedure(halo,domain, halo_width, "ecs_xy_vec_covariant")
    else if(mesh_name == "xy" .and. interp_type == "cubic" .and.&
            components_type == "contravariant") then
        halo_width = 2
        call create_vector_halo_procedure(halo,domain, halo_width, "ecs_xy_vec")
    else if(mesh_name == "xy" .and. interp_type == "cubic" .and.&
        components_type == "covariant") then
        halo_width = 2
        call create_vector_halo_procedure(halo,domain, halo_width, "ecs_xy_vec_covariant")
    else if(mesh_name == "o" .and. &
      (interp_type == "linear" .or. interp_type == "cubic") .and. &
       components_type == "contravariant") then
        halo_width = 2
        call create_vector_halo_procedure(halo,domain, halo_width, "ecs_A_vec")
    else if(mesh_name == "o" .and. &
      (interp_type == "linear" .or. interp_type == "cubic") .and. &
       components_type == "covariant") then
        halo_width = 2
        call create_vector_halo_procedure(halo,domain, halo_width, "ecs_A_vec_covariant")
    else
        call parcomm_global%abort("create_latlon_regrid error - unsupported combination of mesh_name, interp_type and components_type: "// mesh_name//", "//interp_type//", "//components_type)
    end if

end subroutine

subroutine init_transform_to_latlon_vec_tile(q, alpha, beta, latlon_grid, latlon_ind, &
                                  components_type, domain, panel_ind)



    real(kind=8),     intent(out), allocatable :: q(:,:,:)

    real(kind=8),                  intent(in) :: alpha(:), beta(:)
    type(latlon_grid_generator_t), intent(in) :: latlon_grid
    type(indices_container_t),     intent(in) :: latlon_ind
    character(len=*),              intent(in) :: components_type
    type(domain_t),                intent(in) :: domain
    integer(kind=4),               intent(in) :: panel_ind

    integer(kind=4) :: ihor
    real(kind=8)    :: iv(3), jv(3), v1(4), v2(4) !latlon and native basis

    allocate(q(latlon_ind%n,2,2))

    do ihor = 1, latlon_ind%n

        call latlon_grid%get_basis_vectors(iv,jv,latlon_ind%i(ihor),latlon_ind%j(ihor))

        if(components_type == "covariant") then
            v1(1:3) = domain%metric%calculate_b1(panel_ind,alpha(ihor),beta(ihor))
            v2(1:3) = domain%metric%calculate_b2(panel_ind,alpha(ihor),beta(ihor))
        else if(components_type == "contravariant") then
            v1(1:4) = domain%metric%calculate_a1(panel_ind,alpha(ihor),beta(ihor))
            v2(1:4) = domain%metric%calculate_a2(panel_ind,alpha(ihor),beta(ihor))
        else
            call parcomm_global%abort("create vector latlon regrid - unknown vector components type: "// components_type)
        end if

        q(ihor,1,1) = sum(iv(1:3)*v1(1:3))
        q(ihor,1,2) = sum(iv(1:3)*v2(1:3))
        q(ihor,2,1) = sum(jv(1:3)*v1(1:3))
        q(ihor,2,2) = sum(jv(1:3)*v2(1:3))
    end do

end subroutine

end module regrid_factory_mod
