module dep_points_interp_driver_factory_mod

use abstract_dep_points_interp_driver_mod, only : dep_points_interp_driver_t
use dp_interpolator_factory_mod,           only : create_dp_interpolator
use domain_mod,                            only : domain_t
use WD_interp_driver_mod,                  only : WD_interp_driver_t
use parcomm_mod,                           only : parcomm_global
use tiles_mod,                             only : tiles_t
use tile_mod,                              only : tile_t
use string_mod,                            only : string_t
use grid_field_factory_mod,                only : create_grid_field
use generic_config_mod,                    only : generic_config_t
use mpi

implicit none

integer(kind=4), parameter :: MIN_BUFF_ALLOC_SIZE = 10
integer(kind=4), parameter :: MAX_CFL_DEFAULT = 10
contains

subroutine create_dep_points_interp_driver(dep_points_interp_driver, domain, config)

    class(dep_points_interp_driver_t), intent(out), allocatable :: dep_points_interp_driver
    type(domain_t),                    intent(in)    :: domain
    class(generic_config_t),           intent(inout) :: config

    character(len=:), allocatable :: dep_points_interp_driver_name

    call config%get(dep_points_interp_driver_name,"dep_points_interp_driver_name")

    select case (dep_points_interp_driver_name)
    case ("WhiteDongarra")
        call create_WD_dp_interp_driver(dep_points_interp_driver, domain, config)
    case default
        call parcomm_global%abort("DP interp factory, unknown dp interp name "//dep_points_interp_driver_name)
    end select

end subroutine

subroutine create_WD_dp_interp_driver(dep_points_interp_driver, domain, config)

    class(dep_points_interp_driver_t), allocatable, intent(out)   :: dep_points_interp_driver
    type(domain_t),                       intent(in)    :: domain
    class(generic_config_t),              intent(inout) :: config

    integer(kind=4) :: i, j, is, ie, js, je, ks, ke, t, n_my_points, ierr
    type(WD_interp_driver_t), allocatable :: WD_driver
    type(tiles_t) :: tiles

    integer(kind=4) :: max_cfl
    integer(kind=4) :: send_dp_to_tiles(domain%partition%Nt)
    integer(kind=4) :: expected_dp_recv_num, last_mpi_id(domain%partition%ts:domain%partition%te)
    integer(kind=4) :: local_panel, remote_panel, local_tile, remote_tile, s_num
    type(tile_t)    :: temp_tile, intersect_tile
    class(generic_config_t), allocatable :: interp_configs(:)
    logical :: is_intersection
    type(string_t), allocatable :: unique_interp_names(:)
    character(len=:), allocatable :: interp_name, arr_mesh_name
    integer(kind=4) :: num_unique_interp

    call config%get(max_cfl,"max_cfl",default=MAX_CFL_DEFAULT)
    call config%get(interp_configs,"interp_configs")
    call config%get(arr_mesh_name,"arrival_mesh")

    allocate(WD_driver)

    call mpi_comm_dup(domain%parcomm%comm_w, WD_driver%dp_comm, ierr)
    call mpi_comm_dup(domain%parcomm%comm_w, WD_driver%val_comm, ierr)

    call domain%partition%get_tiles("xy", tiles)
    call domain%get_mesh(WD_driver%mesh,arr_mesh_name)

    allocate(WD_driver%tile_map(domain%partition%num_panels))
    do i = 1, domain%partition%num_panels
        allocate(WD_driver%tile_map(i)%ind(tiles%Nx,tiles%Ny))
        WD_driver%tile_map(i)%ind(1:tiles%Nx,1:tiles%Ny) = 0
    end do

    WD_driver%n_my_points = 0
    do t = domain%partition%ts, domain%partition%te
        is = WD_driver%mesh%tile(t)%is
        ie = WD_driver%mesh%tile(t)%ie
        js = WD_driver%mesh%tile(t)%js
        je = WD_driver%mesh%tile(t)%je
        ks = WD_driver%mesh%tile(t)%ks
        ke = WD_driver%mesh%tile(t)%ke
        WD_driver%n_my_points = WD_driver%n_my_points + (ie-is+1)*(je-js+1)*(ke-ks+1)
    end do

    WD_driver%num_sends_dp = 0
    WD_driver%num_recvs_dp = 0
    last_mpi_id(domain%partition%ts:domain%partition%te) = -1000
    do remote_tile = 1, domain%partition%Nt

        remote_panel = domain%partition%panel_map(remote_tile)

        do local_tile = domain%partition%ts, domain%partition%te

            local_panel = domain%partition%panel_map(local_tile)

            call domain%topology%transform_tile_coords(remote_panel, tiles%tile(remote_tile), &
                                                local_panel,  temp_tile,               &
                                                tiles%Nx, tiles%Ny)

            call find_tiles_halo_intersection(tiles%tile(local_tile), max_cfl, max_cfl, &
                                temp_tile, "full", intersect_tile, is_intersection)

            if(is_intersection) then
                WD_driver%num_sends_dp = WD_driver%num_sends_dp+1
                send_dp_to_tiles(WD_driver%num_sends_dp) = remote_tile

                is = tiles%tile(remote_tile)%is
                ie = tiles%tile(remote_tile)%ie
                js = tiles%tile(remote_tile)%js
                je = tiles%tile(remote_tile)%je
                WD_driver%tile_map(remote_panel)%ind(is:ie,js:je) = WD_driver%num_sends_dp
                exit
            end if
        end do

        do local_tile = domain%partition%ts, domain%partition%te

            local_panel = domain%partition%panel_map(local_tile)

            call domain%topology%transform_tile_coords(local_panel, tiles%tile(local_tile), &
                                                       remote_panel,  temp_tile,            &
                                                       tiles%Nx, tiles%Ny)

            call find_tiles_halo_intersection(tiles%tile(remote_tile), MAX_CFL, MAX_CFL, &
                                temp_tile, "full", intersect_tile, is_intersection)

            if(is_intersection .and. &
               domain%partition%proc_map(remote_tile) /= last_mpi_id(local_tile)) then
                   WD_driver%num_recvs_dp = WD_driver%num_recvs_dp + 1
                   last_mpi_id(local_tile) = domain%partition%proc_map(remote_tile)
            end if
        end do
    end do

    WD_driver%send_dp_to_tiles = send_dp_to_tiles(1:WD_driver%num_sends_dp)

    allocate(WD_driver%dp_send_buffer(3*WD_driver%n_my_points+1)) !plus 1 fake point for zero-length messages at the end of buffer
    allocate(WD_driver%dp_send_request(WD_driver%num_sends_dp))

    allocate(WD_driver%gather_map%i(WD_driver%n_my_points))
    allocate(WD_driver%gather_map%j(WD_driver%n_my_points))
    allocate(WD_driver%gather_map%k(WD_driver%n_my_points))
    allocate(WD_driver%gather_map%t(WD_driver%n_my_points))
    allocate(WD_driver%gather_map%start_of_tile(WD_driver%num_sends_dp))
    allocate(WD_driver%dp_recv_buffer(6*WD_driver%n_my_points))
    !Assume that we can receive 2 times more departure points than we have:
    allocate(WD_driver%val_buffer(2*WD_driver%n_my_points*size(interp_configs,1)))
    allocate(WD_driver%val_send_request(WD_driver%num_recvs_dp))



    allocate(WD_driver%interp_fields(size(interp_configs,1)))
    allocate(WD_driver%interpolator_index(size(interp_configs,1)))
    allocate(unique_interp_names(size(interp_configs,1)))
    num_unique_interp = 0

    do i = 1, size(interp_configs,1)

        call interp_configs(i)%get(WD_driver%interp_fields(i),"field_name")

        call interp_configs(i)%get(interp_name,"interpolation_name")
        
        do j = 1, num_unique_interp
            if(interp_name == unique_interp_names(j)%str) exit
        end do
        if(j == num_unique_interp+1) then
            unique_interp_names(j)%str = interp_name
            num_unique_interp = num_unique_interp+1
        end if

        WD_driver%interpolator_index(i) = j

    end do

    WD_driver%n_fld = size(WD_driver%interp_fields,1)
    
    allocate(WD_driver%interpolator_array(1:num_unique_interp))

    do i = 1, num_unique_interp
        call create_dp_interpolator(WD_driver%interpolator_array(i)%interpolator,&
                                    domain,unique_interp_names(i)%str)
    end do

    call move_alloc(WD_driver,dep_points_interp_driver)

end subroutine

subroutine find_tiles_halo_intersection(tile1, halo_width1_x, halo_width1_y, tile2, halo_type, out_tile, is_intersect)

    type(tile_t),     intent(in) :: tile1, tile2
    integer(kind=4),  intent(in) :: halo_width1_x, halo_width1_y
    character(len=*), intent(in) :: halo_type

    type(tile_t),    intent(out) :: out_tile
    logical,         intent(out) :: is_intersect

    is_intersect = .true.

    if (halo_type=='full') then
        out_tile%ks = tile1%ks
        out_tile%ke = tile1%ke

        out_tile%is = max(tile1%is-halo_width1_x, tile2%is)
        out_tile%ie = min(tile1%ie+halo_width1_x, tile2%ie)

        out_tile%js = max(tile1%js-halo_width1_y, tile2%js)
        out_tile%je = min(tile1%je+halo_width1_y, tile2%je)

        if ((out_tile%is<=out_tile%ie) .and. (out_tile%js<=out_tile%je)) return

    else if (halo_type == 'cross') then

        out_tile%ks = tile1%ks
        out_tile%ke = tile1%ke

        out_tile%is = max(tile1%is-halo_width1_x, tile2%is)
        out_tile%ie = min(tile1%ie+halo_width1_x, tile2%ie)

        out_tile%js = max(tile1%js, tile2%js)
        out_tile%je = min(tile1%je, tile2%je)

        if ((out_tile%is<=out_tile%ie) .and. (out_tile%js<=out_tile%je)) return

        out_tile%is = max(tile1%is, tile2%is)
        out_tile%ie = min(tile1%ie, tile2%ie)

        out_tile%js = max(tile1%js-halo_width1_y, tile2%js)
        out_tile%je = min(tile1%je+halo_width1_y, tile2%je)

        if ((out_tile%is<=out_tile%ie) .and. (out_tile%js<=out_tile%je)) return
    else
        print*, 'Error! Wrong halo_type!'
        stop
    end if

    is_intersect = .false.

end subroutine find_tiles_halo_intersection

end module