module partition_factory_mod

use partition_mod, only : partition_t

implicit none

contains

! subroutine create_partition(partition, Nh, Nz, num_panels, myid, Np, staggering_type, strategy)
!
!     type(partition_t),  intent(out) :: partition
!     integer(kind=4),    intent(in)  :: Nh         ! num of points in hor direction
!     integer(kind=4),    intent(in)  :: Nz         ! num of points in z direction
!     integer(kind=4),    intent(in)  :: num_panels !
!     integer(kind=4),    intent(in)  :: myid, Np   ! myid and num of processors
!     character(*),       intent(in)  :: staggering_type, strategy
!
!     integer(kind=4) :: t, Ntiles
!
!     partition%num_panels = num_panels
!
!     allocate(this%tile_o(num_panels*num_tiles))
!     allocate(this%tile_x(num_panels*num_tiles))
!     allocate(this%tile_y(num_panels*num_tiles))
!     allocate(this%tile_xy(num_panels*num_tiles))
!     allocate(this%tile(num_panels*num_tiles))
!     allocate(this%proc_map(num_panels*num_tiles))
!     allocate(this%panel_map(num_panels*num_tiles))
!     this%num_tiles = num_tiles
!     this%Nh = Nh
!     this%Nz = Nz
!
! end subroutine create_partition

subroutine extract_2d_partition(partition_2d, partition)

    type(partition_t), intent(out) :: partition_2d
    type(partition_t), intent(in)  :: partition

    integer(kind=4) :: t, Nt, Nh, Nz, Nz_inter, ts, te

    partition_2d%proc_map   = partition%proc_map
    partition_2d%panel_map  = partition%panel_map
    partition_2d%Nh         = partition%Nh
    partition_2d%Nz         = 1
    partition_2d%num_tiles  = partition%num_tiles
    partition_2d%num_panels = partition%num_panels
    partition_2d%Nt         = partition%Nt
    partition_2d%ts         = partition%ts
    partition_2d%te         = partition%te
    partition_2d%ts_global  = partition%ts_global
    partition_2d%te_global  = partition%te_global
    
    call extract_2d_tiles(partition_2d%tiles_o,   partition%tiles_o)
    call extract_2d_tiles(partition_2d%tiles_x,   partition%tiles_x)
    call extract_2d_tiles(partition_2d%tiles_y,   partition%tiles_y)
    call extract_2d_tiles(partition_2d%tiles_xy,  partition%tiles_xy)
    call extract_2d_tiles(partition_2d%tiles_z,   partition%tiles_z)
    call extract_2d_tiles(partition_2d%tiles_xz,  partition%tiles_xz)
    call extract_2d_tiles(partition_2d%tiles_yz,  partition%tiles_yz)
    call extract_2d_tiles(partition_2d%tiles_xyz, partition%tiles_xyz)
    call extract_2d_tiles(partition_2d%tiles_u,   partition%tiles_u)
    call extract_2d_tiles(partition_2d%tiles_v,   partition%tiles_v)
    call extract_2d_tiles(partition_2d%tiles_p,   partition%tiles_p)

end subroutine extract_2d_partition

subroutine extract_2d_tiles(tiles_2d, tiles)

    use tiles_mod, only : tiles_t

    type(tiles_t), intent(out) :: tiles_2d
    type(tiles_t), intent(in)  :: tiles

    integer(kind=4) :: t

    tiles_2d = tiles
    tiles_2d%Nz = 1

    do t = 1, tiles_2d%Nt
        tiles_2d%tile(t)%ks = 1
        tiles_2d%tile(t)%ke = 1
    end do

end subroutine extract_2d_tiles

end module partition_factory_mod
