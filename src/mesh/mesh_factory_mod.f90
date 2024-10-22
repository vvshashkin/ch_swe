module mesh_factory_mod

use mesh_mod,       only : mesh_t, tile_mesh_t
use tiles_mod,      only : tiles_t
use parcomm_mod,    only : parcomm_global
use orography_mod,  only : orography_t
use partition_mod,  only : partition_t
use metric_mod,     only : metric_t
use grid_field_mod, only : tile_field_t

implicit none

contains

subroutine create_mesh(mesh, partition, metric, halo_width, tiles, &
                       shift_xyz, orography)

    type(partition_t),       intent(in)   :: partition
    class(metric_t),         intent(in)   :: metric
    type(mesh_t),            intent(out)  :: mesh

    integer(kind=4),         intent(in)   :: halo_width
    type(tiles_t),           intent(in)   :: tiles
    real(kind=8),            intent(in)   :: shift_xyz(3)
    type(orography_t),       intent(in), optional :: orography

    integer(kind=4) :: t, pind, i, j, k, ts, te, is, ie, js, je, ks, ke, nh, nx, ny, nz
    real(kind=8)    :: hx, hy, hz

    real(kind=8) :: shift_i, shift_j, shift_k

    ! define horizontal grid parameters
    shift_i = shift_xyz(1); shift_j = shift_xyz(2); shift_k = shift_xyz(3)

    nx = tiles%Nx; ny = tiles%Ny; nz = tiles%Nz
    hx = (metric%alpha1 - metric%alpha0)/(real(nx,8)-1+2.0_8*shift_i)
    hy = (metric%beta1  - metric%beta0 )/(real(ny,8)-1+2.0_8*shift_j)
    hz = 1.0_8 / max(1.0_8,real(nz,8)-1+2.0_8*shift_k)

    ts = partition%ts
    te = partition%te

    allocate(mesh%tile(ts:te))
    mesh%ts = ts
    mesh%te = te

    mesh%scale          = metric%scale
    mesh%vertical_scale = metric%vertical_scale
    mesh%omega          = metric%omega
    mesh%rotation_axis  = metric%rotation_axis

    do t = ts, te

        call tiles%tile(t)%getind(is, ie, js, je, ks, ke)
        pind = partition%panel_map(t)

        call mesh%tile(t)%init(is, ie, js, je, ks, ke, halo_width)

        mesh%tile(t)%nx = nx
        mesh%tile(t)%ny = ny
        mesh%tile(t)%nz = nz

        mesh%tile(t)%hx = hx
        mesh%tile(t)%hy = hy
        mesh%tile(t)%hz = hz

        mesh%tile(t)%shift_i = shift_i
        mesh%tile(t)%shift_j = shift_j
        mesh%tile(t)%shift_k = shift_k

        mesh%tile(t)%alpha_0 = metric%alpha0
        mesh%tile(t)%beta_0  = metric%beta0

        mesh%tile(t)%panel_ind = pind

    end do

    call metric%set_curvilinear_mesh(mesh, orography)

end subroutine create_mesh

end module mesh_factory_mod
