module ecs_halo_vec_Ch_mod

use halo_mod,                     only : halo_vec_t
use domain_mod,                   only : domain_t
use grid_field_mod,               only : grid_field_t
use exchange_abstract_mod,        only : exchange_t
use abstract_particle_interp_mod, only : particle_interp_container_t
use particle_values_mod,          only : particle_values_t
use array_tools_mod,              only : array_container_1d_r8_t
use assembly_particles_mod,       only : assembly_particles_t

implicit none

type, public, extends(halo_vec_t) :: ecs_halo_vec_Ch_t

  class(exchange_t), allocatable :: exchange
  class(particle_interp_container_t), allocatable :: interp_u(:), interp_v(:)
  type(particle_values_t), allocatable :: values_v(:), values_u(:), values_assembled(:)

  type(array_container_1d_r8_t), allocatable :: u_proj(:), v_proj(:)

  type(assembly_particles_t) :: assembly
contains
  procedure :: get_halo_vector => get_halo_vector

end type

contains

subroutine get_halo_vector(this, u, v, domain, halo_width)

    class(ecs_halo_vec_Ch_t), intent(inout) :: this
    class(grid_field_t),      intent(inout) :: u
    class(grid_field_t),      intent(inout) :: v
    type(domain_t),           intent(in)    :: domain
    integer(kind=4),          intent(in)    :: halo_width

    integer(kind=4) :: t, i, j, k, ic, nx, ny

    call this%exchange%do_vec(u, v, domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call this%interp_u(t)%interp%interp_2d_levs(this%values_u(t), u%tile(t), domain)
        call this%interp_v(t)%interp%interp_2d_levs(this%values_v(t), v%tile(t), domain)

        do k = 1, this%values_u(t)%Nfld
          do i = 1, this%values_u(t)%N
            this%values_u(t)%p(i,k) = this%u_proj(t)%a(i)*this%values_u(t)%p(i,k) + &
                                      this%v_proj(t)%a(i)*this%values_v(t)%p(i,k) 
          end do
        end do

    end do

    call this%assembly%do(this%values_assembled, this%values_u, domain)

    do t = domain%partition%ts, domain%partition%te

      do k = domain%mesh_x%tile(t)%ks, domain%mesh_x%tile(t)%ke

        ny = domain%mesh_y%tile(t)%ny
        nx = domain%mesh_x%tile(t)%nx

        ic = 1
        !'south' panel edge
        if(domain%mesh_y%tile(t)%js == 1) then
            do j = 1, halo_width
                do i = domain%mesh_y%tile(t)%is, domain%mesh_y%tile(t)%ie
                  u%tile(t)%p(i,1-j,k) = this%values_assembled(t)%p(ic,k)
                  ic = ic+1
                end do
            end do
        end if

      !'north' panel edge
      if(domain%mesh_y%tile(t)%je == ny) then
          do j = 1, halo_width
              do i = domain%mesh_y%tile(t)%is, domain%mesh_y%tile(t)%ie
                u%tile(t)%p(i,ny+j,k) = this%values_assembled(t)%p(ic,k)
                ic = ic+1
              end do
          end do
      end if

      !'west' panel edge
      if(domain%mesh_x%tile(t)%is == 1) then
          do i = 1, halo_width
              do j = domain%mesh_x%tile(t)%js, domain%mesh_x%tile(t)%je
                v%tile(t)%p(1-i,j,k) = this%values_assembled(t)%p(ic,k)
                ic = ic+1
              end do
          end do
      end if

      !'east' panel edge
      if(domain%mesh_x%tile(t)%ie == nx) then
          do i = 1, halo_width
              do j = domain%mesh_x%tile(t)%js, domain%mesh_x%tile(t)%je
                v%tile(t)%p(nx+i,j,k) = this%values_assembled(t)%p(ic,k)
                ic = ic+1
              end do
          end do
      end if

      ny = domain%mesh_x%tile(t)%ny
      nx = domain%mesh_y%tile(t)%nx

      !'south' panel edge
      if(domain%mesh_x%tile(t)%js == 1) then
          do j = 1, halo_width
              do i = domain%mesh_x%tile(t)%is, domain%mesh_x%tile(t)%ie
                v%tile(t)%p(i,1-j,k) = this%values_assembled(t)%p(ic,k)
                ic = ic+1
              end do
          end do
      end if

    !'north' panel edge
    if(domain%mesh_x%tile(t)%je == ny) then
        do j = 1, halo_width
            do i = domain%mesh_x%tile(t)%is, domain%mesh_x%tile(t)%ie
              v%tile(t)%p(i,ny+j,k) = this%values_assembled(t)%p(ic,k)
              ic = ic+1
            end do
        end do
    end if

    !'west' panel edge
    if(domain%mesh_y%tile(t)%is == 1) then
        do i = 1, halo_width
            do j = domain%mesh_y%tile(t)%js, domain%mesh_y%tile(t)%je
              u%tile(t)%p(1-i,j,k) = this%values_assembled(t)%p(ic,k)
              ic = ic+1
            end do
        end do
    end if

    !'east' panel edge
    if(domain%mesh_y%tile(t)%ie == nx) then
        do i = 1, halo_width
            do j = domain%mesh_y%tile(t)%js, domain%mesh_y%tile(t)%je
              u%tile(t)%p(nx+i,j,k) = this%values_assembled(t)%p(ic,k)
              ic = ic+1
            end do
        end do
    end if

    end do
  end do

end subroutine get_halo_vector

end module ecs_halo_vec_Ch_mod