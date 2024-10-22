module vec_math_mod

use grid_field_mod, only : grid_field_t, tile_field_t
use mesh_mod,       only : mesh_t, tile_mesh_t
use parcomm_mod,    only : parcomm_t, parcomm_global
use mpi

implicit none

contains

function mass(f, mesh, parcomm) result(out)

    type(grid_field_t), intent(in) :: f
    type(mesh_t),       intent(in) :: mesh
    type(parcomm_t),    intent(in) :: parcomm
    real(kind=8)                   :: out

    integer(kind=4) :: t, err
    real(kind=8)    :: out_loc

    out_loc = 0.0_8

    do t = mesh%ts, mesh%te
        out_loc = out_loc + calc_mass_tile(f%tile(t), mesh%tile(t))
    end do

    call mpi_allreduce(out_loc, out, 1, mpi_double, mpi_sum, parcomm%comm_w, err)

end function mass
function calc_mass_tile(f, mesh) result(out)

    type(tile_field_t), intent(in) :: f
    type(tile_mesh_t),  intent(in) :: mesh
    real(kind=8)                   :: out

    integer(kind=4) :: k, j, i

    out = 0.0_8

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                out = out + f%p(i,j,k)*mesh%J(i,j,k)*mesh%hx*mesh%hy
            end do
        end do
    end do

end function calc_mass_tile
function dot_prod(f1, f2, mesh, parcomm) result(out)

    type(grid_field_t), intent(in) :: f1, f2
    type(mesh_t),       intent(in) :: mesh
    type(parcomm_t),    intent(in) :: parcomm
    real(kind=8)                   :: out

    integer(kind=4) :: t, err
    real(kind=8)    :: out_loc

    out_loc = 0.0_8

    do t = mesh%ts, mesh%te
        out_loc = out_loc + calc_dot_prod_on_tile(f1%tile(t), f2%tile(t), mesh%tile(t))
    end do

    call mpi_allreduce(out_loc, out, 1, mpi_double, mpi_sum, parcomm%comm_w, err)

end function dot_prod

function calc_dot_prod_on_tile(f1, f2, mesh) result(out)

    type(tile_field_t), intent(in) :: f1, f2
    type(tile_mesh_t),  intent(in) :: mesh
    real(kind=8)                   :: out

    integer(kind=4) :: k, j, i

    out = 0.0_8

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                out = out + f1%p(i,j,k)*f2%p(i,j,k)*mesh%J(i,j,k)*mesh%hx*mesh%hy
            end do
        end do
    end do

end function calc_dot_prod_on_tile
function l2norm(f, mesh, parcomm) result(out)

    type(grid_field_t), intent(in) :: f
    type(mesh_t),       intent(in) :: mesh
    type(parcomm_t),    intent(in) :: parcomm
    real(kind=8)                   :: out

    integer(kind=4) :: t, err
    real(kind=8)    :: out_loc

    out_loc = 0.0_8

    do t = mesh%ts, mesh%te
        out_loc = out_loc + calc_l2norm_squared_on_tile(f%tile(t), mesh%tile(t))
    end do

    call mpi_allreduce(out_loc, out, 1, mpi_double, mpi_sum, parcomm%comm_w, err)

    out = sqrt(out)

end function l2norm

function calc_l2norm_squared_on_tile(f, mesh) result(out)

    type(tile_field_t), intent(in) :: f
    type(tile_mesh_t),  intent(in) :: mesh
    real(kind=8)                   :: out

    integer(kind=4) :: k, j, i

    out = 0.0_8

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                out = out + f%p(i,j,k)**2*mesh%J(i,j,k)*mesh%hx*mesh%hy
            end do
        end do
    end do

end function calc_l2norm_squared_on_tile

subroutine cube2cart_vec(vx, vy, vz, u, v, mesh, components_type)

    type(grid_field_t), intent(in)    :: u, v
    type(grid_field_t), intent(inout) :: vx, vy, vz
    type(mesh_t),       intent(in)    :: mesh
    character(len=*),   intent(in)    :: components_type

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        if(components_type == "covariant") then
            call calc_cube2cart_vec_tile(vx%tile(t), vy%tile(t), vz%tile(t), &
                                         u%tile(t), v%tile(t), mesh%tile(t) )
        else if(components_type == "contravariant") then
            call calc_cube2cart_vec_tile_contra(vx%tile(t), vy%tile(t), vz%tile(t), &
                                                u%tile(t), v%tile(t), mesh%tile(t) )
        else
            call parcomm_global%abort("cube2cart_vec, unknown components_type: "//components_type//", file: "//__FILE__)
        end if
    end do

end subroutine cube2cart_vec

subroutine calc_cube2cart_vec_tile(vx, vy, vz, u, v, mesh)

    type(tile_field_t), intent(in)    :: u, v
    type(tile_field_t), intent(inout) :: vx, vy, vz
    type(tile_mesh_t),  intent(in)    :: mesh

    real(kind=8)    :: v_xyz(3)
    integer(kind=4) :: k, j, i

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                v_xyz(1:3) = mesh%b1(1:3,i,j,k)*u%p(i,j,k)+mesh%b2(1:3,i,j,k)*v%p(i,j,k)
                vx%p(i,j,k) = v_xyz(1)
                vy%p(i,j,k) = v_xyz(2)
                vz%p(i,j,k) = v_xyz(3)
            end do
        end do
    end do
end subroutine calc_cube2cart_vec_tile

subroutine calc_cube2cart_vec_tile_contra(vx, vy, vz, u, v, mesh)

    type(tile_field_t), intent(in)    :: u, v
    type(tile_field_t), intent(inout) :: vx, vy, vz
    type(tile_mesh_t),  intent(in)    :: mesh

    real(kind=8)    :: v_xyz(3)
    integer(kind=4) :: k, j, i

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                v_xyz(1:3) = mesh%a1(1:3,i,j,k)*u%p(i,j,k)+mesh%a2(1:3,i,j,k)*v%p(i,j,k)
                vx%p(i,j,k) = v_xyz(1)
                vy%p(i,j,k) = v_xyz(2)
                vz%p(i,j,k) = v_xyz(3)
            end do
        end do
    end do
end subroutine calc_cube2cart_vec_tile_contra

subroutine cart2cube_vec(u, v, vx, vy, vz, mesh, components_type)

    type(grid_field_t), intent(inout) :: u, v
    type(grid_field_t), intent(in)    :: vx, vy, vz
    type(mesh_t),       intent(in)    :: mesh
    character(len=*),   intent(in)    :: components_type

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        if(components_type == "covariant") then
            call calc_cart2cube_vec_tile(u%tile(t), v%tile(t), &
                                         vx%tile(t), vy%tile(t), vz%tile(t), &
                                         mesh%tile(t) )
        else if(components_type == "contravariant") then
            call calc_cart2cube_vec_tile_contra(u%tile(t), v%tile(t), &
                                         vx%tile(t), vy%tile(t), vz%tile(t), &
                                         mesh%tile(t) )
        else
            call parcomm_global%abort("cart2cube_vec, unknown components_type: "//components_type//", file: "//__FILE__)
        end if
    end do

end subroutine cart2cube_vec

subroutine calc_cart2cube_vec_tile(u, v, vx, vy, vz, mesh)

    type(tile_field_t), intent(inout) :: u, v
    type(tile_field_t), intent(in)    :: vx, vy, vz
    type(tile_mesh_t),  intent(in)    :: mesh

    integer(kind=4) :: k, j, i

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                u%p(i,j,k) = vx%p(i,j,k)*mesh%b1(1,i,j,k) + &
                             vy%p(i,j,k)*mesh%b1(2,i,j,k) + &
                             vz%p(i,j,k)*mesh%b1(3,i,j,k)

                v%p(i,j,k) = vx%p(i,j,k)*mesh%b2(1,i,j,k) + &
                             vy%p(i,j,k)*mesh%b2(2,i,j,k) + &
                             vz%p(i,j,k)*mesh%b2(3,i,j,k)
            end do
        end do
    end do
end subroutine calc_cart2cube_vec_tile

subroutine calc_cart2cube_vec_tile_contra(u, v, vx, vy, vz, mesh)

    type(tile_field_t), intent(inout) :: u, v
    type(tile_field_t), intent(in)    :: vx, vy, vz
    type(tile_mesh_t),  intent(in)    :: mesh

    integer(kind=4) :: k, j, i

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                u%p(i,j,k) = vx%p(i,j,k)*mesh%a1(1,i,j,k) + &
                             vy%p(i,j,k)*mesh%a1(2,i,j,k) + &
                             vz%p(i,j,k)*mesh%a1(3,i,j,k)

                v%p(i,j,k) = vx%p(i,j,k)*mesh%a2(1,i,j,k) + &
                             vy%p(i,j,k)*mesh%a2(2,i,j,k) + &
                             vz%p(i,j,k)*mesh%a2(3,i,j,k)
            end do
        end do
    end do
end subroutine calc_cart2cube_vec_tile_contra

subroutine multiply_by_J(Jf, f, mesh)

    type(grid_field_t), intent(in)    :: f
    type(grid_field_t), intent(inout) :: Jf
    type(mesh_t),       intent(in)    :: mesh

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        call multiply_by_J_tile(Jf%tile(t), f%tile(t), mesh%tile(t))
    end do

end subroutine multiply_by_J

subroutine multiply_by_J_tile(Jf, f, mesh)

    type(tile_field_t), intent(in)    :: f
    type(tile_field_t), intent(inout) :: Jf
    type(tile_mesh_t),  intent(in)    :: mesh

    integer(kind=4) :: i, j, k

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                Jf%p(i,j,k) = f%p(i,j,k)*mesh%J(i,j,k)
            end do
        end do
    end do

end subroutine multiply_by_J_tile

subroutine devide_by_J(Jf, f, mesh)

    type(grid_field_t), intent(in)    :: f
    type(grid_field_t), intent(inout) :: Jf
    type(mesh_t),       intent(in)    :: mesh

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        call devide_by_J_tile(Jf%tile(t), f%tile(t), mesh%tile(t))
    end do

end subroutine devide_by_J

subroutine devide_by_J_tile(Jf, f, mesh)

    type(tile_field_t), intent(in)    :: f
    type(tile_field_t), intent(inout) :: Jf
    type(tile_mesh_t),  intent(in)    :: mesh

    integer(kind=4) :: i, j, k

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                Jf%p(i,j,k) = f%p(i,j,k)/mesh%J(i,j,k)
            end do
        end do
    end do

end subroutine devide_by_J_tile
subroutine multiply_by_J_self(f, mesh)

    type(grid_field_t), intent(inout) :: f
    type(mesh_t),       intent(in)    :: mesh

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        call multiply_by_J_self_tile(f%tile(t), mesh%tile(t))
    end do

end subroutine multiply_by_J_self

subroutine multiply_by_J_self_tile(f, mesh)

    type(tile_field_t), intent(inout) :: f
    type(tile_mesh_t),  intent(in)    :: mesh

    integer(kind=4) :: i, j, k

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                f%p(i,j,k) = f%p(i,j,k)*mesh%J(i,j,k)
            end do
        end do
    end do

end subroutine multiply_by_J_self_tile

subroutine divide_by_J_self(f, mesh)

    type(grid_field_t), intent(inout) :: f
    type(mesh_t),       intent(in)    :: mesh

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        call divide_by_J_tile_self(f%tile(t), mesh%tile(t))
    end do

end subroutine divide_by_J_self

subroutine divide_by_J_tile_self(f, mesh)

    type(tile_field_t), intent(inout) :: f
    type(tile_mesh_t),  intent(in)    :: mesh

    integer(kind=4) :: i, j, k

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                f%p(i,j,k) = f%p(i,j,k)/mesh%J(i,j,k)
            end do
        end do
    end do

end subroutine divide_by_J_tile_self

subroutine add_random_noise(f,mesh,amplitude)
    type(grid_field_t), intent(inout)        :: f
    type(mesh_t),       intent(in)           :: mesh
    real(kind=8),       intent(in), optional :: amplitude

    integer(kind=4) :: t

    do t = mesh%ts,mesh%te
        call add_random_noise_tile(f%tile(t),mesh%tile(t),amplitude)
    end do
end subroutine

subroutine add_random_noise_tile(f,mesh,amplitude)
    type(tile_field_t), intent(inout)        :: f
    type(tile_mesh_t),  intent(in)           :: mesh
    real(kind=8),       intent(in), optional :: amplitude

    real(kind=8) :: amp_loc
    real(kind=8) :: noise(mesh%is:mesh%ie,mesh%js:mesh%je,mesh%ks:mesh%ke)

    amp_loc = 1.0_8
    if(present(amplitude)) amp_loc = amplitude

    call random_number(noise(mesh%is:mesh%ie,mesh%js:mesh%je,mesh%ks:mesh%ke))
    f%p(mesh%is:mesh%ie,mesh%js:mesh%je,mesh%ks:mesh%ke) = &
        f%p(mesh%is:mesh%ie,mesh%js:mesh%je,mesh%ks:mesh%ke) + &
            amp_loc*(noise(mesh%is:mesh%ie,mesh%js:mesh%je,mesh%ks:mesh%ke)-0.5_8)
end subroutine

subroutine multiply_by_Jvert_self(f,mesh)
    type(grid_field_t), intent(inout) :: f
    type(mesh_t),       intent(in)    :: mesh

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        call multiply_by_Jvert_self_tile(f%tile(t),mesh%tile(t))
    end do
end subroutine

subroutine multiply_by_Jvert_self_tile(f,mesh)
    type(tile_field_t), intent(inout) :: f
    type(tile_mesh_t),  intent(in)    :: mesh

    integer(kind=4) :: i, j, k

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                f%p(i,j,k) = f%p(i,j,k) * mesh%a3(4,i,j,k)
            end do
        end do
    end do
end subroutine

end module vec_math_mod
