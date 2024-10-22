!Module for equiangular cubed sphere metric parameters
module metric_2d_ecs_mod

    use metric_2d_mod,             only : metric_2d_t
    use cubed_sphere_topology_mod, only : cubed_sphere_topology_t
    use parcomm_mod,               only : parcomm_global
    use mesh_mod,                  only : mesh_t

    implicit none

    type, extends(metric_2d_t) :: metric_2d_ecs_t
        class(cubed_sphere_topology_t), allocatable :: topology
        real(kind=8) :: rotation_matrix_panel(3,3,6)
    contains
        procedure :: calc_r  => calculate_ecs_r

        procedure :: calc_a1   => calculate_ecs_a1
        procedure :: calc_a2   => calculate_ecs_a2

        procedure :: calc_b1   => calculate_ecs_b1
        procedure :: calc_b2   => calculate_ecs_b2

        procedure :: calc_Q    => calculate_ecs_Q
        procedure :: calc_Qi   => calculate_ecs_Qi

        procedure :: calc_J   => calculate_ecs_J

        procedure :: calc_G   => calculate_ecs_G
        procedure :: calc_S   => calculate_ecs_S

        procedure :: transform_cart2native => transform_cartesian_to_native_ecs

        procedure :: set_curvilinear_mesh

        procedure :: transform_xyz_to_native
        procedure :: calc_cartesian_hor_wind

    end type metric_2d_ecs_t

    contains

    pure function calculate_ecs_r(this, panel_ind, alpha, beta) result(r)
        class(metric_2d_ecs_t), intent(in) :: this
        integer(kind=4),     intent(in) :: panel_ind
        real(kind=8),        intent(in) :: alpha, beta
        real(kind=8)                    :: r(3)

        r = ecs_point_r_proto(alpha, beta)
        r = ecs_proto2realface(this%topology, this%rotation_matrix, panel_ind, r)

    end function calculate_ecs_r

    pure function calculate_ecs_a1(this, panel_ind, alpha, beta) result(a1)
        class(metric_2d_ecs_t), intent(in) :: this
        integer(kind=4),     intent(in) :: panel_ind
        real(kind=8),        intent(in) :: alpha, beta
        real(kind=8)                    :: a1(3)

        a1 = ecs_a1_proto(alpha, beta)
        a1 = ecs_proto2realface(this%topology, this%rotation_matrix, panel_ind, a1)

    end function calculate_ecs_a1

    pure function calculate_ecs_a2(this, panel_ind, alpha, beta) result(a2)
        class(metric_2d_ecs_t), intent(in) :: this
        integer(kind=4),     intent(in) :: panel_ind
        real(kind=8),        intent(in) :: alpha, beta
        real(kind=8)                    :: a2(3)

        a2 = ecs_a2_proto(alpha, beta)
        a2 = ecs_proto2realface(this%topology, this%rotation_matrix, panel_ind, a2)

    end function calculate_ecs_a2

    pure function calculate_ecs_b1(this, panel_ind, alpha, beta) result(b1)
        class(metric_2d_ecs_t), intent(in) :: this
        integer(kind=4),     intent(in) :: panel_ind
        real(kind=8),        intent(in) :: alpha, beta
        real(kind=8)                    :: b1(3)

        b1 = ecs_b1_proto(alpha, beta)
        b1 = ecs_proto2realface(this%topology, this%rotation_matrix, panel_ind, b1)

    end function calculate_ecs_b1

    pure function calculate_ecs_b2(this, panel_ind, alpha, beta) result(b2)
        class(metric_2d_ecs_t), intent(in) :: this
        integer(kind=4),     intent(in) :: panel_ind
        real(kind=8),        intent(in) :: alpha, beta
        real(kind=8)                    :: b2(3)

        b2 = ecs_b2_proto(alpha, beta)
        b2 = ecs_proto2realface(this%topology, this%rotation_matrix, panel_ind, b2)

    end function calculate_ecs_b2

    pure function ecs_proto2realface(topology, rotation_matrix, panel_ind, r_proto) result(r)
        !transform prototype-face 3d Cartesian coords to real (spherical / elliptical) face
        class(cubed_sphere_topology_t), intent(in) :: topology
        real(kind=8),                   intent(in) :: rotation_matrix(3, 3)
        integer(kind=4),                intent(in) :: panel_ind
        real(kind=8),                   intent(in) :: r_proto(3)
        real(kind=8)                               :: r(3)

        real(kind=8) :: rtemp(3)


        rtemp(1) = topology%ex(1,panel_ind)*r_proto(1) + &
                   topology%ey(1,panel_ind)*r_proto(2) - &
                   topology%n (1,panel_ind)*r_proto(3)

        rtemp(2) = topology%ex(2,panel_ind)*r_proto(1) + &
                   topology%ey(2,panel_ind)*r_proto(2) - &
                   topology%n (2,panel_ind)*r_proto(3)

        rtemp(3) = topology%ex(3,panel_ind)*r_proto(1) + &
                   topology%ey(3,panel_ind)*r_proto(2) - &
                   topology%n (3,panel_ind)*r_proto(3)

        r(1) = sum(rotation_matrix(1:3,1)*rtemp(1:3))
        r(2) = sum(rotation_matrix(1:3,2)*rtemp(1:3))
        r(3) = sum(rotation_matrix(1:3,3)*rtemp(1:3))

    end function ecs_proto2realface

    pure function ecs_point_r_proto(alpha, beta) result(r)
        real(kind=8), intent(in) :: alpha, beta
        real(kind=8)             :: r(3)

        !grid point at proto cube face: (assumes x = tan(alpha), y = tan(beta), z = 1)
        !coordinates at prototype spherical face: vec(r)/||r||, r = (x, y, z)
        r(1) = tan(alpha) / sqrt(1._8+tan(alpha)**2+tan(beta)**2)
        r(2) = tan(beta) / sqrt(1._8+tan(alpha)**2+tan(beta)**2)
        r(3) = 1._8 / dsqrt(1._8+tan(alpha)**2+tan(beta)**2)

    end function ecs_point_r_proto

    pure function ecs_a1_proto(alpha, beta) result(a1)
        !dr / d alpha vector at prototype face for given alpha&beta
        real(kind=8), intent(in)  :: alpha, beta
        real(kind=8)              :: a1(3)

        real(kind=8) :: ta, tb

        ta = tan(alpha)
        tb = tan(beta)
        !d (xyz)^T/ d alpha
        a1(1) = (1d0+ta**2d0)*(1d0+tb**2d0)/(1d0+ta**2d0+tb**2d0)**(3d0/2d0)
        a1(2) = -ta*tb*(1d0+ta**2d0)/(1d0+ta**2d0+tb**2d0)**(3d0/2d0)
        a1(3) = -ta*(1d0+ta**2d0) / (1d0+ta**2d0+tb**2d0)**(3d0/2d0)

    end function ecs_a1_proto

    pure function ecs_a2_proto(alpha, beta) result(a2)
        !dr / d beta vector at prototype face for given alpha&beta
        real(kind=8)              :: a2(3)
        real(kind=8), intent(in)  :: alpha, beta

        real(kind=8) :: ta, tb

        ta = tan(alpha)
        tb = tan(beta)
        !d (xyz)^T/ d beta
        a2(1) = -ta*tb*(1d0+tb**2d0)/(1d0+ta**2d0+tb**2d0)**(3d0/2d0)
        a2(2) = (1d0+ta**2d0)*(1d0+tb**2d0)/(1d0+ta**2d0+tb**2d0)**(3d0/2d0)
        a2(3) = -tb*(1d0+tb**2d0) / (1d0+ta**2d0+tb**2d0)**(3d0/2d0)

    end function ecs_a2_proto

    pure function ecs_b1_proto(alpha, beta) result(b1)
        !contravariant alpha vector at prototype face for given alpha&beta
        real(kind=8), intent(in)  :: alpha, beta
        real(kind=8)              :: b1(3)

        real(kind=8) :: ta, tb, sigm

        ta = tan(alpha)
        tb = tan(beta)
        sigm = sqrt(1._8+ta**2+tb**2)
        !d (xyz)^T/ d alpha
        b1(1) = (1._8+tb**2-(ta*tb)**2/(1._8+ta**2)) / sigm
        b1(2) = 0._8
        b1(3) = -ta*(1._8+tb**2/(1+ta**2))/sigm

    end function ecs_b1_proto

    pure function ecs_b2_proto(alpha, beta) result(b2)
        !contravariant beta vector at prototype face for given alpha&beta
        real(kind=8), intent(in)  :: alpha, beta
        real(kind=8)              :: b2(3)

        real(kind=8) :: ta, tb, sigm

        ta = tan(alpha)
        tb = tan(beta)
        sigm = sqrt(1._8+ta**2+tb**2)
        !d (xyz)^T/ d beta
        b2(1) = 0._8
        b2(2) = ((1+ta**2)-(ta*tb)**2/(1+tb**2))/sigm
        b2(3) =-tb*(ta**2/(1+tb**2)+1._8)/sigm

    end function ecs_b2_proto

    pure function calculate_ecs_Q(this, panel_ind, alpha, beta) result(Q)
        !Calculate metric tensor
        !Q = |a1*a1  a1*a2| = |Q(1) Q(2)|
        !    |a1*a2  a2*a2|   |Q(2) Q(3)|
        class(metric_2d_ecs_t), intent(in) :: this
        integer(kind=4),     intent(in) :: panel_ind
        real(kind=8),        intent(in) :: alpha, beta
        real(kind=8)                    :: Q(3)

        real(kind=8) :: ta, tb, sigm

        ta = tan(alpha)
        tb = tan(beta)
        sigm = sqrt(1._8+ta**2+tb**2)

        Q(1) = (1._8+ta**2)**2*(1._8+tb**2)/sigm**4
        Q(2) =-(1._8+ta**2)*(1._8+tb**2)*ta*tb/sigm**4
        Q(3) = (1._8+ta**2)*(1._8+tb**2)**2/sigm**4

    end function calculate_ecs_Q

    pure function calculate_ecs_Qi(this,panel_ind, alpha, beta) result(Qi)
        !Calculate inverse metric tensor
        !QI = inv |a1*a1  a1*a2| = |QI(1) QI(2)|
        !         |a1*a2  a2*a2|   |QI(2) QI(3)|
        class(metric_2d_ecs_t), intent(in) :: this
        integer(kind=4),     intent(in) :: panel_ind
        real(kind=8),        intent(in) :: alpha,beta
        real(kind=8)                    :: Qi(3)
        !local
        real(kind=8) ta, tb, sigm

        ta = tan(alpha);   tb = tan(beta)
        sigm = sqrt(1._8+ta**2+tb**2)

        Qi(1) = sigm**2 / (1._8+ta**2)
        Qi(2) = sigm**2*ta*tb / ((1+ta**2)*(1+tb**2))
        Qi(3) = sigm**2 / (1._8+tb**2)

    end function calculate_ecs_Qi

    pure function calculate_ecs_G(this, panel_ind, alpha, beta) result(G)
        !Calculate christoffel symbols
        class(metric_2d_ecs_t), intent(in) :: this
        integer(kind=4),     intent(in) :: panel_ind
        real(kind=8),        intent(in) :: alpha,beta
        real(kind=8)                    :: G(2, 2, 2)

        real(kind=8) :: X, Y, d2

        X = tan(alpha)
        Y = tan(beta)
        d2 = 1._8+X**2+Y**2

        G(1,1,1) = 2.0_8*X*Y**2 / d2
        G(1,2,1) =-Y*(1+Y**2) / d2
        G(2,1,1) =-Y*(1+Y**2) / d2
        G(2,2,1) = 0.0_8

        G(1,1,2) = 0.0_8
        G(1,2,2) =-X*(1+X**2) / d2
        G(2,1,2) =-X*(1+X**2) / d2
        G(2,2,2) = 2.0_8*X**2*Y / d2

    end function calculate_ecs_G

    pure function calculate_ecs_S(this, panel_ind, alpha, beta) result(S)
        !Calculate christoffel symbols - 1st kind
        class(metric_2d_ecs_t), intent(in) :: this
        integer(kind=4),       intent(in) :: panel_ind
        real(kind=8),          intent(in) :: alpha,beta
        real(kind=8)                      :: S(3, 2)

        real(kind=8) :: X, Y, sigma

        X = tan(alpha)
        Y = tan(beta)
        sigma = 1._8+X**2+Y**2

        S(1,1) = 2.0_8*x*y**2*(1.0_8+x**2)**2*(1.0_8+y**2) / sigma**3
        S(2,1) =-0.5_8*y*(1.0_8+y**2)*(1.0_8+x**2)*(-x**4+3.0_8*x**2*y**2+y**2+1.0_8) / sigma**3
        S(3,1) = x*(1.0_8+x**2)*(1.0_8+y**2)**2*(-x**2+y**2-1.0_8) / sigma**3

        S(1,2) = y*(1.0_8+y**2)*(1.0_8+x**2)**2*(-y**2+x**2-1.0_8) / sigma**3
        S(2,2) =-0.5_8*x*(1.0_8+x**2)*(1.0_8+y**2)*(-y**4+3.0_8*y**2*x**2+x**2+1.0_8) / sigma**3
        S(3,2) = 2.0_8*y*x**2*(1.0_8+y**2)**2*(1.0_8+x**2) / sigma**3

    end function calculate_ecs_S

    pure function calculate_ecs_J(this, panel_ind, alpha, beta) result(J)
        !Compute sqrt of metric tensor determinant
        class(metric_2d_ecs_t), intent(in) :: this
        integer(kind=4),     intent(in) :: panel_ind
        real(kind=8),        intent(in) :: alpha, beta
        real(kind=8)                    :: J

        real(kind=8) :: ta, tb, sigm

        ta = tan(alpha)
        tb = tan(beta)
        sigm = sqrt(1._8+ta**2+tb**2)

        J = (1._8+ta**2)*(1._8+tb**2) / sigm**3

    end function calculate_ecs_J

    subroutine transform_cartesian_to_native_ecs(this, panel_ind, alpha, beta, r)
        class(metric_2d_ecs_t), intent(in)  :: this
        integer(kind=4),     intent(out) :: panel_ind
        real(kind=8),        intent(out) :: alpha, beta
        real(kind=8),        intent(in)  :: r(3)

        real(kind=8)    :: r_dot_n, r_dot_n_max, r_loc(3)
        integer(kind=8) :: ipanel

        !find panel with maximum projection of r onto cube outer normal
        !this will be the panel the point r belongs to
        panel_ind = 1
        r_dot_n_max = 0.0_8
        do ipanel = 1, 6
            !minus sign because of inner normal stored in topology
            r_dot_n =-sum(r(1:3)*this%topology%n(1:3,ipanel))
            if(r_dot_n > r_dot_n_max) then
                panel_ind = ipanel
                r_dot_n_max = r_dot_n
            end if
        end do

        !it is assumed that ||r|| == 1
        alpha = atan( sum(this%topology%ex(1:3,panel_ind)*r(1:3)) / r_dot_n_max)
        beta  = atan( sum(this%topology%ey(1:3,panel_ind)*r(1:3)) / r_dot_n_max)

    end subroutine transform_cartesian_to_native_ecs

    subroutine set_curvilinear_mesh(this, mesh)
        class(metric_2d_ecs_t), intent(in)    :: this
        type(mesh_t),        intent(inout) :: mesh

        integer(kind=4) :: t, i, j, k, ks, ke, panel_ind, halo_width
        real(kind=8)    :: alpha, beta, ta, tb, sigm, r(3), a1(3), a2(3), b1(3), b2(3)



        do t = mesh%ts, mesh%te
            halo_width = mesh%tile(t)%halo_width
            ks = mesh%tile(t)%ks
            ke = mesh%tile(t)%ke
            panel_ind = mesh%tile(t)%panel_ind
            do j = mesh%tile(t)%js-halo_width, mesh%tile(t)%je+halo_width
                beta = mesh%tile(t)%get_beta(j)
                tb = tan(beta)
                do i = mesh%tile(t)%is-halo_width, mesh%tile(t)%ie+halo_width
                    alpha = mesh%tile(t)%get_alpha(i)
                    ta = tan(alpha)
                    sigm = sqrt(1._8+ta**2+tb**2)

                    r(1) = ta/sigm
                    r(2) = tb/sigm
                    r(3) = 1._8/sigm

                    r(1:3) = ecs_proto2realface(this%topology, this%rotation_matrix, panel_ind, r)
                    mesh%tile(t)%rx(i,j,ks) = r(1)
                    mesh%tile(t)%ry(i,j,ks) = r(2)
                    mesh%tile(t)%rz(i,j,ks) = r(3)

                    mesh%tile(t)%h(i,j,ks) = 0.0_8

                    a1(1) = (1d0+ta**2d0)*(1d0+tb**2d0)/(1d0+ta**2d0+tb**2d0)**(3d0/2d0)
                    a1(2) = -ta*tb*(1d0+ta**2d0)/(1d0+ta**2d0+tb**2d0)**(3d0/2d0)
                    a1(3) = -ta*(1d0+ta**2d0) / (1d0+ta**2d0+tb**2d0)**(3d0/2d0)
                    mesh%tile(t)%a1(1:3,i,j,ks) = ecs_proto2realface(this%topology, this%rotation_matrix, panel_ind, a1)
                    mesh%tile(t)%a1(4,i,j,ks) = 0.0_8

                    a2(1) = -ta*tb*(1d0+tb**2d0)/(1d0+ta**2d0+tb**2d0)**(3d0/2d0)
                    a2(2) = (1d0+ta**2d0)*(1d0+tb**2d0)/(1d0+ta**2d0+tb**2d0)**(3d0/2d0)
                    a2(3) = -tb*(1d0+tb**2d0) / (1d0+ta**2d0+tb**2d0)**(3d0/2d0)
                    mesh%tile(t)%a2(1:3,i,j,ks) = ecs_proto2realface(this%topology, this%rotation_matrix, panel_ind, a2)
                    mesh%tile(t)%a2(4,i,j,ks) = 0.0_8

                    mesh%tile(t)%a3(1:4,i,j,ks) = [0.0_8, 0.0_8, 0.0_8, 1.0_8]

                    b1(1) = (1._8+tb**2-(ta*tb)**2/(1._8+ta**2)) / sigm
                    b1(2) = 0._8
                    b1(3) =-ta*(1._8+tb**2/(1+ta**2))/sigm
                    mesh%tile(t)%b1(1:3,i,j,ks) = ecs_proto2realface(this%topology, this%rotation_matrix, panel_ind, b1)

                    b2(1) = 0._8
                    b2(2) = ((1+ta**2)-(ta*tb)**2/(1+tb**2))/sigm
                    b2(3) =-tb*(ta**2/(1+tb**2)+1._8)/sigm
                    mesh%tile(t)%b2(1:3,i,j,ks) = ecs_proto2realface(this%topology, this%rotation_matrix, panel_ind, b2)

                    mesh%tile(t)%b3(1:4,i,j,ks) = [0.0_8, 0.0_8, 0.0_8, 1.0_8]

                    mesh%tile(t)%Q(1,i,j,ks) = (1._8+ta**2)**2*(1._8+tb**2)/sigm**4
                    mesh%tile(t)%Q(2,i,j,ks) =-(1._8+ta**2)*(1._8+tb**2)*ta*tb/sigm**4
                    mesh%tile(t)%Q(3,i,j,ks) = (1._8+ta**2)*(1._8+tb**2)**2/sigm**4

                    mesh%tile(t)%Qi(1,i,j,ks) = sigm**2 / (1._8+ta**2)
                    mesh%tile(t)%Qi(2,i,j,ks) = sigm**2*ta*tb / ((1+ta**2)*(1+tb**2))
                    mesh%tile(t)%Qi(3,i,j,ks) = sigm**2 / (1._8+tb**2)

                    mesh%tile(t)%J(i,j,ks) = (1._8+ta**2)*(1._8+tb**2) / sigm**3

                end do
            end do
        end do

    end subroutine set_curvilinear_mesh

    subroutine transform_xyz_to_native(this,alpha,beta, panel_ind,x,y,z,mesh)
        use grid_field_mod, only : grid_field_t, tile_field_t
        use mesh_mod,       only : mesh_t
    
        class(metric_2d_ecs_t), intent(in)    :: this
        type(grid_field_t),     intent(inout) :: alpha, beta, panel_ind
        type(grid_field_t),     intent(in)    :: x, y, z 
        type(mesh_t),           intent(in)    :: mesh
    
        real(kind=8)    :: r(3), r_dot_n(6), r_dot_n_max
        integer(kind=4) :: t, i, j, k, p_ind, ipanel
    
        do t = mesh%ts, mesh%te
            do k = mesh%tile(t)%ks, mesh%tile(t)%ke
                do j = mesh%tile(t)%js, mesh%tile(t)%je
                    do i = mesh%tile(t)%is, mesh%tile(t)%ie
                        r(1:3) = [x%tile(t)%p(i,j,k), y%tile(t)%p(i,j,k), z%tile(t)%p(i,j,k)]

                        do ipanel = 1, 6
                            !minus sign because of inner normal stored in topology
                            r_dot_n(ipanel) =-sum(r(1:3)*this%topology%n(1:3,ipanel))
                        end do
                        p_ind = maxloc(r_dot_n,1)
                        r_dot_n_max = r_dot_n(p_ind)
                        panel_ind%tile(t)%p(i,j,k) = 1.0_8*p_ind
                        
                        !it is assumed that ||r|| == 1
                        alpha%tile(t)%p(i,j,k) = atan( sum(this%topology%ex(1:3,p_ind)*r(1:3)) / r_dot_n_max)
                        beta%tile(t)%p(i,j,k)  = atan( sum(this%topology%ey(1:3,p_ind)*r(1:3)) / r_dot_n_max)
                    end do
                end do
            end do
        end do
    
    end subroutine

    subroutine calc_cartesian_hor_wind(this, vx, vy, vz, u, v, &
                                       panel_ind, alpha, beta, &
                                       mesh, native_components_type)

        use grid_field_mod, only : grid_field_t, tile_field_t
        use mesh_mod,       only : mesh_t

        class(metric_2d_ecs_t), intent(in)    :: this
        type(grid_field_t),     intent(inout) :: vx, vy, vz
        type(grid_field_t),     intent(in)    :: u, v, alpha, beta, panel_ind
        type(mesh_t),           intent(in)    :: mesh
        character(len=*),       intent(in)    :: native_components_type

        integer(kind=4) :: t

        select case (native_components_type)
        case ("contravariant")
            do t = mesh%ts, mesh%te
                call calc_cartesian_hor_wind_contra(vx%tile(t), vy%tile(t), vz%tile(t), &
                                                    u%tile(t), v%tile(t),               &
                                                    panel_ind%tile(t), alpha%tile(t),   &
                                                    beta%tile(t), mesh%tile(t), &
                                                    this%topology,this%rotation_matrix_panel)
            end do
        case ("covariant")
            call parcomm_global%abort("ecs_2d_metric, calc cartesian hor wind, covariant native components not implemented")
        case default
            call parcomm_global%abort("ecs_2d_metric, calc cartesian hor wind, unknown native components type: "//native_components_type)
        end select

    end subroutine

    subroutine calc_cartesian_hor_wind_contra(vx,vy,vz,u,v,panel_ind,alpha,beta,mesh,topology, rotation_matrix)
        use grid_field_mod, only : tile_field_t
        use mesh_mod,       only : tile_mesh_t

        type(tile_field_t),            intent(inout) :: vx, vy, vz
        type(tile_field_t),            intent(in)    :: u, v, alpha, beta, panel_ind
        type(tile_mesh_t),             intent(in)    :: mesh
        type(cubed_sphere_topology_t), intent(in)    :: topology
        real(kind=8),                  intent(in)    :: rotation_matrix(3,3,6)

        integer(kind=4) :: i, j, k, p_ind
        real(kind=8) :: v_tmp1(3), a1(3), a2(3)
        real(kind=8) :: ta, tb, sigma

        do k = mesh%ks, mesh%ke
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie

                    ta = tan(alpha%p(i,j,k))
                    tb = tan(beta%p(i,j,k))
                    sigma = (1d0+ta**2d0+tb**2d0)**(3d0/2d0)
            
                    a1(1) = (1d0+ta**2d0)*(1d0+tb**2d0) / sigma
                    a1(2) = -ta*tb*(1d0+ta**2d0) / sigma
                    a1(3) = -ta*(1d0+ta**2d0) / sigma

                    a2(1) = -ta*tb*(1d0+tb**2d0)/ sigma
                    a2(2) = (1d0+ta**2d0)*(1d0+tb**2d0)/ sigma
                    a2(3) = -tb*(1d0+tb**2d0) / sigma

                    v_tmp1(1:3) = u%p(i,j,k)*a1(1:3) + v%p(i,j,k)*a2(1:3)
                    
                    p_ind = int(panel_ind%p(i,j,k))
 
                    vx%p(i,j,k) = sum(rotation_matrix(1:3,1,p_ind)*v_tmp1(1:3))
                    vy%p(i,j,k) = sum(rotation_matrix(1:3,2,p_ind)*v_tmp1(1:3))
                    vz%p(i,j,k) = sum(rotation_matrix(1:3,3,p_ind)*v_tmp1(1:3))
                end do
            end do
        end do        

    end subroutine

end module metric_2d_ecs_mod
