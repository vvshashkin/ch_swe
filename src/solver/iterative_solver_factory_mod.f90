module iterative_solver_factory_mod

use iterative_solver_mod,  only : iterative_solver_t
use generic_config_mod,    only : generic_config_t
use domain_mod,            only : domain_t
use cg_solver_mod,         only : cg_solver_t
use bicgstab_solver_mod,   only : bicgstab_solver_t
! use multigrid_factory_mod, only : create_multigrid
use chebyshev_solver_mod,  only : chebyshev_solver_t
use richardson_solver_mod, only : richardson_solver_t

implicit none

contains

subroutine create_iterative_solver(solver, solver_config, domain)

    class(iterative_solver_t), allocatable, intent(out)    :: solver
    class(generic_config_t),                intent(inout)  :: solver_config
    type(domain_t),                         intent(in)     :: domain

    character(len = :), allocatable :: solver_name
    integer(kind=4) :: max_iter
    real(kind=8)    :: relative_tolerance
    logical         :: is_verbose

    call solver_config%get(solver_name, "solver_name")
    call solver_config%get(max_iter, "max_iter")
    call solver_config%get(relative_tolerance, "rel_tol")
    call solver_config%get(is_verbose, "verbose", default = .false.)

    select case (solver_name)
    case ("cg")
        allocate(solver, source = cg_solver_t(max_iter = max_iter,    &
                             relative_tolerance = relative_tolerance, &
                             is_verbose = is_verbose))
    case ("bicgstab")
        allocate(solver, source = bicgstab_solver_t(max_iter = max_iter, &
                             relative_tolerance = relative_tolerance,    &
                             is_verbose = is_verbose))
    case ("richardson")
        richardson_solver: block
            real(kind=8) :: omega
            call solver_config%get(omega, "richardson_omega")
            allocate(solver, source = richardson_solver_t(omega = omega, &
                             max_iter = max_iter, &
                             relative_tolerance = relative_tolerance,    &
                             is_verbose = is_verbose))
        end block richardson_solver
    case ("multigrid")
        call create_multigrid(solver, solver_config, domain)
        solver%max_iter = max_iter
        solver%relative_tolerance = relative_tolerance
        solver%is_verbose = is_verbose
    case ("chebyshev")
        parse_config: block
            real(kind=8) :: lambda_max, lambda_min
            call solver_config%get(lambda_min, "lambda_min")
            call solver_config%get(lambda_max, "lambda_max")

            allocate(solver, source = chebyshev_solver_t(max_iter = max_iter, &
                                 relative_tolerance = relative_tolerance,    &
                                 is_verbose = is_verbose, &
                                 lambda_min = lambda_min, &
                                 lambda_max = lambda_max) )

        end block parse_config
    case default
        call domain%parcomm%abort("unknown solver_name: "//solver_name)
    end select

end subroutine

subroutine create_multigrid(solver, config, domain)

    use domain_factory_mod,          only : create_domain
    use grid_field_based_vector_mod, only : grid_field_based_vector_t
    use abstract_quadrature_mod,     only : quadrature_t
    use quadrature_factory_mod,      only : create_quadrature
    use prolongation_factory_mod,    only : create_prolongation_operator
    use restriction_factory_mod,     only : create_restriction_operator
    use linear_operator_factory_mod, only : create_linear_operator
    use preconditioner_factory_mod,  only : create_preconditioner
    use multigrid_solver_mod,        only : multigrid_solver_t
    use domain_factory_mod,          only : create_domain

    class(iterative_solver_t), allocatable, intent(out)   :: solver
    type(domain_t),                         intent(in)    :: domain
    class(generic_config_t),                intent(inout) :: config !this is whole model config


    class(generic_config_t), allocatable :: mg_config ! specific part of solver config
    type(multigrid_solver_t), allocatable :: mg

    integer(kind=4) :: l

    call config%get(mg_config, "multigrid")

    allocate(mg)

    call mg_config%get(mg%n_levels, "n_levels")
    call mg_config%get(mg%n_pre_smooth, "n_pre_smooth")
    call mg_config%get(mg%n_post_smooth, "n_post_smooth")
    call mg_config%get(mg%n_coarse_smooth, "n_coarse_smooth")

    allocate(mg%mg_data(mg%n_levels))

    create_domains: block

        class(generic_config_t), allocatable :: domain_config

        domain_config = domain%config

        call create_domain(mg%mg_data(1)%domain, domain_config, domain%parcomm, domain%partition)
        do l = 2, mg%n_levels
            call create_coarse_domain(mg%mg_data(l)%domain, mg%mg_data(l-1)%domain)
            print*, mg%mg_data(l)%domain%partition%te_global
        end do
    end block create_domains

    create_vectors: block

        class(quadrature_t),             allocatable :: quadrature
        type(grid_field_based_vector_t), allocatable :: tmp
        integer(kind=4) :: hw_xy

        character(len=:), allocatable :: quadrature_name

        call mg_config%get(quadrature_name,  "quadrature_name")

        hw_xy = 8

        do l = 1, mg%n_levels

            call create_quadrature(quadrature, quadrature_name, mg%mg_data(l)%domain%mesh_o)

            allocate(tmp)
            call tmp%init(hw_xy, 0, mg%mg_data(l)%domain%mesh_o, quadrature)

            allocate(mg%mg_data(l)%x,      source = tmp)
            allocate(mg%mg_data(l)%b,      source = tmp)
            allocate(mg%mg_data(l)%r,      source = tmp)
            allocate(mg%mg_data(l)%r_prec, source = tmp)

            deallocate(tmp)
            deallocate(quadrature)

        end do

    end block create_vectors

    create_intergrid_operators: block
        character(len=:), allocatable :: restriction_name, prolongation_name

        call mg_config%get(restriction_name,  "restriction_name" )
        call mg_config%get(prolongation_name, "prolongation_name")

        do l = 1, mg%n_levels-1
            call create_restriction_operator(mg%mg_data(l)%restriction_op, &
                                             restriction_name,             &
                                             mg%mg_data(l)%domain)

            call create_prolongation_operator(mg%mg_data(l)%prolongation_op, &
                                              prolongation_name,             &
                                              mg%mg_data(l+1)%domain)
        end do
    end block create_intergrid_operators

    create_linear_operators: block
        class(generic_config_t), allocatable :: subconfig

        call config%get(subconfig, "linear_operator")

        do l = 1, mg%n_levels
            call create_linear_operator(mg%mg_data(l)%linear_operator, &
                                        mg%mg_data(l)%domain, subconfig)
        end do
    end block create_linear_operators

    create_smoother: block

        use geosci_config_mod, only : geosci_config_t

        character(len=:), allocatable :: smoother_precond_name, smoother_name, coarse_solver_name
        class(generic_config_t), allocatable :: smoother_cfg
        type(geosci_config_t) :: coarse_solver_cfg

        call mg_config%get(smoother_precond_name, "smoother_precond")
        call mg_config%get(coarse_solver_name, "coarse_solver_name")

        call mg_config%get(smoother_cfg, "smoother")
        call smoother_cfg%get(smoother_name, "smoother_name")
        call smoother_cfg%get(smoother_precond_name, "precond")

        call smoother_cfg%set("rel_tol", 1.0_8*10**(-16))
        call smoother_cfg%set("verbose", .true.)
        call smoother_cfg%set("max_iter", -1)
        call smoother_cfg%set("solver_name", smoother_name)

        do l = 1, mg%n_levels-1
            call create_iterative_solver(mg%mg_data(l)%smoother, smoother_cfg, mg%mg_data(l)%domain)
        end do

        do l = 1, mg%n_levels
            call create_preconditioner(                          &
                                  mg%mg_data(l)%smooth_precond,  &
                                  mg%mg_data(l)%linear_operator, &
                                  smoother_precond_name,         &
                                  mg%mg_data(l)%domain)
        end do

        call coarse_solver_cfg%set("rel_tol", 1.0_8*10.0_8**(-12))
        call coarse_solver_cfg%set("verbose", .true.)
        call coarse_solver_cfg%set("max_iter", mg%n_coarse_smooth)
        call coarse_solver_cfg%set("solver_name", coarse_solver_name)

        call create_iterative_solver(mg%mg_data(mg%n_levels)%coarse_solver, coarse_solver_cfg, mg%mg_data(mg%n_levels)%domain)

    end block create_smoother

    call move_alloc(mg, solver)

end subroutine create_multigrid

subroutine create_coarse_domain(domain_coarse, domain_fine)

    use domain_factory_mod, only : create_domain
    use generic_config_mod, only : generic_config_t
    use partition_mod,      only : partition_t

    type(domain_t), intent(out)   :: domain_coarse
    type(domain_t), intent(inout) :: domain_fine

    type(partition_t) :: partition_coarse
    class(generic_config_t), allocatable :: domain_config
    integer(kind=4) :: N

    allocate(domain_config, source = domain_fine%config)

    call create_coarse_partition(partition_coarse, domain_fine%partition, &
                                   domain_fine%horizontal_staggering)

    call domain_fine%config%get(N, "N")
    call domain_config%set("N", N / 2, overwrite = .true.)

    call create_domain(domain_coarse, domain_config, &
                       domain_fine%parcomm, partition_coarse)

end subroutine create_coarse_domain

subroutine create_coarse_partition(partition_coarse, partition_fine, staggering_type)

    use parcomm_mod,   only : parcomm_global
    use partition_mod, only : partition_t

    type(partition_t), intent(out) :: partition_coarse
    type(partition_t), intent(in)  :: partition_fine
    character(len=*),  intent(in)  :: staggering_type

    integer(kind=4) :: t, is, ie, js, je, ks, ke, Nz, Nz_inter, Nh, Nt

    if (mod(partition_fine%Nh, 2) /= 0) then
        call parcomm_global%abort("Error! Number of fine points is not even")
    end if

    partition_coarse%Nh = partition_fine%Nh / 2
    partition_coarse%Nz = partition_fine%Nz

    partition_coarse%num_panels = partition_fine%num_panels
    partition_coarse%num_tiles  = partition_fine%num_tiles
    partition_coarse%proc_map   = partition_fine%proc_map
    partition_coarse%panel_map  = partition_fine%panel_map
    partition_coarse%ts         = partition_fine%ts
    partition_coarse%te         = partition_fine%te
    partition_coarse%ts_global  = partition_fine%ts_global
    partition_coarse%te_global  = partition_fine%te_global
    partition_coarse%Nt         = partition_fine%Nt

    Nz = partition_coarse%Nz
    Nz_inter = partition_fine%tiles_z%Nz
    Nh = partition_coarse%Nh
    Nt = partition_coarse%Nt

    call partition_coarse%tiles_o%init(Nt, Nz, Nh, Nh)

    do t = 1, partition_fine%Nt

        is = (partition_fine%tiles_o%tile(t)%is + 1) / 2
        ie = (partition_fine%tiles_o%tile(t)%ie    ) / 2

        js = (partition_fine%tiles_o%tile(t)%js + 1) / 2
        je = (partition_fine%tiles_o%tile(t)%je    ) / 2

        ks = partition_fine%tiles_o%tile(t)%ks
        ke = partition_fine%tiles_o%tile(t)%ke

        if ((ie-is+1<2) .or. (je-js+1<2)) then
            call parcomm_global%abort("Too many levels of multigrid for this"// &
                                       " number of mpi ranks!")
        end if
        call partition_coarse%tiles_o%tile(t)%init(is, ie, js, je, ks, ke)

    end do

    call partition_coarse%tiles_x%init (Nt, Nz, Nh,   Nh+1, partition_coarse%tiles_o%tile)
    call partition_coarse%tiles_y%init (Nt, Nz, Nh+1, Nh,   partition_coarse%tiles_o%tile)
    call partition_coarse%tiles_xy%init(Nt, Nz, Nh+1, Nh+1, partition_coarse%tiles_o%tile)

    call partition_coarse%tiles_z%init  (Nt, Nz_inter, Nh,   Nh,   partition_coarse%tiles_o%tile)
    call partition_coarse%tiles_xz%init (Nt, Nz_inter, Nh,   Nh+1, partition_coarse%tiles_o%tile)
    call partition_coarse%tiles_yz%init (Nt, Nz_inter, Nh+1, Nh,   partition_coarse%tiles_o%tile)
    call partition_coarse%tiles_xyz%init(Nt, Nz_inter, Nh+1, Nh+1, partition_coarse%tiles_o%tile)

    do t = 1, Nt
        if (partition_coarse%tiles_x%tile(t)%ie == nh) partition_coarse%tiles_x%tile(t)%ie = nh+1
        if (partition_coarse%tiles_y%tile(t)%je == nh) partition_coarse%tiles_y%tile(t)%je = nh+1

        if (partition_coarse%tiles_xy%tile(t)%ie == nh) partition_coarse%tiles_xy%tile(t)%ie = nh+1
        if (partition_coarse%tiles_xy%tile(t)%je == nh) partition_coarse%tiles_xy%tile(t)%je = nh+1

        if (partition_coarse%tiles_xz%tile(t)%ie == nh) partition_coarse%tiles_xz%tile(t)%ie = nh+1
        if (partition_coarse%tiles_yz%tile(t)%je == nh) partition_coarse%tiles_yz%tile(t)%je = nh+1

        if (partition_coarse%tiles_xyz%tile(t)%ie == nh) partition_coarse%tiles_xyz%tile(t)%ie = nh+1
        if (partition_coarse%tiles_xyz%tile(t)%je == nh) partition_coarse%tiles_xyz%tile(t)%je = nh+1

        if(partition_coarse%tiles_z%tile(t)%ke  == Nz)  partition_coarse%tiles_z%tile(t)%ke  = Nz_inter
        if(partition_coarse%tiles_xz%tile(t)%ke == Nz) partition_coarse%tiles_xz%tile(t)%ke = Nz_inter
        if(partition_coarse%tiles_yz%tile(t)%ke == Nz) partition_coarse%tiles_yz%tile(t)%ke = Nz_inter
        if(partition_coarse%tiles_xyz%tile(t)%ke == Nz) partition_coarse%tiles_xyz%tile(t)%ke = Nz_inter
    end do

    if (staggering_type == 'A') then
        partition_coarse%tiles_u = partition_coarse%tiles_o
        partition_coarse%tiles_v = partition_coarse%tiles_o
        partition_coarse%tiles_p = partition_coarse%tiles_o
    else if (staggering_type == 'Ah') then
        partition_coarse%tiles_u = partition_coarse%tiles_xy
        partition_coarse%tiles_v = partition_coarse%tiles_xy
        partition_coarse%tiles_p = partition_coarse%tiles_xy
    else if (staggering_type == 'C') then
        partition_coarse%tiles_u = partition_coarse%tiles_x
        partition_coarse%tiles_v = partition_coarse%tiles_y
        partition_coarse%tiles_p = partition_coarse%tiles_o
    else if (staggering_type == 'Ch') then
        partition_coarse%tiles_u = partition_coarse%tiles_y
        partition_coarse%tiles_v = partition_coarse%tiles_x
        partition_coarse%tiles_p = partition_coarse%tiles_xy
    else
        call parcomm_global%abort('Unknown staggering_type in partition initialization!')
    end if


end subroutine create_coarse_partition

end module iterative_solver_factory_mod
