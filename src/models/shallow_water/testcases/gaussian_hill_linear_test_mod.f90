module gaussian_hill_linear_test_mod

use test_fields_mod,     only : scalar_field_generator_t, &
                                set_vector_test_field, set_scalar_test_field, &
                                gaussian_hill_scalar_field_generator_t
use domain_mod,          only : domain_t
use generic_config_mod,  only : generic_config_t
use stvec_flexible_mod,  only : stvec_flexible_t, stvec_t
use grid_field_mod,      only : grid_field_t
use abstract_scorer_mod, only : scorer_t

use const_mod,           only : pi, Earth_omega, Earth_radii, &
                                Earth_grav, Earth_sidereal_T

implicit none

real(kind=8), private, parameter :: phi0_default = 0.0_8, &
                                    lam0_default = pi,    &
                                    d0_default   = 0.25_8 


type, extends(scorer_t) :: gaussian_hill_scorer_t
    character(len=:), allocatable :: refsol_path
    contains
        procedure print_scores
end type

contains

subroutine setup_Gaussian_hill_linear_test(state, scorer, config, v_components_type, domain)

    class(generic_config_t), intent(inout) :: config
    character(len=*),        intent(in)    :: v_components_type
    type(domain_t),          intent(in)    :: domain

    class(stvec_t),               intent(inout) :: state
    class(scorer_t), allocatable, intent(out)   :: scorer

    real(kind=8) :: phi0, lam0, d0
    type(grid_field_t), pointer :: h, u, v
    type(gaussian_hill_scalar_field_generator_t) :: height_field

    call config%get(phi0, "phi0", default=phi0_default)
    call config%get(lam0, "lam0", default=lam0_default)
    call config%get(d0,   "d0",   default=d0_default)

    height_field = gaussian_hill_scalar_field_generator_t(phi0=phi0,lam0=lam0,d0=d0)

    select type(state)
        class is (stvec_flexible_t)
    
            call state%get_field(h,"h")
            call set_scalar_test_field(h, height_field, domain%mesh_p, 0)
    
            call state%get_field(u,"u")
            call state%get_field(v,"v")
            call u%assign(0.0_8, domain%mesh_u)
            call v%assign(0.0_8, domain%mesh_v)

    end select

    scorer = gaussian_hill_scorer_t()
    select type (scorer)
    type is (gaussian_hill_scorer_t)
        call config%get(scorer%refsol_path,"refsol_path",default="ref_solutions/gauss_linear/gauss_h_cub128_xy.dat")
    end select

end subroutine

subroutine print_scores(this, state, domain, time)

    use mpi
    use stvec_swm_mod, only : stvec_swm_t

    class(gaussian_hill_scorer_t), intent(inout) :: this
    class(stvec_t),                intent(in)    :: state
    type(domain_t),                intent(in)    :: domain
    real(kind=8),                  intent(in)    :: time

    integer(kind=4) :: i, j, t, inc, ix, iy, np, it, ierr, np_glob, panel_ind
    real(kind=8)    :: l2, linf, e, l2_glob, linf_glob

    integer(kind=4), parameter :: Nx_ref_sol = 192
    real(kind=8),    parameter :: ref_sol_dt = 3600.0_8

    real(kind=8) :: buffer(Nx_ref_sol+1,Nx_ref_sol+1,6)
    character(len=256) :: str_buff, fmt_str

    it = nint(time/ref_sol_dt)
    inc = Nx_ref_sol / domain%partition%Nh

    open(117,file=this%refsol_path,access="direct",recl=2*6*(Nx_ref_sol+1)**2)
    read(117,rec=it+1) buffer
    close(117)

    select type(state)
    type is (stvec_swm_t)

        l2 = 0.0_8
        linf = 0.0_8

        np = 0
        do t = domain%mesh_xy%ts, domain%mesh_xy%te
            panel_ind = domain%mesh_xy%tile(t)%panel_ind
            do j = domain%mesh_xy%tile(t)%js, domain%mesh_xy%tile(t)%je
                iy = (j-1)*inc+1
                do i = domain%mesh_xy%tile(t)%is, domain%mesh_xy%tile(t)%ie
                    ix = (i-1)*inc+1
                    e = state%h%tile(t)%p(i,j,1)-buffer(ix,iy,panel_ind)
                    linf = max(abs(e),linf)
                    l2 = l2+e**2
                    np = np+1
                end do
            end do
        end do
    end select

    call mpi_allreduce(l2,   l2_glob,   1, MPI_DOUBLE, MPI_SUM, domain%parcomm%comm_w, ierr)
    call mpi_allreduce(np,   np_glob,   1, MPI_INT,    MPI_SUM, domain%parcomm%comm_w, ierr)
    call mpi_allreduce(linf, linf_glob, 1, MPI_DOUBLE, MPI_MAX, domain%parcomm%comm_w, ierr)

    l2_glob = sqrt(l2_glob/np_glob)
    write(str_buff,*) "Errors, t=", time/3600.0_8, " hours, l2 = ", l2_glob, ", linf=", linf_glob

    if(domain%parcomm%myid == 0) then
        write(fmt_str,"(A,I8.8,A)") "(A",len(trim(str_buff)),")"
        print trim(fmt_str), str_buff
    end if

end subroutine

end module gaussian_hill_linear_test_mod