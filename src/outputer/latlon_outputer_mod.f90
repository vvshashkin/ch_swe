module latlon_outputer_mod

use outputer_abstract_mod, only : outputer_t, outputer_vector_t
use grid_field_mod,        only : grid_field_t
use exchange_abstract_mod, only : exchange_t
use domain_mod,            only : domain_t
use tiles_mod,             only : tiles_t
use abstract_regrid_mod,   only : regrid_t, regrid_vec_t

implicit none

type, public, extends(outputer_t) :: latlon_outputer_t

    class(regrid_t),   allocatable :: regrid
    integer(kind=4)                :: Nlon, Nlat, Nz

contains
    procedure, public :: write => latlon_write
end type latlon_outputer_t

type, public, extends(outputer_vector_t) :: latlon_vec_outputer_t

    class(regrid_vec_t),   allocatable :: regrid
    integer(kind=4)                    :: Nlon, Nlat, Nz

contains
    procedure, public :: write => latlon_vec_write
end type latlon_vec_outputer_t

contains

subroutine latlon_write(this, f, domain, file_name, rec_num)

    class(latlon_outputer_t),         intent(inout) :: this
    type(grid_field_t),               intent(inout) :: f
    character(*),                     intent(in)    :: file_name
    type(domain_t),                   intent(in)    :: domain
    integer(kind=4),  optional,       intent(in)    :: rec_num

    integer(kind=4) :: i, j, k, t, ks, ke, panel_ind
    integer(kind=4) :: ts, te, js, je, is, ie

    real(kind=8), allocatable :: f_latlon(:,:,:)

    if(domain%parcomm%myid == 0) allocate(f_latlon(this%Nlon,this%Nlat,this%Nz))

    call this%regrid%do_regrid(f_latlon, f, domain)

    if(domain%parcomm%myid == 0) then

        open(newunit=this%out_stream, file = trim(file_name), &
             access="direct", recl = this%Nlat*this%Nlon*this%Nz)

        write(this%out_stream,rec = rec_num) real(f_latlon,4)

        close(this%out_stream)

    end if

end subroutine latlon_write

subroutine latlon_vec_write(this, u, v, domain, file_name_u, file_name_v, rec_num)

    class(latlon_vec_outputer_t), intent(inout) :: this
    type(grid_field_t),           intent(inout) :: u, v
    character(*),                 intent(in)    :: file_name_u, file_name_v
    type(domain_t),               intent(in)    :: domain
    integer(kind=4),              intent(in)    :: rec_num

    integer(kind=4) :: i, j, k, t, ks, ke, panel_ind
    integer(kind=4) :: ts, te, js, je, is, ie
    real(kind=8), allocatable :: u_latlon(:,:,:), v_latlon(:,:,:)

    if(domain%parcomm%myid == 0) then
        allocate(u_latlon(this%Nlon,this%Nlat,this%Nz))
        allocate(v_latlon(this%Nlon,this%Nlat,this%Nz))
    end if

    call this%regrid%do_regrid(u_latlon, v_latlon, u, v, domain)

    if(domain%parcomm%myid == 0) then

        open(newunit=this%out_stream_u, file = trim(file_name_u), &
             access="direct", recl = this%Nlat*this%Nlon*this%Nz)

        write(this%out_stream_u,rec = rec_num) real(u_latlon,4)

        close(this%out_stream_u)

        open(newunit=this%out_stream_v, file = trim(file_name_v), &
             access="direct", recl = this%Nlat*this%Nlon*this%Nz)

        write(this%out_stream_v,rec = rec_num) real(v_latlon,4)

        close(this%out_stream_v)

    end if

end subroutine latlon_vec_write

end module latlon_outputer_mod
