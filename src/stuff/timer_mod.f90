module timer_mod

use str_util_mod, only : to_lower

implicit none

private

integer(kind=4), parameter :: max_string_length   = 128
integer(kind=4), parameter :: max_timers_count    = 300

integer(kind=4) :: active_timers_count = 0

character(len=max_string_length) :: timer_name(max_timers_count)

real(kind=8),    dimension(max_timers_count) :: total_time
real(kind=8),    dimension(max_timers_count) :: min_time
real(kind=8),    dimension(max_timers_count) :: mean_time
real(kind=8),    dimension(max_timers_count) :: max_time
real(kind=8),    dimension(max_timers_count) :: start_time
integer(kind=4), dimension(max_timers_count) :: calls_count
logical,         dimension(max_timers_count) :: is_timer_active

integer(kind=4) :: comm_w

public :: init_timer_mod, start_timer, stop_timer, print_timer

contains

subroutine init_timer_mod(mpi_comm)

    integer(kind=4),  intent(in) :: mpi_comm
    is_timer_active(:) = .false.

    comm_w = mpi_comm

    call start_timer("__total__")

end subroutine init_timer_mod

subroutine add_new_timer(section_name)

    character(len=*), intent(in) :: section_name

    active_timers_count = active_timers_count + 1

    if (active_timers_count > max_timers_count) then
        print*, "Error in timer_mod! Run out of timers, increase max_timers_count"
        stop
    end if

    timer_name (active_timers_count) = to_lower(section_name)
    total_time (active_timers_count) = 0.0_8
    calls_count(active_timers_count) = 0

end subroutine add_new_timer

subroutine start_timer(section_name)

    use mpi, only : mpi_wtime

    character(len=*), intent(in)  :: section_name

    integer(kind=4) :: timer_idx

    timer_idx = get_timer_idx(section_name)

    if (timer_idx > active_timers_count) call add_new_timer(section_name)

    if (is_timer_active(timer_idx)) then
        print*, "Error in start_timer! Timer "//section_name// " is already started!"
        stop
    end if

    is_timer_active(timer_idx) = .true.
    calls_count(timer_idx) = calls_count(timer_idx) + 1
    start_time(timer_idx) = mpi_wtime()

end subroutine start_timer

subroutine stop_timer(section_name)

    use mpi, only : mpi_wtime

    character(len=*), intent(in)  :: section_name

    integer(kind=4)                  :: timer_idx
    real(kind=8)                     :: stop_time

    stop_time = mpi_wtime()

    timer_idx = get_timer_idx(section_name)

    if (timer_idx > active_timers_count) then
        print*, "Error in stop_timer!!! Timer "//section_name// " is not found!!!"
        stop
    end if

    if (.not. is_timer_active(timer_idx)) then
        print*, "Error in stop_timer!!! "//section_name// " timer is not started!"
        stop
    end if

    total_time(timer_idx) = total_time(timer_idx) + (stop_time - start_time(timer_idx))
    is_timer_active(timer_idx) = .false.

end subroutine stop_timer

subroutine print_timer(net_section_name)

    use mpi, only : mpi_comm_rank

    character(len=*), optional, intent(in) :: net_section_name

    real(kind=8)    :: pct_time
    integer(kind=4) :: net_timer_idx, timer_idx, myid, err

    call stop_timer("__total__")

    net_timer_idx = 1

    if (present(net_section_name)) then
        net_timer_idx = get_timer_idx(net_section_name)

        if (net_timer_idx > active_timers_count) then
            print*, "Error in print_timer!!! Net_timer is not found!!!"
            stop
        end if
    end if

    call calculate_timer_stats()

    call mpi_comm_rank(comm_w, myid, err)

    if (myid == 0) then
        ! Write out timer information in wiki formatted table
        write(*,'(A2,A32,7(A2,A21),A2)')                          &
                '||','=           Routine            =',          &
                '||','=   min time(s)     =',                     &
                '||','=   mean time(s)    =',                     &
                '||','=   max time(s)     =',                     &
                '||','=     No. calls     =',                     &
                '||','=       %time       =',                     &
                '||','= time per call(s)  =',                     &
                '||'

        do timer_idx = 1, active_timers_count
            pct_time = mean_time(timer_idx)/mean_time(net_timer_idx)*100.0_8
            write(*,                                                &
                ('(A2,A32,A2,3(f21.8,A2),i21,2(A2,f21.8),A2)'))   &
                '||', trim(timer_name(timer_idx)),                      &
                '||', min_time(timer_idx),                         &
                '||', mean_time(timer_idx),                        &
                '||', max_time(timer_idx),                         &
                '||', calls_count(timer_idx),                               &
                '||', pct_time,                                            &
                '||', mean_time(timer_idx)/calls_count(timer_idx),           &
                '||'
        end do
     end if

    call start_timer("__total__")

end subroutine print_timer

subroutine calculate_timer_stats()

    use mpi

    integer(kind=4) :: timer_idx, proc_count, err

    call mpi_comm_size(comm_w, proc_count, err)

    ! check all timers are closed
    do timer_idx = 1, active_timers_count
        if (is_timer_active(timer_idx)) then
            write(*, '(A,A,A)') 'Error! Timer for section ',trim(timer_name(timer_idx)), &
                    ' not closed.'
        stop
        end if
    end do

    do timer_idx = 1, active_timers_count
        call mpi_allreduce(total_time(timer_idx), min_time(timer_idx),  1, MPI_DOUBLE_PRECISION, mpi_min, comm_w, err)
        call mpi_allreduce(total_time(timer_idx), max_time(timer_idx),  1, MPI_DOUBLE_PRECISION, mpi_max, comm_w, err)
        call mpi_allreduce(total_time(timer_idx), mean_time(timer_idx), 1, MPI_DOUBLE_PRECISION, mpi_sum, comm_w, err)
        mean_time(timer_idx) = mean_time(timer_idx) / proc_count
    end do

end subroutine calculate_timer_stats

function get_timer_idx(section_name) result(timer_idx)

    character(len=*), intent(in)  :: section_name
    integer(kind=4) :: timer_idx

    character(len=len(section_name)) :: lowname

    lowname = to_lower(section_name)

    do timer_idx = 1, active_timers_count
        if ( trim(lowname) == trim(timer_name(timer_idx)) ) return
    end do

end function get_timer_idx

end module timer_mod
