!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module CPUTime
#ifdef HAVE_MPI
    use mpi
#endif
    use types

    public :: start_times, get_cputime, get_walltime

    private

    real(dp), save :: cputime_start
    integer, save :: walltime_start

    contains

    subroutine start_times
        implicit none

        call cpu_time(cputime_start)
        call system_clock(walltime_start)
    end subroutine

    ! returns cpu time in seconds
    function get_cputime()
        implicit none
        real(dp) :: get_cputime
        real(dp) :: cputime_end
#ifdef HAVE_MPI
        integer :: ierr
        include 'mpicommon.f'
        real(dp) :: retval
#elif HAVE_COARRAY
        real(dp) :: retval
#else
        real(dp) :: retval
#endif

        call cpu_time(cputime_end)

        retval = cputime_end - cputime_start
#ifdef HAVE_MPI
        if (world_size > 1) then
            call mpi_allreduce(mpi_in_place, retval, 1, mpi_double_precision, &
                mpi_sum, mpi_comm_world, ierr)
        endif
#elif HAVE_COARRAY
        call co_sum(retval, result_image=1)
#endif

        get_cputime = retval
    end function

    ! returns walltime in seconds
    function get_walltime()
        implicit none
        real(dp) :: get_walltime
        integer :: walltime_end, count_rate

        call system_clock(walltime_end, count_rate)

        get_walltime = real(walltime_end - walltime_start,dp)/real(count_rate,dp)

    end function


end module
