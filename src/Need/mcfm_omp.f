!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      program mcfm
          use omp_lib
          use types
#ifdef HAVE_MPI
          use mpi
#endif
          use CPUTime
          use MCFMBenchmark
          use m_config, only : cfg_var_configadded
          use parseinput, only : cfg
          implicit none
          real(dp) :: r0,er0
          integer :: threadmax
#ifdef HAVE_MPI
          integer :: ierr,mylen,support
          character (len=mpi_max_processor_name) :: procname
#endif
          include 'mpicommon.f'

          integer:: exitCode

          threadmax = omp_get_max_threads()
          call omp_set_num_threads(threadmax)

#ifdef HAVE_MPI
          call mpi_init_thread(mpi_thread_funneled,support,ierr)
          call mpi_comm_rank(mpi_comm_world,rank,ierr)
          call mpi_comm_size(mpi_comm_world,world_size,ierr)
          call mpi_get_processor_name(procname,mylen,ierr)
#elif HAVE_COARRAY
          world_size = num_images()
          rank = this_image() - 1
#else
          world_size = 1
          rank = 0
#endif

          call start_times()

          call mcfmmain(r0,er0)

          if (cfg_var_configadded(cfg, "extra%benchmark")) then
              exitCode = comparisonCode(2)
          else
              exitCode = 0
          endif

#ifdef HAVE_MPI
          ! this is to have a consistent exit code for all mpi ranks
          call mpi_bcast(exitCode, 1, mpi_integer, 0, mpi_comm_world, ierr)
          call mpi_finalize(ierr)
#endif

          if (exitCode /= 0) then
              error stop 1
          endif

      end program
