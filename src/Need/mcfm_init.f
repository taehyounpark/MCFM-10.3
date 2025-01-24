!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine mcfm_init()
          use omp_lib
          use mod_qcdloop_c
          use avh_olo
          use parseinput
          use MCFMStorage
          use MCFMBenchmark, only : setupBenchmark
          use MCFMSetupPlots, only: setup_plots
          use MCFMPlotting, only : plots_allocate
          use iso_fortran_env
      implicit none
c***********************************************************************
c                                                                      *
c  This routine should initialize any necessary variables and          *
c  perform the usual print-outs                                        *
c                                                                      *
c***********************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'cutoff.f'
      include 'limits.f'
      include 'npart.f'
      include 'phasemin.f'
      include 'xmin.f'
      include 'facscale.f'
      include 'scale.f'
      include 'verbose.f'
      include 'phot_dip.f'
      include 'includect.f'
      include 'frag.f'
      include 'energy.f'
      include 'masses.f'
      include 'nflav.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'mpicommon.f'
      include 'nproc.f'
      include 'runname.f'
      include 'rtsmin.f'
      include 'lib/TensorReduction/Include/TRtensorcontrol.f' 
      include 'lib/TensorReduction/Include/TRmaxindex.f' 

      real(dp):: p1ext(4),p2ext(4),p(mxpart,4),val
      common/pext/p1ext,p2ext
      data p/mxpart*3._dp,mxpart*4._dp,mxpart*0._dp,mxpart*5._dp/

      integer :: nplotmax
      common/nplotmax/nplotmax

      character(len=:), allocatable :: cmd
      character(len=cfg_string_len) :: benchmark
      character(len=cfg_string_len) :: part

      logical :: newPlotting

!$omp threadprivate(/pext/)

c initialize qcdloop cache for all omp threads
!$omp parallel
      call qlcachesize(100)
!$omp end parallel
c initialize OneLoop
      call olo_onshell(1d-8)
      call olo_unit(-1,'error')

cInitialise data in commonblock /pvmaxindex/
      maxcindex=3
      maxdindex=4
      maxeindex=5

      if (rank >= 1) verbose=.false.

      ! TODO: add this to configuration file
c Initialize parameter settings
      call mdata

      epinv=1.e1_dp
      epinv2=1.e1_dp

      if (rank >= 1) verbose=.false.

      call default_config()

      call cfg_update_from_arguments(cfg)
      call cfg_sort(cfg)

      ! the following are already needed for the banner
      call cfg_get(cfg, "general%nproc", nproc)
      call parse_ewcorr()
      call cfg_get(cfg, "general%part", part)
      call parse_part(part)

      if (rank == 0) call banner

      if (rank == 0 .and. world_size > 1) then
          write (*,*) ""
#ifdef HAVE_MPI
          write (*,*) "Running MCFM with ", world_size, " MPI processes"
#elif HAVE_COARRAY
          write (*,*) "Running MCFM with ", world_size, " Coarray images"
#endif
          write (*,*) ""
      endif
      if (rank == 0) then
          write (*,*) ""
          write (*,*) "Running MCFM with ", omp_get_max_threads(), " OMP threads"
          write (*,*) ""
      endif

      if (rank == 0) then
          write (*,*) ""
          write (*,*) "MCFM compiled with ", compiler_version(), " using the options ",
     &          compiler_options()
          write (*,*) ""
      endif

      if (rank ==0) then
          allocate(character(len=1024) :: cmd)
          call get_command(cmd)
          write (*,*) ""
          write (*,*) "Running MCFM as ", trim(cmd)
          write (*,*) ""
          deallocate(cmd)
      endif


      if (cfg_var_configadded(cfg, "extra%benchmark")) then
          call cfg_get(cfg, "extra%benchmark", benchmark)
          write(6,*)
          write(6,*) '****************************************'
          write(6,*) '*      Running in benchmark mode       *'
          write(6,*) '* bench = ', benchmark
          write(6,*) '****************************************'
          write(6,*)
          call setupBenchmark()
      endif

      call read_config()
      !call cfg_write(cfg, "stdout")

      if (verbose .and. rank == 1) then
      write(6,*)
      write(6,*) '****************************************'
      write(6,*) '*     Cross section in femtobarns      *'
      write(6,*) '****************************************'
      write(6,*)
      endif

c Counter-terms for radiation in top decay should be included
      includect=.true.

c Set-up incoming beams and PS integration cut-offs
c--- Note: version 6.4 onwards, scale cutoff with c.o.m. energy
c--- Note: since Sep. 2018 no longer scale since it is dimensionless
c      cutoff=cutoff*(sqrts/2000._dp)**2
      rtsmin=min(rtsmin,sqrt(wsqmin+cutoff))
      rtsmin=min(rtsmin,sqrt(bbsqmin+cutoff))
      taumin=(rtsmin/sqrts)**2
      xmin=1.e-8_dp

      p1ext(4)=-half*sqrts
      p1ext(1)=0._dp
      p1ext(2)=0._dp
      p1ext(3)=-half*sqrts

      p2ext(4)=-half*sqrts
      p2ext(1)=0._dp
      p2ext(2)=0._dp
      p2ext(3)=+half*sqrts

c Set-up run name
      call setrunname(scale,facscale,cfg%config_directory)

c npart=9 is a dummy value, to ensure that all histograms are included
      npart=9
      val=1.e-15_dp
      nprocbelow = nproc

      call cfg_get(cfg,"histogram%newstyle",newPlotting)

      if (newPlotting .or. nproc == 1610 .or. nproc == 1650) then
          call setup_plots()
          call plots_allocate()
      else

      ! This first call has two purposes:
      ! 1. determine nplotmax for histogram allocation.
      ! 2. must be run in parallel so that each thread initializes
      !    in "tagbook" mode and takes runs fully in "tagplot"
      !    for subsequent calls in the integration routines
!$omp parallel
      call nplotter(p,val,val**2,1)
!$omp end parallel

!$omp parallel
      call initHistogramStorage(nplotmax)
!$omp end parallel
      call initMasterStorage(nplotmax)

      endif

      p(:,:) = 0._dp

c Initialize flag for photon fragmentation dipoles
      phot_dip(:)=.false.
      fragint_mode=.false.
c Initialize integer:: used in TensorReduction to zero
      TRtensorcontrol=0

      end

