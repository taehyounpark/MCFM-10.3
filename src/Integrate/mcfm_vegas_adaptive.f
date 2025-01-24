!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine mcfm_vegas_adaptive(integ,integ_err)
      use types
      use iso_fortran_env
      use Integration
      use mod_sobseq
      use MCFMStorage
      use MCFMPrint
      use PDFerrors
      use Scalevar
      use MCFMPlotting, only: plots_allocate
      use SCET, only : doMultitaucut, scetreweight, includeTaucutgrid
      use parseinput
      use MCFMSettings
      use ResummationIntegration
      use SingletopPrint, only: finalizeStorageFixedOrder,
     &        printcross_singletop, writeAllHistogramsSingletop
      use CPUTime, only: get_walltime
      use SCET, only: useQT
      use ptveto, only: usept
      use ggHwilson, only: expansionorder
#ifdef HAVE_MPI
      use mpi
#endif
c This is to avoid compiler version dependence of mpi.mod
c#ifdef HAVE_MPI
c      implicit none
c      include 'mpif.h'
c#endif
      implicit none
      include 'ipsgen.f'
      include 'reset.f'
      include 'kpart.f'
      include 'kprocess.f'
      include 'vegas_common.f'
      include 'nproc.f'
      include 'taucut.f'
      include 'mpicommon.f'
      include 'frag.f'

      include 'vegas_damp.f'

      real(dp), intent(out) :: integ, integ_err

      real(dp) :: region(2*mxdim)
      real(dp) :: lowint, virtint, realint, scetint, fragint, qtint, ptint
      external lowint, virtint, realint, scetint, fragint, qtint, ptint

      ! these hold values for all maxParts integration parts and a maximum of
      ! maxIps ipsgen contributions
      real(dp) :: totalsig, totalsd
      ! this array holds modified sd values to give
      ! easier integrations slight preference, see below how it's used
      real(dp) :: sdMod(maxParts,maxIps)
      integer :: initcalls(maxParts,maxIps)
      logical :: enableIps(maxParts,maxIps)
      logical :: computePart(maxParts)
      logical :: doFirstCall(maxParts,maxIps)
      logical :: warmupComplete(maxParts,maxIps)

      integer(kind=8) :: strides
      real(dp) :: skiprnd

      ! this sets the vegas grid dampening parameter for all parts
      ! in warmup(1) and final(2) runs
      ! setting it to 0 will disable grid adjustments
      real(dp) :: vegasDampPart(maxParts,2)

      ! saves whether initcalls has been specified through configuration
      logical :: callsConfigAdded(maxParts)

      real(dp) :: xreal,xreal2
      common/xreal/xreal,xreal2

      integer :: stage
c      integer :: hist_readsnapshot, hist_writesnapshot
      integer :: nprocabove
      logical :: origCoeffonly

      ! number of iterations to compute in one batch
      integer :: iterBatchWarmup
      integer :: iterBatch, iterBatch1, iterBatch2
      ! multiplicator for number of calls if precision per batch is
      ! too low in warmup; also multiplier for full run
      real(dp) :: iterCallMult
      ! precision goal in percent of one batch for warmup
      real(dp) :: warmupPrecisionGoal
      ! precision goal of final result
      real(dp) :: resultPrecisionGoal
      real(dp) :: absResultPrecisionGoal
      real(dp) :: maxWalltime

      real(dp) :: globalCallMult

      ! factor to multiply number of calls after warmup for full run
      real(dp) :: callBoost

      logical :: usesobol
      logical :: writeIntermediate

      ! when this number of calls per iteration has been reached,
      ! switch iterBatch to 1
      integer :: maxCallsPerIter

      integer :: nextLoc(2)
      integer :: i,j,k
      integer :: part,ips
      logical :: warmdone

      integer :: nprocextra,ndimextra

      logical :: bin
      common/bin/bin

      logical :: dryrun
      common/dryrun/dryrun

      logical :: readin

      character(len=1024):: runname
      common/runname/runname

      integer(kind=8), parameter, dimension(1:30)   :: s = [
     &      1,2,3,3,4,
     &      4,5,5,5,5,
     &      5,5,6,6,6,
     &      6,6,6,7,7,
     &      7,7,7,7,7,
     &      7,7,7,7,7
     &      ]
      integer(kind=8), parameter, dimension(1:30)   :: a = [
     &      0,1,1,2,1,
     &      4,2,4,7,11,
     &      13,14,1,13,16,
     &      19,22,25,1,4,
     &      7,8,14,19,21,
     &      28,31,32,37,41
     &      ]
      integer(kind=8), parameter, dimension(7,1:30) :: m = reshape([
     &              1,0,0,0,0,0,0,
     &              1,3,0,0,0,0,0,
     &              1,3,1,0,0,0,0,
     &              1,1,1,0,0,0,0,
     &              1,1,3,3,0,0,0,
     &              1,3,5,13,0,0,0,
     &              1,1,5,5,17,0,0,
     &              1,1,5,5,5,0,0,
     &              1,1,7,11,19,0,0,
     &              1,1,5,1,1,0,0,
     &              1,1,1,3,11,0,0,
     &              1,3,5,5,31,0,0,
     &              1,3,3,9,7,49,0,
     &              1,1,1,15,21,21,0,
     &              1,3,1,13,27,49,0,
     &              1,1,1,15,7,5,0,
     &              1,3,1,15,13,25,0,
     &              1,1,5,5,19,61,0,
     &              1,3,7,11,23,15,103,
     &              1,3,7,13,13,15,69,
     &              1,1,3,13,7,35,63,
     &              1,3,5,9,1,25,53,
     &              1,3,1,13,9,35,107,
     &              1,3,1,5,27,61,31,
     &              1,1,5,11,19,41,61,
     &              1,3,5,3,3,13,69,
     &              1,1,7,13,1,19,1,
     &              1,3,7,5,13,19,59,
     &              1,1,3,9,25,29,41,
     &              1,3,5,13,23,1,55
     &      ], [7,30])
      integer(kind=8) :: toskip

      interface
          function fxn(vector,wgt)
              use types
              implicit none
              include 'mxdim.f'
              real(dp) :: fxn
              real(dp), intent(in) :: vector(mxdim)
              real(dp), intent(in) :: wgt
          end function
      end interface

      procedure (fxn), pointer :: ifun => null ()

      logical :: test

      call cfg_get(cfg, "integration%precisiongoal", resultPrecisionGoal)
      call cfg_get_add(cfg, "integration%absprecisiongoal", absResultPrecisionGoal, 0._dp, "")
      call cfg_get(cfg, "integration%warmupprecisiongoal", warmupPrecisionGoal)
          call cfg_get_add(cfg, "integration%maxwalltime", maxWalltime, 86400*365d0, "")
      call cfg_get(cfg, "integration%warmupchisqgoal", warmupChisqGoal)
      call cfg_get(cfg, "integration%usesobol", usesobol)
      call cfg_get(cfg, "integration%writeintermediate", writeIntermediate)

      call cfg_get_add(cfg, "integration%ndmx", ndmx, 100, "Number of Vegas grid subdivisions.")
      call cfg_get_add(cfg, "integration%callboost", callBoost, 4._dp, "")
      call cfg_get_add(cfg, "integration%itercallmult", iterCallMult, 1.4_dp,
     &      "Multiply calls/it by this factor after every iterBatch number of iterations.")

      call cfg_get_add(cfg, "integration%iterbatch1", iterBatch1, 10,
     &      "Batch of iterations using same calls/it just after warmup.")

      call cfg_get_add(cfg, "integration%iterbatch2", iterBatch2, 2,
     &      "Batch of iterations using same calls/it after first iterbatch1.")

      call cfg_get_add(cfg, "integration%iterbatchwarmup", iterBatchWarmup, 5,
     &      "Number of iterations for warmup using same calls/it.")
      call cfg_get_add(cfg, "integration%maxcallsperiter", maxCallsPerIter, 400000000,
     &      "Switch to iterBatch=1 after this number of calls/it.")

      call cfg_get_add(cfg, "integration%globalcallmult", globalCallMult,
     &      1.0_dp, "global call multiplier")

      call cfg_get(cfg, "integration%readin", readin)

      ! no binning for warumup stage
      bin = .false.

      ! we don't need or want this here
      dryrun = .false.
      nprn = 0

      do part=1,maxParts; do ips=1,maxIps
          iterationStorage(part,ips)%used = .false.
          associate (info => iterationStorage(part,ips)%vinfo)
              info%si = 0._dp
              info%swgt = 0._dp
              info%schi = 0._dp
              info%lastIter = 0
              info%doFirstCall = .false. ! ensure a first full run?
              info%warmupComplete = .false.
              info%damp = [1.5_dp, 0.8_dp]
              info%useSobol = usesobol
              info%sobolInitialized = .false.
              info%sobolSkip = 1000
          end associate
      enddo; enddo

c---      setup inital number of calls for warmup phase
c---      these could also be tuned to specific processes
c---
c---      note that the MC integration error is only accurate if the
c---      number of calls is "sufficiently" large, especially for
c---      difficult integrands

      callsConfigAdded(:) = .false.

      call cfg_get_add(cfg, "integration%initcallslord", initcalls(lord,1),
     &                  100000, "initial calls LO")
      initcalls(lord,2:maxIps) = initcalls(lord,1)
      if (cfg_var_configadded(cfg, "integration%initcallslord")) then
          callsConfigAdded(lord) = .true.
      endif

      call cfg_get_add(cfg, "integration%initcallsnlovirt", initcalls(nloVirt,1),
     &                  100000, "initial calls NLO virt")
      initcalls(nloVirt,2:maxIps) = initcalls(nloVirt,1)
      if (cfg_var_configadded(cfg, "integration%initcallsnlovirt")) then
          callsConfigAdded(nloVirt) = .true.
      endif

      call cfg_get_add(cfg, "integration%initcallsnloreal", initcalls(nloReal,1),
     &                  500000, "initial calls NLO real")
      initcalls(nloReal,2:maxIps) = initcalls(nloReal,1)
      if (cfg_var_configadded(cfg,"integration%initcallsnloreal")) then
          callsConfigAdded(nloReal) = .true.
      endif

      call cfg_get_add(cfg, "integration%initcallsnlofrag", initcalls(nloFrag,1),
     &                  100000, "initial calls NLO frag")
      initcalls(nloFrag,2:maxIps) = initcalls(nloFrag,1)
      if (cfg_var_configadded(cfg,"integration%initcallsnlofrag")) then
          callsConfigAdded(nloFrag) = .true.
      endif

      call cfg_get_add(cfg, "integration%initcallsnlorealextra", initcalls(nloRealExtra,1),
     &                  500000, "initial calls NLO real extra")
      initcalls(nloRealExtra,2:maxIps) = initcalls(nloRealExtra,1)
      if (cfg_var_configadded(cfg,"integration%initcallsnlorealextra")) then
          callsConfigAdded(nloRealExtra) = .true.
      endif

      call cfg_get_add(cfg, "integration%initcallssnlobelow", initcalls(snloBelow,1),
     &                  200000, "initial calls sNLO below")
      initcalls(snloBelow,2:maxIps) = initcalls(snloBelow,1)
      if (cfg_var_configadded(cfg, "integration%initcallssnlobelow")) then
          callsConfigAdded(snloBelow) = .true.
      endif

          call cfg_get_add(cfg, "integration%initcallsqtnlobelow", initcalls(qtnloBelow,1),
     &                              200000, "initial calls sNLO below")
          initcalls(qtnloBelow,2:maxIps) = initcalls(qtnloBelow,1)
          if (cfg_var_configadded(cfg, "integration%initcallsqtnlobelow")) then
              callsConfigAdded(qtnloBelow) = .true.
          endif

      call cfg_get_add(cfg, "integration%initcallssnloabove", initcalls(snloAbove,1),
     &                  400000, "initial calls sNLO above")
      initcalls(snloAbove,2:maxIps) = initcalls(snloAbove,1)
      if (cfg_var_configadded(cfg, "integration%initcallssnloabove")) then
          callsConfigAdded(snloAbove) = .true.
      endif

      call cfg_get_add(cfg, "integration%initcallsnnlobelow", initcalls(nnloBelow,1),
     &                  200000, "initial calls NNLO below")
      initcalls(nnloBelow,2:maxIps) = initcalls(nnloBelow,1)
      if (cfg_var_configadded(cfg, "integration%initcallsnnlobelow")) then
          callsConfigAdded(nnloBelow) = .true.
      endif

          call cfg_get_add(cfg, "integration%initcallsqtnnlobelow", initcalls(qtnnloBelow,1),
     &                              200000, "initial calls NNLO below")
          initcalls(qtnnloBelow,2:maxIps) = initcalls(qtnnloBelow,1)
          if (cfg_var_configadded(cfg, "integration%initcallsqtnnlobelow")) then
              callsConfigAdded(qtnnloBelow) = .true.
          endif

      call cfg_get_add(cfg, "integration%initcallsnnlovirtabove", initcalls(nnloVirtAbove,1),
     &                  1000000, "initial calls nnlo virt above")
      initcalls(nnloVirtAbove,2:maxIps) = initcalls(nnloVirtAbove,1)
      if (cfg_var_configadded(cfg, "integration%initcallsnnlovirtabove")) then
          callsConfigAdded(nnloVirtAbove) = .true.
      endif

      call cfg_get_add(cfg, "integration%initcallsnnlorealabove", initcalls(nnloRealAbove,1),
     &                  10000000, "initial calls nnlo real above")
      initcalls(nnloRealAbove,2:maxIps) = initcalls(nnloRealAbove,1)
      if (cfg_var_configadded(cfg, "integration%initcallsnnlorealabove")) then
          callsConfigAdded(nnloRealAbove) = .true.
      endif

      ! resummation

          call cfg_get_add(cfg, "integration%initcallsqtn3lobelow", initcalls(qtn3loBelow,1),
     &                              200000, "initial calls N3LO below")
          initcalls(qtn3loBelow,2:maxIps) = initcalls(qtn3loBelow,1)
          if (cfg_var_configadded(cfg, "integration%initcallsqtn3lobelow")) then
              callsConfigAdded(qtn3loBelow) = .true.
          endif
      call cfg_get_add(cfg, "integration%initcallsnloresummed", initcalls(nloResummed,1),
     &                  4000, "initial calls for qt below cut")
      initcalls(nloResummed,2:maxIps) = initcalls(nloResummed,1)
      if (cfg_var_configadded(cfg, "integration%initcallsnloresummed")) then
          callsConfigAdded(nloResummed) = .true.
      endif

      call cfg_get_add(cfg, "integration%initcallsnloresexp", initcalls(nloResexp,1),
     &                  10000, "initial calls for resummed f.o. expansion")
      initcalls(nloResexp,2:maxIps) = initcalls(nloResexp,1)
      if (cfg_var_configadded(cfg, "integration%initcallsnloresexp")) then
          callsConfigAdded(nloResexp) = .true.
      endif


      call cfg_get_add(cfg, "integration%initcallsnloresabove", initcalls(nloResAbove,1),
     &                  50000, "initial calls for nlo res above")
      initcalls(nloResAbove,2:maxIps) = initcalls(nloResAbove,1)
      if (cfg_var_configadded(cfg, "integration%initcallsnloresabove")) then
          callsConfigAdded(nloResAbove) = .true.
      endif

      call cfg_get_add(cfg, "integration%initcallsnnloresvirtabove", initcalls(nnloResVirtAbove,1),
     &                  1000000, "initial calls for nnlo res virt above")
      initcalls(nnloResVirtAbove,2:maxIps) = initcalls(nnloResVirtAbove,1)
      if (cfg_var_configadded(cfg, "integration%initcallsnnloresvirtabove")) then
          callsConfigAdded(nnloResVirtAbove) = .true.
      endif

      call cfg_get_add(cfg, "integration%initcallsnnloresrealabove", initcalls(nnloResRealAbove,1),
     &                  4000000, "initial calls for nnlo res real above")
      initcalls(nnloResRealAbove,2:maxIps) = initcalls(nnloResRealAbove,1)
      if (cfg_var_configadded(cfg, "integration%initcallsnnloresrealabove")) then
          callsConfigAdded(nnloResRealAbove) = .true.
      endif

! pt-veto resummation

      call cfg_get_add(cfg, "integration%initcallsnloresvetovirt", initcalls(nloresvetovirt,1),
     &                  50000, "initial calls for nlo veto res above")
      initcalls(nloresvetovirt,2:maxIps) = initcalls(nloresvetovirt,1)
      if (cfg_var_configadded(cfg, "integration%initcallsnloresvetovirt")) then
          callsConfigAdded(nloresvetovirt) = .true.
      endif

      call cfg_get_add(cfg, "integration%initcallsnloresvetoreal", initcalls(nloresvetoreal,1),
     &                  50000, "initial calls for nlo veto res above")
      initcalls(nloresvetoreal,2:maxIps) = initcalls(nloresvetoreal,1)
      if (cfg_var_configadded(cfg, "integration%initcallsnloresvetoreal")) then
          callsConfigAdded(nloresvetoreal) = .true.
      endif

      call cfg_get_add(cfg, "integration%initcallsnnloresvetobelow", initcalls(nnloresvetobelow,1),
     &                  100000, "initial calls for nnlo veto res below")
      initcalls(nnloresvetobelow,2:maxIps) = initcalls(nnloresvetobelow,1)
      if (cfg_var_configadded(cfg, "integration%initcallsnnloresvetobelow")) then
          callsConfigAdded(nnloresvetobelow) = .true.
      endif

      call cfg_get_add(cfg, "integration%initcallsnnloresvetorealabove", initcalls(nnloresvetorealabove,1),
     &                  4000000, "initial calls for nnlo veto res real above")
      initcalls(nnloresvetorealabove,2:maxIps) = initcalls(nnloresvetorealabove,1)
      if (cfg_var_configadded(cfg, "integration%initcallsnnloresvetorealabove")) then
          callsConfigAdded(nnloresvetorealabove) = .true.
      endif

      call cfg_get_add(cfg, "integration%initcallsnnloresvetovirtabove", initcalls(nnloresvetovirtabove,1),
     &                  1000000, "initial calls for nnlo veto res virt above")
      initcalls(nnloresvetovirtabove,2:maxIps) = initcalls(nnloresvetovirtabove,1)
      if (cfg_var_configadded(cfg, "integration%initcallsnnloresvetovirtabove")) then
          callsConfigAdded(nnloresvetovirtabove) = .true.
      endif
      ipsgen = 1

      origKpart = kpart
      origCoeffonly = coeffonly

c--- =====================================
c--- process specific initializations here
c--- =====================================

#define SETUP(part,num) if (.not. callsConfigAdded(part)) initcalls(part,:) = num

      enableIps(:,:) = .false.
      enableIps(:,1) = .true.

      if (nproc == 1 .or. nproc == 6) then
          SETUP(lord,5d4)
          SETUP(nloReal,1d6)
          SETUP(nloVirt,2d5)
          SETUP(nnloBelow,1d6)
          SETUP(nnloVirtAbove,4d6)
          SETUP(nnloRealAbove,4d7)
      endif

      if (nproc == 2 .or. nproc == 7) then
          SETUP(lord,5d4)
          SETUP(nloReal,1d6)
          SETUP(nloVirt,2d5)
          enableIps(nloReal,2) = .true.
      endif

      if (nproc >= 61 .and. nproc <= 90) then
        if ((nproc == 61) .or. (nproc == 81) .or. (nproc == 82)) then
 ! Additional region for gg->VV
          enableIps(nnloBelow,2) = .true.
          if (kresorder == 6) enableIps(nloResummed,2)=.true.
        endif
        SETUP(lord,4d5)
        SETUP(nloReal,4d6)
        SETUP(nloVirt,1d6)
        SETUP(nnloBelow,4d5)
        SETUP(nnloVirtAbove,1d6)
        SETUP(nnloRealAbove,2d7)
        SETUP(nloResummed,4d5)
        SETUP(nloResVetoReal,2d6)
        SETUP(nloResVetoVirt,4d5)
        SETUP(nnloResVetoBelow,1d6)
        SETUP(nnloResVetoRealAbove,1d7)
        SETUP(nnloResVetoVirtAbove,4d6)
      endif

      if (nproc >= 461 .and. nproc <= 490) then
        SETUP(lord,4d5)
        SETUP(nloReal,4d6)
        SETUP(nloVirt,1d6)
      endif

      if (nproc == 126) then
          enableIps(:,2) = .true.
        SETUP(lord,4d5)
      endif

      if (nproc == 8211) then
          enableIps(:,2) = .true.
      endif

      if (nproc == 290 .or. nproc == 295) then
        enableIps(:,2) = .true.
        SETUP(lord,2d5)
        SETUP(nloReal,2d6)
        SETUP(nloVirt,4d5)
        SETUP(nnloVirtAbove,1d7)
        SETUP(nnloRealAbove,2d7)
      endif

      if (nproc == 292 .or. nproc == 297) then
        enableIps(:,2) = .true.
        SETUP(lord,4d5)
        SETUP(nloReal,4d6)
        SETUP(nloVirt,2d6)
      endif

      if (nproc == 291 .or. nproc == 296) then
        enableIps(:,2) = .true.
        enableIps(nloReal,3:4) = .true.
        SETUP(nloReal,2d6) ! d6
        SETUP(nloVirt,1d5) ! d5
      endif

      if (nproc == 2921 .or. nproc == 2971) then
        enableIps(:,2) = .true.
        enableIps(nloReal,3:4) = .true.
        SETUP(nloReal,2d7) ! d7
        SETUP(nloVirt,1d6) ! d6
      endif

      if (nproc == 294 .or. nproc == 299) then
        enableIps(:,2) = .true.
        SETUP(nloReal,1d7)
        SETUP(nloVirt,1d6)
      endif

      if (nproc == 2941 .or. nproc == 2991) then
        enableIps(:,2) = .true.
        SETUP(nloReal,4d6) ! d7
        SETUP(nloVirt,2d6) ! d6
      endif

       if (nproc == 165) then
           SETUP(lord,2000000)
           SETUP(nloVirt,2000000)
           SETUP(nloReal,10000000)
       endif

      if (nproc == 300 .or. nproc == 3000) then
        if (origKpart == kresummed) then
            enableIps(:,1:2) = .true.
        else
            enableIps(:,1:2) = .true.
        endif
        SETUP(nloReal,1d6)
        SETUP(nloVirt,2d5)
      endif

      if (nproc == 302 .or. nproc == 3002) then
        if (origKpart == kresummed) then
            continue
        else
            enableIps([nloVirt,nloReal,lord],1:2) = .true.
        endif
      endif

      if (nproc == 304 .or. nproc == 309 .or. nproc == 3004 .or. nproc == 3009) then
        enableIps(lord,1:2) = .true.
        SETUP(lord,1d7)
      endif

      if (nproc == 164 .or. nproc==169) then
        iterationStorage(lord,1)%vinfo%damp = [1.5, 0.8]
        SETUP(lord,250000)

        iterationStorage(nloVirt,1)%vinfo%damp = [1.5, 0.8]
        SETUP(nloVirt,100000)

        enableIps(nloReal, 1:2) = .true.
        iterationStorage(nloReal,1)%vinfo%damp = [1.5, 0.8]
        SETUP(nloReal,500000)
      endif

      if (nproc == 1281 .or. nproc == 1291 .or. nproc == 1301 .or.
     &        nproc == 1311 .or. nproc == 1321 .or. nproc == 1282 .or.
     &        nproc == 1292 .or. nproc == 1302 .or. nproc == 1302 .or.
     &        nproc == 1312 .or. nproc == 1322) then
          enableIps(:,1:2) = .true.
      endif

      if (nproc == 226 .or. nproc == 2261) then
          enableIps(:,1:2) = .true.
      endif

      if ((nproc == 273) .or. (nproc == 274)) then
          SETUP(lord,200000)
          SETUP(nloVirt,400000)
          SETUP(nloReal,2000000)
      endif

      if ((nproc == 278) .or. (nproc == 279)) then
          SETUP(lord,800000)
      endif

      if ((nproc == 544) .or. (nproc == 547)) then
          SETUP(lord,1000000)
          SETUP(nloVirt,1000000)
          SETUP(nloReal,4000000)
      endif

      if ((nproc >= 609) .and. (nproc <= 623)) then
          SETUP(lord,200000)
          SETUP(nloVirt,200000)
          SETUP(nloReal,2000000)
      endif

      if ((nproc == 613) .or. (nproc == 618) .or. (nproc == 623)) then
          SETUP(lord,4000000)
          SETUP(nloVirt,8000000)
          SETUP(nloReal,20000000)
      endif

      if (nproc >= 640 .and. nproc <= 659) then
          SETUP(lord,1000000)
      endif

      if (nproc >= 661 .and. nproc <= 669) then
          SETUP(lord,5000000)
      endif

      if (nproc > 900 .and. nproc < 920) then
          SETUP(lord,1d6)
      endif
c---
c--- setup calls/it
c---
      do j=1,maxParts
          do k=1,maxIPS
            iterationStorage(j,:)%vinfo%callsPerIt = initcalls(j,k)
          enddo
      enddo

c---
c--- determine parts to calculate
c---
      computePart(:) = .false.

      if (origKpart == kresummed .and. kresorder == 2) then
          computePart(nloResummed) = .true.
      endif

      if (origKpart == kresummed .and. kresorder == 3) then
          computePart(nloResummed) = .true.
          computePart(nloResexp) = .true.
          computePart(nloResAbove) = .true.
      endif

      if (origKpart == kresummed .and. kresorder == 4) then
          computePart(nloResummed) = .true.
          computePart(nloResexp) = .true.
          computePart(nloResAbove) = .true.
      endif

      if (origKpart == kresummed .and. kresorder == 5) then
          ! order = 2 is handled inside these routines
          computePart(nloResummed) = .true.
          computePart(nloResexp) = .true.

          !computePart(nloResAbove) = .true.
          computePart(nnloResVirtAbove) = .true.
          computePart(nnloResRealAbove) = .true.
      endif

      if (origKpart == kresummed .and. kresorder == 6) then
          ! order = 2 is handled inside these routines
          computePart(nloResummed) = .true.
          computePart(nloResexp) = .true.

          !computePart(nloResAbove) = .true.
          computePart(nnloResVirtAbove) = .true.
          computePart(nnloResRealAbove) = .true.
      endif

      if (origKpart == kresummed .and. kresorder == 7) then
          ! order = 2 is handled inside these routines
          computePart(nloResummed) = .true.
          computePart(nloResexp) = .true.

          !computePart(nloResAbove) = .true.
          computePart(nnloResVirtAbove) = .true.
          computePart(nnloResRealAbove) = .true.
      endif

      if (krespart == kresexp) then
          computePart(:) = .false.
          computePart(nloResexp) = .true.
      elseif (krespart == kresonly) then
          computePart(:) = .false.
          computePart(nloResummed) = .true.
      elseif (krespart == kresmatchcorr .and. kresorder == 4 .and. usept) then
          computePart(:) = .false.
          computePart(nloResexp) = .true.
          computePart(nloResVetoVirt) = .true.
          computePart(nloResVetoReal) = .true.
      elseif (krespart == kresmatchcorr .and. kresorder == 4) then
          computePart(:) = .false.
          computePart(nloResexp) = .true.
          computePart(nloResAbove) = .true.
      elseif (krespart == kresmatchcorr .and. kresorder == 6 .and. usept) then
          computePart(:) = .false.
          computePart(nloResexp) = .true.
          computePart(nnloResVetoBelow) = .true.
          computePart(nnloResVetoRealAbove) = .true.
          computePart(nnloResVetoVirtAbove) = .true.
      elseif (krespart == kresmatchcorr .and. kresorder == 6) then
          computePart(:) = .false.
          computePart(nloResexp) = .true.
          computePart(nnloResVirtAbove) = .true.
          computePart(nnloResRealAbove) = .true.
      elseif (krespart == kresabove .and. kresorder == 4) then
          computePart(:) = .false.
          computePart(nloResAbove) = .true.
      elseif (krespart == kresabove .and. kresorder == 6) then
          computePart(:) = .false.
          computePart(nnloResVirtAbove) = .true.
          computePart(nnloResRealAbove) = .true.
      endif

      if (origKpart == klord) then
          computePart(lord) = .true.
      endif

      if (origKpart == kvirt) then
c          origCoeffonly = .true.
          computePart(nloVirt) = .true.
      endif

      if (origKpart == kreal) then
c          origCoeffonly = .true.
          computePart(nloReal) = .true.
      endif

      if (origKpart == kfrag) then
          origCoeffonly = .true.
          computePart(nloFrag) = .true.
      endif

      if (origKpart == ktota) then
          if (nproc == 1610 .and. (origCoeffonly .eqv. .false.)) then
              ! for this process the lord piece is calculated separately
              ! and nloVirt with coeffonly
              computePart(lord) = .true.
          endif
          computePart(nloVirt) = .true.
          computePart(nloReal) = .true.

          if (frag) then
              computePart(nloFrag) = .true.
          endif
      endif

      if (origKpart == ktodk) then
          computePart(nloVirt) = .true.
          computePart(nloReal) = .true.
          computePart(nloRealExtra) = .true.
          call setuprealextra(nprocextra)
      endif


      if (origKpart ==  ksnlo) then
          if (nproc == 1610 .and. (origCoeffonly .eqv. .false.)) then
              ! for this process the lord piece is calculated separately
              ! and nloVirt with coeffonly
              computePart(lord) = .true.
          endif
              if (useQT) then
                  computePart(qtnloBelow) = .true.
              else
                  computePart(snloBelow) = .true.
              endif
          computePart(snloAbove) = .true.
          if (onlypowcorr) computePart(snloAbove) = .false.

          if (ksnlopart == ksnloV) then
              computePart(:) = .false.
                  if (useQT) then
                      computePart(qtnloBelow) = .true.
                  else
                      computePart(snloBelow) = .true.
                  endif
              endif

          if (ksnlopart == ksnloR) then
              computePart(:) = .false.
              computePart(snloAbove) = .true.
          endif
      endif

      if (origKpart == knnlo) then
        if (origCoeffonly .eqv. .false.) then
          ! alternative enable snloBelow and snloAbove
          computePart(nloVirt) = .true.
          computePart(nloReal) = .true.

          if (nproc == 1610) then
              ! we set coeffonly to true for 161 below
              ! so include this
              computePart(lord) = .true.
          endif
        endif

        if (useQT) then
            computePart(qtnnloBelow) = .true.
        else
            computePart(nnloBelow) = .true.
        endif
        computePart(nnloVirtAbove) = .true.
        computePart(nnloRealAbove) = .true.
        if (onlypowcorr) then
          computePart(nnloVirtAbove) = .false.
          computePart(nnloRealAbove) = .false.
        endif
      endif

      if (knnlopart == knnloRR) then
        computePart(:) = .false.
        computePart(nnloRealAbove) = .true.
      endif

      if (knnlopart == knnloRV) then
        computePart(:) = .false.
        computePart(nnloVirtAbove) = .true.
      endif

      if (knnlopart == knnloVV) then
        computePart(:) = .false.
        origcoeffonly = .true.
        if (useQT) then
            computePart(qtnnloBelow) = .true.
        else
            computePart(nnloBelow) = .true.
        endif
      endif

      if (origkpart == kn3lo) then
        ! for now we just compute the below cut part
        computePart(:) = .false.
        origcoeffonly = .true.
        computePart(qtn3loBelow) = .true.
      endif

      if (any(origKpart == [ksnlo,knnlo,kn3lo])) then
        call setupscet(nprocabove)
      elseif (origKpart == kResummed) then
          call setupscet(nprocabove)
          !if (nproc == 31) then
              !nprocabove = 41
          !else
              !error stop __FILE__//": this process is not (yet) set up for resummation"
          !endif
      endif

          ! special case for Z+jet NNLO
      if (nprocbelow == 41) then

          ! when splitting real PS into three pieces
          !enableIps(nnloRealAbove,1:3) = .true.
          !initcalls(nnloRealAbove,1:3) = 1000000

          enableIps(nnloBelow,1) = .true.
          initcalls(nnloBelow,1) = 5000000

          ! use just one real emission PS
          enableIps(nnloRealAbove,1) = .true.
          initcalls(nnloRealAbove,1) = 20000000

          enableIps(nnloVirtAbove,1) = .true.
          initcalls(nnloVirtAbove,1) = 8000000
      endif

      if (nproc == 44) then
          enableIps(nloReal,1:1) = .true.
          initcalls(nloReal,1:1) = 2000000
      endif

      if ((nprocbelow == 204) .or. (nprocbelow == 210))then
          enableIps(nnloRealAbove,1) = .true.
          initcalls(nnloRealAbove,1) = 10000000

          enableIps(nnloVirtAbove,1) = .true.
          initcalls(nnloVirtAbove,1) = 2000000
      endif

      ! special case for nnlo singletop
      if (nproc == 1650) then
          enableIps(:,:) = .false.

          call cfg_get(cfg, "singletop%nnlo_enable_light", test)
          if (test) then
              enableIps(nloVirt,1) = .true.
              initcalls(nloVirt,:) = 2000000
              enableIps(nloReal,1:4) = .true.
              initcalls(nloReal,:) = 8000000
          endif

          call cfg_get(cfg, "singletop%nnlo_enable_heavy_prod", test)
          if (test) then
              enableIps(lord,5:6) = .true.
              initcalls(lord,5:6) = 100000

              enableIps(nloVirt,2) = .true.
              initcalls(nloVirt,:) = 2000000
              enableIps(nloReal,5:8) = .true.
              initcalls(nloReal,:) = 8000000
          endif

          call cfg_get(cfg, "singletop%nnlo_enable_heavy_decay", test)
          if (test) then
              enableIps(lord,7) = .true.
              initcalls(lord,7) = 2000000

              enableIps(nloVirt,3) = .true.
              initcalls(nloVirt,:) = 2000000

              enableIps(nloReal,9) = .true.
              initcalls(nloReal,:) = 8000000
          endif

      elseif (nproc == 1610) then
        enableIps(:,:) = .false.

        enableIps(lord,1) = .true.
        initcalls(lord,1) = 2000000

        call cfg_get(cfg, "singletop%nnlo_enable_interf_lxh", test)
        if (test .and. origKpart == knnlo) then
            ! real on light, real on heavy
            enableIps(nloReal,4) = .true.
            initcalls(nloReal,4) = 20000000

            ! real on light, virt on heavy line
            enableIps(nloReal,5) = .true.
            initcalls(nloReal,5) = 8000000

            ! virt on light line, real on heavy
            enableIps(nloVirt,4) = .true.
            initcalls(nloVirt,4) = 4000000

            ! virt on light, virt on heavy
            enableIps(nloVirt,5) = .true.
            initcalls(nloVirt,5) = 4000000
        endif

        call cfg_get(cfg, "singletop%nnlo_enable_interf_lxd", test)
        if (test .and. origKpart == knnlo) then
c           ! RR
            enableIps(nloReal,6) = .true.
            initcalls(nloReal,6) = 20000000

            ! RV
            enableIps(nloReal,7) = .true.
            initcalls(nloReal,7) = 4000000

            ! VR
            enableIps(nloVirt,6) = .true.
            initcalls(nloVirt,6) = 4000000

            ! VV
            enableIps(nloVirt,7) = .true.
            initcalls(nloVirt,7) = 4000000
        endif

        call cfg_get(cfg, "singletop%nnlo_enable_interf_hxd", test)
        if (test .and. origKpart == knnlo) then
            ! RR
            ! Run this with alpha=0.1 fully inclusively.
            ! With alpha=1 this is a cancellation nightmare and would
            ! in principle require separation into quark and gluon
            ! channels
            enableIps(nloReal,8) = .true.
            initcalls(nloReal,8) = 20000000

            ! RV
            enableIps(nloReal,9) = .true.
            initcalls(nloReal,9) = 4000000

            ! VR
            enableIps(nloVirt,8) = .true.
            initcalls(nloVirt,8) = 4000000

            ! VV
            enableIps(nloVirt,9) = .true.
            initcalls(nloVirt,9) = 4000000
        endif

        call cfg_get(cfg, "singletop%nnlo_enable_light", test)
        if (test) then
            enableIps(snloBelow,1:2) = .true.
            initcalls(snloBelow,1:2) = 1000000

            enableIps(snloAbove,1:4) = .true.
            initcalls(snloAbove,1:4) = 4000000

            enableIps(nloVirt,1) = .true.
            initcalls(nloVirt,1) = 4000000
            enableIps(nloReal,1) = .true.
            initcalls(nloReal,1) = 8000000

            enableIps(nnloBelow,1:2) = .true.
            initcalls(nnloBelow,1:2) = 4000000

            enableIps(nnloVirtAbove,1:2) = .true.
            initcalls(nnloVirtAbove,1:2) = 8000000

#define SINGLETOP_NNLO_REALCHANNELS 12
#if SINGLETOP_NNLO_REALCHANNELS == 66
            enableIps(nnloRealAbove,8:73) = .true.
            initcalls(nnloRealAbove,8:73) = 5000000

            ! difficult
            initcalls(nnloRealAbove,13) = 8000000
            initcalls(nnloRealAbove,15) = 8000000
#else
            enableIps(nnloRealAbove,8:19) = .true.
            initcalls(nnloRealAbove,8:19) = 15000000

            ! these are difficult
            initcalls(nnloRealAbove,9) = 20000000
            initcalls(nnloRealAbove,12) = 20000000
            initcalls(nnloRealAbove,15) = 20000000

#endif
        endif

        call cfg_get(cfg, "singletop%nnlo_enable_heavy_prod", test)
        if (test) then
            enableIps(snloBelow,3:4) = .true.
            initcalls(snloBelow,3:4) = 1000000

            enableIps(snloAbove,5:6) = .true.
            initcalls(snloAbove,5:6) = 4000000

            enableIps(nloVirt,2) = .true.
            initcalls(nloVirt,2) = 4000000
            enableIps(nloReal,2) = .true.
            initcalls(nloReal,2) = 4000000

            enableIps(nnloBelow,3:4) = .true.
            initcalls(nnloBelow,3:4) = 4000000

            enableIps(nnloVirtAbove,3:4) = .true.
            initcalls(nnloVirtAbove,3:4) = 4000000

            enableIps(nnloRealAbove,1:6) = .true.
            initcalls(nnloRealAbove,1:6) = 20000000
        endif

        call cfg_get(cfg, "singletop%nnlo_enable_heavy_decay", test)
        if (test) then
            enableIps(snloBelow,5) = .true.
            initcalls(snloBelow,5) = 1000000

            enableIps(snloAbove,7) = .true.
            initcalls(snloAbove,7) = 4000000

            enableIps(nloVirt,3) = .true.
            initcalls(nloVirt,3) = 4000000
            enableIps(nloReal,3) = .true.
            initcalls(nloReal,3) = 5000000

            enableIps(nnloBelow,5) = .true.
            initcalls(nnloBelow,5) = 4000000

            enableIps(nnloVirtAbove,5) = .true.
            initcalls(nnloVirtAbove,5) = 4000000

            enableIps(nnloRealAbove,7) = .true.
            initcalls(nnloRealAbove,7) = 20000000
        endif
      endif

      where (initcalls(:,:) /= 0) initcalls = initcalls * globalCallMult

      ! special case for nnlo singletop
      if (nproc == 1610) then
          ! for the rest compute coefficients
          coeffonly = .true.
          origCoeffonly = .true.
      endif

      ! additional contributions for Wilson coefficient calculation in Higgs+jet NNLO
      if (((nproc == 204) .or. (nprocbelow == 210))
     &    .and. origkpart == knnlo .and. knnlopart == 0) then
          computePart(nnloWilsonVirt) = .true.
          computePart(nnloWilsonReal) = .true.
          initcalls(nnloWilsonVirt,:) = initcalls(nloVirt,:)
          initcalls(nnloWilsonReal,:) = initcalls(nloReal,:)
      endif

c---
c--- setup calls/it
c---
      do j=1,maxParts
          do k=1,maxIPS
              iterationStorage(j,k)%vinfo%callsPerIt = initcalls(j,k)
          enddo
      enddo


c---      setup sobol generator
      ! for now we stride across mpi processes and use omp critical
      if (usesobol) then
          if (rank == 0) then
              write (*,*) "Using sobol with world_size = ", world_size
          endif
c         if (world_size > 1) then
c             write (*,*) "Striding sobol with world_size = ", world_size

c             if (world_size <= 2) then
c             strides = 1
c             else if (world_size == 4) then
c             strides = 2
c             else if (world_size == 8) then
c             strides = 3
c             else if (world_size == 16) then
c             strides = 4
c             else if (world_size == 32) then
c             strides = 5
c             else if (world_size == 64) then
c             strides = 6
c             else if (world_size == 128) then
c             strides = 7
c             else if (world_size == 256) then
c             strides = 8
c             else
c             write (*,*) "WARNING: Sobol sequence only usable with 2^n mpi processes."
c             write (*,*) "Falling back to using pseudo-random numbers."
c             usesobol = .false.
c             do part=1,maxParts; do ips=1,maxIps
c                 associate (info => iterationStorage(part,ips)%vinfo)
c                     info%useSobol = usesobol
c                 end associate
c             enddo; enddo
c             endif
c         else
c             strides = 0
c         endif

          if (ndim > 22) then
              write (*,*) "sobol direction vectors for ndim > ", ndim, " not implemented"
              error stop
          endif
      endif ! usesobol

      if (readin) then
          ! for now let's just read the file from all processes
          !if (rank == 0) then
              call deserializeMCFM()
          !endif
          !call mpi_broadcast_iterationStorage()
      endif

      ! to save whether power of alphas has been found for each part
      ! used in scale variation in lowint/realint/fragint
      foundPow(:) = .false.
c ================
c === warmup phase
c ================

      warmdone = .false.

      do part=1,maxParts
          if (computePart(part)) then
              ipsLordLoop: do ipsgen=1,maxIPS
              currentPart = part
              currentIps = ipsgen

              if (.not. enableIps(part,ipsgen)) cycle ipsLordLoop

              if (iterationStorage(currentPart,currentIps)%vinfo%warmupComplete) then
                  cycle ipsLordLoop
              endif

              warmdone = .true.

              fragint_mode = .false.
              expansionorder = 0
              ndimextra = 0

              select case (part)
                case (lord)
                    nproc = nprocbelow
                    abovecut = .false.
                    usescet = .false.
                    kpart = klord
                    coeffonly = origCoeffonly
                    ifun => lowint
                case (nloReal)
                    nproc = nprocbelow
                    abovecut = .false.
                    usescet = .false.
                    kpart = kreal
                    coeffonly = .false.
                    ifun => realint
                case (nloVirt)
                    nproc = nprocbelow
                    abovecut = .false.
                    usescet = .false.
                    kpart = kvirt
                    coeffonly = origCoeffonly
                    ifun => virtint
                case (nloFrag)
                    nproc = nprocbelow
                    abovecut = .false.
                    usescet = .false.
                    kpart = kfrag
                    coeffonly = origCoeffonly
                    fragint_mode = .true.
                    ifun => fragint
                case (nloRealExtra)
                    nproc = nprocextra
                    abovecut = .false.
                    usescet = .false.
                    kpart = kreal
                    coeffonly = .false.
                    ifun => realint
                case (snloBelow)
                    nproc = nprocbelow
                    abovecut = .false.
                    usescet = .true.
                    kpart = ksnlo
                    coeffonly = origCoeffonly
                    if (useQT_nnlo) then
                      ifun => qtint
                    elseif (usept) then
                      ifun => ptint
                    else
                      ifun => scetint
                    endif
                case (qtnloBelow)
                    nproc = nprocbelow
                    abovecut = .false.
                    usescet = .true.
                    kpart = ksnlo
                    coeffonly = origCoeffonly
                    ifun => qtsubint
                case (snloAbove)
                    nproc = nprocabove
                    abovecut = .true.
                    usescet = .true.
                    kpart = klord
                    coeffonly = origCoeffonly
                    ifun => lowint
                case (nloResAbove)
                    nproc = nprocabove
                    abovecut = .true.
                    usescet = .false.
                    kpart = klord
                    coeffonly = origCoeffonly
                    ifun => lowint
                case (nloResVetoVirt)
                    nproc = nprocbelow
                    abovecut = .false.
                    usescet = .false.
                    kpart = kvirt
                    coeffonly = origCoeffonly
                    ifun => virtint
                case (nloResVetoReal)
                    nproc = nprocbelow
                    abovecut = .false.
                    usescet = .false.
                    kpart = kreal
                    coeffonly = origCoeffonly
                    ifun => realint
                case (nnloResVetoBelow)
                    nproc = nprocbelow
                    abovecut = .false.
                    usescet = .true.
                    kpart = knnlo
                    ifun => ptint
                case (nnloResVetoVirtAbove)
                    nproc = nprocabove
                    abovecut = .true.
                    usescet = .true.
                    kpart = kvirt
!                    coeffonly = .true.
                    ifun => virtint
                case (nnloResVetoRealAbove)
                    nproc = nprocabove
                    abovecut = .true.
                    usescet = .true.
                    kpart = kreal
!                    coeffonly = .true.
                    ifun => realint
                case (nnloResVirtAbove)
                    nproc = nprocabove
                    abovecut = .true.
                    usescet = .false.
                    kpart = kvirt
                    coeffonly = .false.
                    ifun => virtint
                case (nnloResRealAbove)
                    nproc = nprocabove
                    abovecut = .true.
                    usescet = .false.
                    kpart = kreal
                    coeffonly = .true.
                    ifun => realint
                case (nnloBelow)
                    nproc = nprocbelow
                    abovecut = .false.
                    usescet = .true.
                    kpart = knnlo
                    coeffonly = .true.
                    if (useQT_nnlo) then
                      ifun => qtint
                    elseif (usept) then
                      ifun => ptint
                    else
                      ifun => scetint
                    endif
                case (qtnnloBelow)
                    nproc = nprocbelow
                    abovecut = .false.
                    usescet = .true.
                    kpart = knnlo
                    coeffonly = .true.
                    ifun => qtsubint
                case (nnloVirtAbove)
                    nproc = nprocabove
                    abovecut = .true.
                    usescet = .true.
                    kpart = kvirt
                    coeffonly = .true.
                    ifun => virtint
                case (nnloRealAbove)
                    nproc = nprocabove
                    abovecut = .true.
                    usescet = .true.
                    kpart = kreal
                    coeffonly = .true.
                    ifun => realint
                case (nloResummed)
                    nproc = nprocbelow
                    usescet = .false.
                    abovecut = .false.
                    kpart = klord
                    ifun => resint
                    if (usept) then
                      ndimextra=2
                    endif
                case (nloResexp)
                    nproc = nprocbelow
                    usescet = .false.
                    abovecut = .false.
                    kpart = klord
                    ifun => resexpint
                    if (usept) then
                      ndimextra = 2
                    endif
                case (qtn3loBelow)
                    nproc = nprocbelow
                    abovecut = .false.
                    usescet = .true.
                    kpart = kn3lo
                    coeffonly = origCoeffonly
                    ifun => qtsubint
                case (nnloWilsonReal)
                    nproc = nprocbelow
                    abovecut = .false.
                    usescet = .false.
                    kpart = kreal
                    coeffonly = .true.
                    expansionOrder = 2
                    ifun => realint
                case (nnloWilsonVirt)
                    nproc = nprocbelow
                    abovecut = .false.
                    usescet = .false.
                    kpart = kvirt
                    coeffonly = .true.
                    expansionOrder = 2
                    ifun => virtint
                case default
                    error stop "unknown part"
                end select

                reset = .true.
                scalereset = .true.

                call chooser
                ndim = ndim + ndim_incr(part) + ndimextra

                associate (info => iterationStorage(currentPart,currentIps)%vinfo )
                    if (info%useSobol) then
                    if (.not. info%sobolInitialized) then
                                    do i=1,ndim+3
                        call info%sstate(i)%initialize(s(i),a(i),m(:,i))
                        ! instead of striding we now skip 2^40 elements for each rank
                        do k=0, rank-1
                            skiprnd = info%sstate(i)%skip_ahead(info%sobolSkip + int(k,8)*2_int64**40_int64)
                        enddo
                        enddo
                        info%sobolInitialized = .true.
                    endif
                    endif
                end associate

                do
                  if (rank == 0) then
                  write (*,*) ""
                  write (*,'(A,A,I2)') trim(partStrings(part)), " warmup integration, contribution ", currentIps
                  write (*,*) ""
                  endif

                  iterationStorage(currentPart,currentIps)%used = .true.
                  call integrate(ifun, 0, iterBatchWarmup, ndim, iterationStorage(currentPart,currentIps)%vinfo)

                  associate (info => iterationStorage(currentPart,currentIps)%vinfo)
                  if ( info%si == 0._dp ) then
                      if (rank == 0) then
                          write (*,*) "Integral zero. Skipping this contribution."
                      endif
                      exit
                  else if ( (info%sd() / abs(info%sig()) > warmupPrecisionGoal) ) then
                      info%callsPerIt = info%callsPerIt * iterCallMult
                      if (rank == 0) then
                          write (*,'(A,I3,A)') "Relative warmup precision goal of ",
     &                          nint(warmupPrecisionGoal*100),
     &                          " percent not reached"
                          write (*,*) "Increasing calls to ", info%callsPerIt
                      endif
                  else if ( info%chisq() > warmupChisqGoal ) then
                      info%callsPerIt = info%callsPerIt * iterCallMult
                      if (rank == 0) then
                          write (*,'(A,F5.3,A)') "Chisq/it warmup precision goal of ",
     &                          warmupChisqgoal, " not reached"
                          write (*,*) "Increasing calls to ", info%callsPerIt
                      endif
                  else
                      if (rank == 0) then
                          write (*,*) "Reached warmup precisionGoal with ",
     &                          info%callsPerIt, " calls per iteration"
                      endif
                      ! warmup complete
                      exit
                  endif
                  end associate
                enddo

                ndim = ndim - ndim_incr(part) - ndimextra

                if ( iterationStorage(currentPart,currentIps)%vinfo%si == 0._dp ) then
                    ! skip main integration
                    iterationStorage(currentPart,currentIps)%vinfo%doFirstCall = .false.
                    iterationStorage(currentPart,currentIps)%used = .false.
                else
                    iterationStorage(currentPart,currentIps)%vinfo%doFirstCall = .true.
                endif

                iterationStorage(currentPart,currentIps)%vinfo%warmupComplete = .true.

                if (rank == 0) then
                    call serializeMCFM()
                endif


              enddo ipsLordLoop

          endif
      enddo

      if (rank == 0) then
          totalsig = 0._dp
          totalsd = 0._dp
          do part=1,maxParts
            do ips=1,maxIps
            if (iterationStorage(part,ips)%used) then
                totalsig = totalsig + iterationStorage(part,ips)%vinfo%sig()
                totalsd = totalsd + iterationStorage(part,ips)%vinfo%sd()**2
            endif
            enddo
          enddo
          totalsd = sqrt(totalsd)
          if (warmdone) then
            write (*,*) ""
            write (*,*) "WARMUP phase complete"
            write (*,*) "Intermediate warmup result"
          else
            write (*,*) "Result loaded from snapshot"
          endif
            call printcross(totalsig, totalsd, chisqMax())
            write (*,*) ""
      endif

c ===========================
c === final integration phase
c ===========================


      ! enable binning of histograms
      bin = .true.

      integrationLoop: do
        nextLoc(:) = 0

        ! first check if there has been a contribution
        ! with no final run, but only warmup integration
        loopContrib: do i=1,maxParts
            loopIpsgen: do j=1,maxIps
            if (iterationStorage(i,j)%vinfo%doFirstCall) then
                nextLoc(1) = i
                nextLoc(2) = j
                ! do not load result from saved file
                stage = 1
                if (rank == 0) then
                    write (*,'(A,A,A,I2)') "first full integration for ", trim(partStrings(i)), " contribution ", j
                endif
                iterBatch = iterBatch1
                exit loopContrib
            endif
            enddo loopIpsgen
        enddo loopContrib

        ! otherwise pick contribution with biggest sd
        if (nextLoc(1) == 0) then
            sdMod(:,:) = 0._dp
            do i=1,maxParts
            do j=1,maxIps
                if (iterationStorage(i,j)%used) then
                            ! updated approach:
                            ! weight errors with squared chisq
                            ! to prioritize largest chisq contributions
                            sdMod(i,j) = iterationStorage(i,j)%vinfo%sd() *
     &                              max(1d0,iterationStorage(i,j)%vinfo%chisq())
                endif
            enddo
            enddo

            ! we fake the "real emission" contributions uncertainty
            ! for the prioritization system here
            ! to boost the easier to compute contributions
            sdMod(nloReal,:) = sdMod(nloReal,:) / 1.5_dp
            sdMod(snloAbove,:) = sdMod(snloAbove,:) / 1.5_dp
            sdMod(nnloRealAbove,:) = sdMod(nnloRealAbove,:) / 1.5_dp
            sdMod(nloResummed,:) = sdMod(nloResummed,:) / 15._dp

            nextLoc = maxloc(sdMod(:,:))
            ! load previous result from saved file
            stage = 2

            if (rank == 0) then
            write (*,'(A,A,I2)') trim(partStrings(nextLoc(1))), " full integration, contribution ", nextLoc(2)
            endif

            iterBatch = iterBatch2
        endif

        currentPart = nextLoc(1)
        currentIps = nextLoc(2)

        if (.not. iterationStorage(currentPart,currentIps)%used) then
            error stop "something went wrong, this part should not be selected for integration"
        endif

        ! boost number of calls w.r.t. to warmup phase
        associate (info => iterationStorage(currentPart,currentIps)%vinfo)
            if (info%doFirstCall) then
            !info%callsPerIt = info%callsPerIt * 4
            info%callsPerIt = int(real(info%callsPerIt,dp) * callBoost)
            endif
                if (info%callsPerIt*iterCallMult > 500000000_int64) then
                  if (.not. info%doFirstCall) then
                      iterBatch = 1
                  endif
                endif

                if (info%doFirstCall) then
                    info%doFirstCall = .false.
                endif
        end associate

        ! Once we have reached maxCallsPerIter we just run one iteration
        ! per batch to have frequent snapshots and updated results.
        ! We set it here again for resumed runs.
        if (iterationStorage(currentPart,currentIps)%vinfo%callsPerIt*iterCallMult
     &              > 20000000000_int64) then
              iterationStorage(currentPart,currentIps)%vinfo%callsPerIt = 20000000000_int64
          iterBatch = 1
        endif

        ! now perform the full integration, similar to the warmup above

        reset = .true.
        scalereset = .true.
        fragint_mode = .false.
        expansionorder = 0

        ndimextra = 0
        select case (currentPart)
          case (lord)
              nproc = nprocbelow
              abovecut = .false.
              usescet = .false.
              kpart = klord
              coeffonly = origCoeffonly
              ifun => lowint
          case (nloReal)
              nproc = nprocbelow
              abovecut = .false.
              usescet = .false.
              kpart = kreal
              coeffonly = .false.
              ifun => realint
              ! bin this contribution to all tau bins
              if (origKpart == knnlo) then
!$omp parallel
              scetreweight(:) = 1._dp
!$omp end parallel
              endif
          case (nloVirt)
              nproc = nprocbelow
              abovecut = .false.
              usescet = .false.
              kpart = kvirt
              coeffonly = origCoeffonly
              ifun => virtint
              ! bin this contribution to all tau bins
              if (origKpart == knnlo) then
!$omp parallel
              scetreweight(:) = 1._dp
!$omp end parallel
              endif
          case (nloFrag)
              nproc = nprocbelow
              abovecut = .false.
              usescet = .false.
              kpart = kfrag
              coeffonly = origCoeffonly
              fragint_mode = .true.
              ifun => fragint
          case (nloRealExtra)
              nproc = nprocextra
              abovecut = .false.
              usescet = .false.
              kpart = kreal
              coeffonly = .false.
              ifun => realint
          case (snloBelow)
              nproc = nprocbelow
              abovecut = .false.
              usescet = .true.
              kpart = ksnlo
              coeffonly = origCoeffonly
              if (useQT_nnlo) then
                ifun => qtint
              elseif (usept) then
                ifun => ptint
              else
                ifun => scetint
              endif
          case (qtnloBelow)
              nproc = nprocbelow
              abovecut = .false.
              usescet = .true.
              kpart = ksnlo
              coeffonly = origCoeffonly
              ifun => qtsubint
          case (snloAbove)
              nproc = nprocabove
              abovecut = .true.
              usescet = .true.
              kpart = klord
              coeffonly = origCoeffonly
              ifun => lowint
          case (nloResAbove)
              nproc = nprocabove
              abovecut = .true.
              usescet = .false.
              kpart = klord
              coeffonly = origCoeffonly
              ifun => lowint
          case (nloResVetoVirt)
              nproc = nprocbelow
              abovecut = .false.
              usescet = .false.
              kpart = kvirt
              coeffonly = origCoeffonly
              ifun => virtint
          case (nloResVetoReal)
              nproc = nprocbelow
              abovecut = .false.
              usescet = .false.
              kpart = kreal
              coeffonly = origCoeffonly
              ifun => realint
          case (nnloResVetoBelow)
              nproc = nprocbelow
              abovecut = .false.
              usescet = .true.
              kpart = knnlo
              ifun => ptint
          case (nnloResVetoVirtAbove)
              nproc = nprocabove
              abovecut = .true.
              usescet = .true.
              kpart = kvirt
!              coeffonly = .true.
              ifun => virtint
          case (nnloResVetoRealAbove)
              nproc = nprocabove
              abovecut = .true.
              usescet = .true.
              kpart = kreal
!              coeffonly = .true.
              ifun => realint
          case (nnloResVirtAbove)
              nproc = nprocabove
              abovecut = .true.
              usescet = .false.
              kpart = kvirt
              coeffonly = .false.
              ifun => virtint
          case (nnloResRealAbove)
              nproc = nprocabove
              abovecut = .true.
              usescet = .false.
              kpart = kreal
              coeffonly = .true.
              ifun => realint
          case (nnloBelow)
              nproc = nprocbelow
              abovecut = .false.
              usescet = .true.
              kpart = knnlo
              coeffonly = .true.
              if (useQT_nnlo) then
                ifun => qtint
              elseif (usept) then
                ifun => ptint
              else
                ifun => scetint
              endif
          case (qtnnloBelow)
              nproc = nprocbelow
              abovecut = .false.
              usescet = .true.
              kpart = knnlo
              coeffonly = .true.
              ifun => qtsubint
          case (nnloVirtAbove)
              nproc = nprocabove
              abovecut = .true.
              usescet = .true.
              kpart = kvirt
              coeffonly = .true.
              ifun => virtint
          case (nnloRealAbove)
              nproc = nprocabove
              abovecut = .true.
              usescet = .true.
              kpart = kreal
              coeffonly = .true.
              ifun => realint
          case (nloResummed)
              nproc = nprocbelow
              abovecut = .false.
              usescet = .false.
              kpart = klord
              ifun => resint
              if (usept) then
                ndimextra=2
              endif
          case (nloResexp)
              nproc = nprocbelow
              usescet = .false.
              abovecut = .false.
              kpart = klord
              ifun => resexpint
              if (usept) then
                ndimextra=2
              endif
          case (qtn3loBelow)
              nproc = nprocbelow
              abovecut = .false.
              usescet = .true.
              kpart = kn3lo
              coeffonly = origCoeffonly
              ifun => qtsubint
          case (nnloWilsonReal)
              nproc = nprocbelow
              abovecut = .false.
              usescet = .false.
              kpart = kreal
              coeffonly = .true.
              expansionOrder = 2
              ifun => realint
          case (nnloWilsonVirt)
              nproc = nprocbelow
              abovecut = .false.
              usescet = .false.
              kpart = kvirt
              coeffonly = .true.
              expansionOrder = 2
              ifun => virtint
          case default
              error stop "unknown part"
          end select

          ! sanity setup. this includeTaucutgrid is initialized to .true.
          ! in parseinput.f. But when run in a full nnlo calculation
          ! after some scet runs the real integration might run again
          ! with includeTaucutgrid modified.
c          if (usescet) then
!$omp parallel
          includeTaucutgrid(:) = .true.
!$omp end parallel

c          endif

          reset = .true.
          scalereset = .true.

          ipsgen = currentIps

          call chooser
          ndim = ndim + ndim_incr(currentPart) + ndimextra

          associate (info => iterationStorage(currentPart,currentIps)%vinfo )
              if (info%useSobol) then
              if (.not. info%sobolInitialized) then
                          do i=1,ndim+3
                  call info%sstate(i)%initialize(s(i),a(i),m(:,i))
                  ! instead of striding we now skip 2^40 elements for each rank
                  do k=0, rank-1
                      skiprnd = info%sstate(i)%skip_ahead(info%sobolSkip + int(k,8)*2_int64**40_int64)
                  enddo
                  enddo
              endif
              endif
          end associate

          iterationStorage(currentPart,currentIps)%used = .true.
          call integrate(ifun, stage, iterBatch, ndim, iterationStorage(currentPart,currentIps)%vinfo)
          ndim = ndim - ndim_incr(currentPart) - ndimextra

          iterationStorage(currentPart,currentIps)%vinfo%callsPerIt =
     &          iterationStorage(currentPart,currentIps)%vinfo%callsPerIt * iterCallMult

          if (rank == 0) then
              call serializeMCFM()
          endif


          ! determine if we had at least one full integration for all parts
          ! then show final results, save histogram
          ! and exit if precision goal reached
          if (all(iterationStorage(:,:)%vinfo%doFirstCall .eqv. .false.)) then
              if (rank ==0) then
              call finalizeStorage
              if (nprocbelow == 1610) then
                  call finalizeStorageFixedOrder
              endif
              ! generate and save user readable output
              endif

              totalsig = 0._dp
              totalsd = 0._dp
              do part=1,maxParts
            do ips=1,maxIps
                if (iterationStorage(part,ips)%used) then
                    totalsig = totalsig + iterationStorage(part,ips)%vinfo%sig()
                    totalsd = totalsd + iterationStorage(part,ips)%vinfo%sd()**2
                endif
            enddo
              enddo
              totalsd = sqrt(totalsd)

              if ((abs(totalsd/totalsig) < resultPrecisionGoal) .or.
     &              abs(totalsd) < absResultPrecisionGoal ) then
                      exit integrationLoop
              elseif (get_walltime() > maxWalltime) then
                  if (rank == 0) then
                      write (*,*) "Maximum walltime of ", maxWalltime, " seconds reached"
                  endif
                  exit integrationLoop
              else ! print intermediate results
              if (rank == 0 .and. writeintermediate) then
                  write (*,*) "Intermediate full result"

                  if (nprocbelow == 1610) then
                  call printcross_singletop()
                  endif
                  call printallcross()

                  if (doPDFerrors) then
                    call printPDFuncertainties()
                  endif


                  if (doScalevar) then
                  call printScaleuncertainties()
                  endif

                  if (doMultitaucut) then
                  call printTaucuts()
                  endif

                  if (nprocbelow == 1610) then
                  call writeAllHistogramsSingletop()
                  !call writeAllHistograms()
                  else
                  call writeAllHistograms()
                  endif
              endif

              if (rank == 0) then
                  ! intermediate benchmark result
                  benchmarkinter: block
                  use MCFMBenchmark
                  integer :: exitCode
                  if (cfg_var_configadded(cfg, "extra%benchmark")) then
                      exitCode = comparisonCode(1)
                  endif
                  end block benchmarkinter
              endif
              endif

          endif

          ! save snapshots

      enddo integrationLoop

      integ = totalsig
      integ_err = totalsd

      end subroutine

