!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module nplotter_ZH_WW
      use types
      use MCFMPlotting
      use qtResummation_params, only: transitionSwitch
      implicit none

      public :: setup, book
      private

      integer, save, allocatable :: histos(:)

      contains

      subroutine setup()
          use types
          use parseinput
          implicit none

          include 'mpicommon.f'

          allocate(histos(2))

          if (rank == 0) then
              write (*,*) "RESUMMATION: Using transition with switch ", transitionSwitch
          endif

          histos(1) = plot_setup_custom([0.0010d0,0.0013d0,0.0016d0,0.0020d0, &
                0.0025d0,0.0032d0,0.0040d0,0.0050d0,0.0063d0,0.0079d0, &
                0.0100d0,0.0126d0,0.0158d0,0.0200d0,0.0251d0,0.0316d0, &
                0.0398d0,0.0501d0,0.0631d0,0.0794d0,0.1000d0,0.1259d0, &
                0.1585d0,0.1995d0,0.2512d0,0.3162d0,0.3981d0,0.5012d0, &
                0.6310d0,0.7943d0,1.0000d0,1.2589d0,1.5849d0,1.9953d0, &
                2.5119d0,3.1623d0,3.9811d0,5.0119d0,6.3096d0,7.9433d0, &
                10.0000d0,12.5893d0,15.8489d0,19.9526d0,25.1189d0, &
                31.6228d0,39.8107d0,50.1187d0,63.0957d0,79.4328d0,100.0000d0, &
                125.8925d0], 'pt34_fine')
          histos(2) = plot_setup_uniform(0d0,80d0,1d0,'pt34_uni')

      end subroutine

      subroutine book(p,wt,ids,vals,wts)
          use types
          use ResummationTransition
          use SCET, only: currentNd
          use ieee_arithmetic
          implicit none
          include 'mxpart.f'
          include 'kpart.f'
          include 'taucut.f'! abovecut
          include 'frag.f'

          real(dp), intent(in) :: p(mxpart,4), wt

          integer, allocatable, intent(out) :: ids(:)
          real(dp), allocatable, intent(out) :: vals(:)
          real(dp), allocatable, intent(out) :: wts(:)

          real(dp) :: puremass,ptpure

          real(dp) :: pparton(mxpart,4)
          real(dp) :: ptZH, mZH, trans

          call getptilde(currentNd,pparton)
          ptZH = ptpure(pparton(3,:)+pparton(4,:)+pparton(5,:)+ &
                        pparton(6,:)+pparton(7,:)+pparton(8,:))

          mZH = puremass(pparton(3,:)+pparton(4,:)+pparton(5,:)+ &
                        pparton(6,:)+pparton(7,:)+pparton(8,:))

          if (origKpart == kresummed) then
              if (abovecut .eqv. .false.) then
                  trans = transition((ptZH/mZH)**2d0,0.001d0,transitionSwitch,0.001d0)
              else
                  ! fo piece without transition
                  trans = 1._dp
              endif
          else
              trans = 1._dp
          endif

          if (ieee_is_nan(ptZH)) then
              ptZH = -1._dp
          endif

          ids = histos

          ! DEBUG
          !trans = 1._dp

          vals = [ptZH,ptZH]
          wts = [wt*trans,wt*trans]

      end subroutine

end module
