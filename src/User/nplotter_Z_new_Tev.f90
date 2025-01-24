!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

module nplotter_Z
      use types
      use MCFMPlotting
      use ResummationTransition, only: transition
      use qtResummation_params, only: transitionSwitch
      implicit none

      public :: setup, book
      private

      integer, save, allocatable :: histos(:)

      logical, save :: resabove

      contains

      subroutine setup()
          use types
          use parseinput
          implicit none
          include 'src/Inc/mpicommon.f'
          include 'src/Inc/runstring.f'

          real(dp), allocatable :: qtvals(:)
          real(dp), allocatable :: qtvals_large(:)

          character(len=255) :: rundir
          common/rundir/rundir

          allocate(histos(2))

          if (rank == 0) then
              write (*,*) "RESUMMATION: Using transition with switch ", transitionSwitch
          endif

          qtvals = [1.0000d0,1.2589d0,1.5849d0,1.9953d0, &
                2.5119d0,3.1623d0,3.9811d0,5.0119d0,6.3096d0,7.9433d0, &
                10.0000d0,12.5893d0,15.8489d0,19.9526d0,25.1189d0, &
                31.6228d0,39.8107d0,50.1187d0,63.0957d0,79.4328d0,100.0000d0, &
                125.89d0, 158.49d0, 199.53d0, 251.19d0, 316.23d0, 398.11d0, &
                501.19d0, 630.96d0, 794.33d0]

          histos(1) = plot_setup_custom(qtvals,'pt34_fine')

          histos(2) = plot_setup_custom([00.0d0, 000.5d0, 001.0d0, 001.5d0, 002.0d0, 002.5d0, 003.0d0, 003.5d0, 004.0d0, 004.5d0, &
                    005.0d0, 005.5d0, 006.0d0, 006.5d0, 007.0d0, 007.5d0, 008.0d0, 008.5d0, 009.0d0, 009.5d0, 010.0d0, 010.5d0, &
                    011.0d0, 011.5d0, 012.0d0, 012.5d0, 013.0d0, 013.5d0, 014.0d0, 014.5d0, 015.0d0, 015.5d0, 016.0d0, 016.5d0, &
                    017.0d0, 017.5d0, 018.0d0, 018.5d0, 019.0d0, 019.5d0, 020.0d0, 020.5d0, 021.0d0, 021.5d0, 022.0d0, 022.5d0, &
                    023.0d0, 023.5d0, 024.0d0, 024.5d0, 025.0d0, 026.0d0, 027.0d0, 028.0d0, 029.0d0, 030.0d0, 032.0d0, 034.0d0, &
                    036.0d0, 038.0d0, 040.0d0],'pt34_tev')            

      end subroutine

      subroutine book(p,wt,ids,vals,wts)
          use types
          use ResummationTransition
          use ieee_arithmetic
          implicit none
          include 'src/Inc/mxpart.f'
          include 'src/Inc/kpart.f'
          include 'src/Inc/taucut.f' ! abovecut

          real(dp), intent(in) :: p(mxpart,4)
          real(dp), intent(in) :: wt

          integer, allocatable, intent(out) :: ids(:)
          real(dp), allocatable, intent(out) :: vals(:)
          real(dp), allocatable, intent(out) :: wts(:)

          real(dp) :: pt,pttwo, ptpure, puremass, twomass, delphi, etarap

          real(dp) :: pt3, pt34
          real(dp) :: trans04, trans06
          real(dp) :: phistar, phiacop, costhetastar, delphi34
          real(dp) :: costhetaCS, sinthetaCS
          real(dp) :: r, m34

!         pt3 = pt(3,p)
          pt34 = pttwo(3,4,p)

!         delphi34 = delphi(p(3,:),p(4,:))
!         phiacop = 2._dp*atan(sqrt((1._dp+cos(delphi34))/(1._dp-cos(delphi34))))
!         costhetastar = tanh((etarap(3,p)-etarap(4,p))/2._dp)
!         phistar = tan(phiacop/2._dp)*sin(acos(costhetastar))


!         m34 = puremass(p(3,:)+p(4,:))

!         if (origKpart == kresummed) then
!             if (abovecut .eqv. .false.) then
!                 trans04 = transition((pt34/twomass(3,4,p))**2d0,0.001d0, 0.4d0 ,0.001d0)
!                 trans06 = transition((pt34/twomass(3,4,p))**2d0,0.001d0, 0.6d0 ,0.001d0)
!             else
!                 ! fo piece without transition
!                 trans04 = 1._dp
!                 trans06 = 1._dp
!             endif
!         else
!             trans04 = 1._dp
!             trans06 = 1._dp
!         endif

          if (ieee_is_nan(pt34)) then
              pt34 = -1._dp
          endif

!         if (ieee_is_nan(phistar)) then
!             phistar = -1._dp
!         endif

          ids = histos
          vals = [pt34,pt34]
          wts = [wt,wt]

      end subroutine
      
      
end module
