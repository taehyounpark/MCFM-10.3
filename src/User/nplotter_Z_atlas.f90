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

          allocate(histos(3))

          if (rank == 0) then
              write (*,*) "RESUMMATION: Using transition with switch ", transitionSwitch
          endif

          qtvals = [1.0000d0,1.2589d0,1.5849d0,1.9953d0, &
                2.5119d0,3.1623d0,3.9811d0,5.0119d0,6.3096d0,7.9433d0, &
                10.0000d0,12.5893d0,15.8489d0,19.9526d0,25.1189d0, &
                31.6228d0,39.8107d0,50.1187d0,63.0957d0,79.4328d0,100.0000d0, &
                125.89d0, 158.49d0, 199.53d0, 251.19d0, 316.23d0, 398.11d0, &
                501.19d0, 630.96d0, 794.33d0]

          histos(1) = plot_setup_custom(qtvals,'pt34_fine_notrans')

          histos(2) = plot_setup_custom([0d0,2.5d0,5.0d0,8.0d0,11.4d0, &
              14.9d0,18.5d0,22.0d0,25.5d0,29.0d0,32.6d0,36.4d0,40.40d0, &
              44.9d0, 50.2d0, 56.4d0, 63.9d0, 73.4d0, 85.4d0, 105.0d0, &
              132.0d0, 173.0d0, 253.0d0, 4000d0], 'pt34_atlas13_trans04')

          histos(3) = plot_setup_custom([0d0,2.5d0,5.0d0,8.0d0,11.4d0, &
              14.9d0,18.5d0,22.0d0,25.5d0,29.0d0,32.6d0,36.4d0,40.40d0, &
              44.9d0, 50.2d0, 56.4d0, 63.9d0, 73.4d0, 85.4d0, 105.0d0, &
              132.0d0, 173.0d0, 253.0d0, 4000d0], 'pt34_atlas13_trans06')

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

          !pt3 = pt(3,p)
          pt34 = pttwo(3,4,p)

          !delphi34 = delphi(p(3,:),p(4,:))
          !phiacop = 2._dp*atan(sqrt((1._dp+cos(delphi34))/(1._dp-cos(delphi34))))
          !costhetastar = tanh((etarap(3,p)-etarap(4,p))/2._dp)
          !phistar = tan(phiacop/2._dp)*sin(acos(costhetastar))


!         m34 = puremass(p(3,:)+p(4,:))

          if (origKpart == kresummed) then
              if (abovecut .eqv. .false.) then
                  trans04 = transition((pt34/twomass(3,4,p))**2d0,0.001d0, 0.4d0 ,0.001d0)
                  trans06 = transition((pt34/twomass(3,4,p))**2d0,0.001d0, 0.6d0 ,0.001d0)
              else
                  ! fo piece without transition
                  trans04 = 1._dp
                  trans06 = 1._dp
              endif
          else
              trans04 = 1._dp
              trans06 = 1._dp
          endif

          if (ieee_is_nan(pt34)) then
              pt34 = -1._dp
          endif

          if (ieee_is_nan(phistar)) then
              phistar = -1._dp
          endif

          ids = histos
          vals = [pt34,pt34,pt34]
          wts = [wt,wt*trans04,wt*trans06]

      end subroutine
      
      
end module
