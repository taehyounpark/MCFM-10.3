!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module nplotter_ZGamma
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

          include 'mpicommon.f'

          allocate(histos(6))

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
                31.6228d0,39.8107d0,50.1187d0,63.0957d0,79.4328d0,100.0000d0], &
                'pt345_fine')
          histos(2) = plot_setup_uniform(0._dp,60._dp,1.0_dp,'pt345_coarse')
          histos(3) = plot_setup_custom([0d0,5d0,10d0,15d0,20d0,30d0, &
              40d0,50d0,70d0,100d0,200d0,1500d0], 'ptllgamma_cms')
          histos(4) = plot_setup_custom([0d0,0.1d0,0.2d0,0.3d0,0.4d0, &
              0.5d0,0.6d0,0.7d0,0.8d0,0.9d0,0.95d0,1d0],'delphi345_cms')
          histos(5) = plot_setup_uniform(0._dp,500._dp,10._dp,'m345')
          histos(6) = plot_setup_uniform(1d0,81d0,4d0,'pt345_andrin')

      end subroutine

      subroutine book(p,wt,ids,vals,wts)
          use types
          use ResummationTransition
          use ieee_arithmetic
          implicit none
          include 'mxpart.f'
          include 'kpart.f'
          include 'taucut.f'! abovecut
          include 'constants.f'

          real(dp), intent(in) :: p(mxpart,4)
          real(dp), intent(in) :: wt

          integer, allocatable, intent(out) :: ids(:)
          real(dp), allocatable, intent(out) :: vals(:)
          real(dp), allocatable, intent(out) :: wts(:)

          real(dp) :: pttwo,delphi,threemass
          real(dp) :: ptthree

          real(dp) :: pt34, pt345, delphi345, trans, m345

          pt34 = pttwo(3,4,p)
          if (ieee_is_nan(pt34)) pt34 = -1._dp
          pt345 = ptthree(3,4,5,p)
          if (ieee_is_nan(pt345)) pt345 = -1._dp
          delphi345 = delphi(p(3,:)+p(4,:), p(5,:))
          if (ieee_is_nan(delphi345)) delphi345 = -1._dp
          m345 = threemass(3,4,5,p)
          if (ieee_is_nan(m345)) m345 = -1._dp
!         delphi34 = delphi(p(3,:),p(4,:))
!         phiacop = 2._dp*atan(sqrt((1._dp+cos(delphi34))/(1._dp-cos(delphi34))))
!         costhetastar = tanh((etarap(3,p)-etarap(4,p))/2._dp)
!         phistar = tan(phiacop/2._dp)*sin(acos(costhetastar))

          if (origKpart == kresummed) then
              if (abovecut .eqv. .false.) then
                  trans = transition((pt345/threemass(3,4,5,p))**2d0,0.001d0,transitionSwitch,0.001d0)
              else
                  ! fo piece without transition
                  trans = 1._dp
              endif
          else
              trans = 1._dp
          endif

           !DEBUG
          !trans = 1d0

          ids = histos
          vals = [pt345,pt345,pt345,delphi345/pi,m345,pt345]
          wts = [trans*wt,wt*trans,trans*wt,trans*wt,trans*wt,trans*wt]

      end subroutine

end module
