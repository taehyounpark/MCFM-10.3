!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module nplotter_ZZ
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
          implicit none

          allocate(histos(10))

          histos(1) = plot_setup_uniform(0._dp,200._dp,2._dp,'pt3456_lin')
          histos(2) = plot_setup_uniform(0._dp,500._dp,10._dp,'m3456')

          histos(3) = plot_setup_custom([ &
                0.0100d0,0.0126d0,0.0158d0,0.0200d0,0.0251d0,0.0316d0, &
                0.0398d0,0.0501d0,0.0631d0,0.0794d0,0.1000d0,0.1259d0, &
                0.1585d0,0.1995d0,0.2512d0,0.3162d0,0.3981d0,0.5012d0, &
                0.6310d0,0.7943d0,1.0000d0,1.2589d0,1.5849d0,1.9953d0, &
                2.5119d0,3.1623d0,3.9811d0,5.0119d0,6.3096d0,7.9433d0, &
                10.0000d0,12.5893d0,15.8489d0,19.9526d0,25.1189d0, &
                31.6228d0,39.8107d0,50.1187d0,63.0957d0,79.4328d0,100.0000d0, &
                125.893d0,158.489d0,199.526d0,251.189d0, &
                316.228d0,398.107d0,501.187d0,630.957d0,794.328d0,1000.000d0], &
                'pt3456')

! ZZ Fig. 5 (left) from arXiv:2009.01186
! note that last bin is extended (from 300) to include overflow
          histos(4) = plot_setup_custom([0d0,25d0,50d0,75d0,100d0,150d0,200d0,14000d0],'ptZZ')

! ZZ Fig. 4 (left) from arXiv:2009.01186
          histos(5) = plot_setup_uniform(0._dp,150._dp,15._dp,'ptl1')
          histos(6) = plot_setup_uniform(0._dp,150._dp,15._dp,'ptl2')
          histos(7) = plot_setup_uniform(0._dp,150._dp,15._dp,'ptl3')
          histos(8) = plot_setup_uniform(0._dp,150._dp,15._dp,'ptl4')

! ZZ Fig. 15(a) from arXiv:2103.01918
          histos(9) = plot_setup_custom([ 20d0,86d0,96d0,104d0,122d0,128d0,137d0,149d0, &
                       163d0,182d0,191d0,199d0,216d0,234d0,253d0,274d0,297d0,322d0, &
                       350d0,381d0,436d0,2000d0 ],'m4lpta')

! Cross-section in the onshell region, 180 < m4l < 2000
          histos(10) = plot_setup_custom([ 180d0,2000d0 ],'m4lonshell')

      end subroutine

      subroutine book(p,wt,ids,vals,wts)
          use types
          use ResummationTransition
          use ieee_arithmetic
          use ptveto, only: usept
          implicit none
          include 'mxpart.f'
          include 'kpart.f'
          include 'taucut.f'! abovecut
          include 'constants.f'
          include 'masses.f'
          include 'nwz.f'

          real(dp), intent(in) :: p(mxpart,4)
          real(dp), intent(in) :: wt
          
          integer, allocatable, intent(out) :: ids(:)
          real(dp), allocatable, intent(out) :: vals(:)
          real(dp), allocatable, intent(out) :: wts(:)
          real(dp) :: pt3456,pt34,ptl
          real(dp) :: ptfour,pt,delphi,wta,wtb,wtc,wtd,wte
          real(dp) :: puremass, mZZ,m12,m34,trans,ptvals(mxpart),tmp
!          real(dp) :: transm3456(9)
          integer :: id(mxpart)

          pt3456 = ptfour(3,4,5,6,p)
          mZZ = puremass(p(3,:)+p(4,:)+p(5,:)+p(6,:))
          ptvals(1) = pt(3,p)
          ptvals(2) = pt(4,p)
          ptvals(3) = pt(5,p)
          ptvals(4) = pt(6,p)
          call arraysort(4,ptvals,id)

! pt slices in Fig. 15 from arXiv:2103d01918
          wta=0d0
!          wtb=0d0
!          wtc=0d0
!          wtd=0d0
!          wte=0d0
          if (pt3456 < 10d0) wta=wt
!          if ((pt3456 > 10d0) .and. (pt3456 < 20d0)) wtb=wt
!          if ((pt3456 > 20d0) .and. (pt3456 < 50d0)) wtc=wt
!          if ((pt3456 > 50d0) .and. (pt3456 < 100d0)) wtd=wt
!          if ((pt3456 > 100d0) .and. (pt3456 < 600d0)) wte=wt

!          transm3456(:)=1._dp
! transition function is only applied when computing 'resLO', 'resNLO', 'resNNLO' (and +'p')
          if ((origKpart == kresummed) .and. (krespart == 0)) then
              if (abovecut .eqv. .false.) then
                  trans = transition((pt3456/mZZ)**2d0,0.001d0,transitionSwitch,0.001d0)
!                  transm3456(1) = transition((pt3456/mZZ)**2d0,0.001d0,0.05d0,0.001d0)
!                  transm3456(2) = transition((pt3456/mZZ)**2d0,0.001d0,0.20d0,0.001d0)
              else
                  ! fo piece without transition
                  trans = 1._dp
              endif
          else
              trans = 1._dp
          endif

! do not apply any transition function for jet veto resummation
          if (usept) trans = 1._dp

          if (ieee_is_nan(pt3456)) then
              pt3456 = -1._dp
          endif

          ids = histos
          vals = [pt3456,mZZ,pt3456,pt3456,ptvals(id(1)),ptvals(id(2)),ptvals(id(3)),ptvals(id(4)),mZZ,mZZ]
          wts = [wt*trans,wt*trans,wt*trans,wt*trans,wt*trans,wt*trans,wt*trans,wt*trans,wta*trans,wt*trans]

      end subroutine

end module
