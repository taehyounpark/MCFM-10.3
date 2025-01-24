!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module nplotter_WZ
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

          allocate(histos(4))

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

! WZ Fig. 16 (top row) from arXiv:2110.11231
          histos(4) = plot_setup_custom([ 25d0,35d0,50d0,70d0,90d0,110d0,130d0,160d0,200d0,300d0 ],'ptl')

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
          include 'masses.f'
          include 'nwz.f'

          real(dp), intent(in) :: p(mxpart,4)
          real(dp), intent(in) :: wt
          
          integer, allocatable, intent(out) :: ids(:)
          real(dp), allocatable, intent(out) :: vals(:)
          real(dp), allocatable, intent(out) :: wts(:)
          real(dp) :: pt3456,pt34,ptl
          real(dp) :: ptfour,pt
          real(dp) :: puremass, mWZ,trans
          integer :: id(mxpart)

          pt3456 = ptfour(3,4,5,6,p)
          mWZ = puremass(p(3,:)+p(4,:)+p(5,:)+p(6,:))
          if (nwz == +1) then
            ptl=pt(4,p)
          else
            ptl=pt(3,p)
          endif

! transition function is only applied when computing 'resLO', 'resNLO', 'resNNLO' (and +'p')
          if ((origKpart == kresummed) .and. (krespart == 0)) then
              if (abovecut .eqv. .false.) then
                  trans = transition((pt3456/mWZ)**2d0,0.001d0,transitionSwitch,0.001d0)
              else
                  ! fo piece without transition
                  trans = 1._dp
              endif
          else
              trans = 1._dp
          endif

          if (ieee_is_nan(pt3456)) then
              pt3456 = -1._dp
          endif

          ids = histos
          vals = [pt3456,mWZ,pt3456,ptl]
          wts = [wt*trans,wt*trans,wt*trans,wt*trans]

      end subroutine

end module
