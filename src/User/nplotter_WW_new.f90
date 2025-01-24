!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module nplotter_WW
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

          allocate(histos(18))

          histos(1) = plot_setup_uniform(0._dp,200._dp,2._dp,'pt3456_lin')
          histos(2) = plot_setup_uniform(0._dp,500._dp,10._dp,'m3456')
          histos(3) = plot_setup_uniform(0._dp,500._dp,10._dp,'m34m56')
          histos(4) = plot_setup_custom([ &
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
          histos(5) = plot_setup_custom([0.0010d0,0.0013d0,0.0016d0,0.0020d0, &
                0.0025d0,0.0032d0,0.0040d0,0.0050d0,0.0063d0,0.0079d0, &
                0.0100d0,0.0126d0,0.0158d0,0.0200d0,0.0251d0,0.0316d0, &
                0.0398d0,0.0501d0,0.0631d0,0.0794d0,0.1000d0,0.1259d0, &
                0.1585d0,0.1995d0,0.2512d0,0.3162d0,0.3981d0,0.5012d0, &
                0.6310d0,0.7943d0,1.0000d0,1.2589d0,1.5849d0,1.9953d0, &
                2.5119d0,3.1623d0,3.9811d0,5.0119d0,6.3096d0,7.9433d0, &
                10.0000d0,12.5893d0,15.8489d0,19.9526d0,25.1189d0, &
                31.6228d0,39.8107d0,50.1187d0,63.0957d0,79.4328d0,100.0000d0], &
                'dlogpt3456')

         histos(6) = plot_setup_uniform(0._dp,200._dp,2._dp,'pt3456_transc')
         histos(7) = plot_setup_uniform(0._dp,200._dp,2._dp,'pt3456_trans1')
         histos(8) = plot_setup_uniform(0._dp,200._dp,2._dp,'pt3456_trans2')

         histos(9) = plot_setup_custom([0.0010d0,0.0013d0,0.0016d0,0.0020d0, &
                0.0025d0,0.0032d0,0.0040d0,0.0050d0,0.0063d0,0.0079d0, &
                0.0100d0,0.0126d0,0.0158d0,0.0200d0,0.0251d0,0.0316d0, &
                0.0398d0,0.0501d0,0.0631d0,0.0794d0,0.1000d0,0.1259d0, &
                0.1585d0,0.1995d0,0.2512d0,0.3162d0,0.3981d0,0.5012d0, &
                0.6310d0,0.7943d0,1.0000d0],'pt3456overm3456sq_transc')
         histos(10) = plot_setup_custom([0.0010d0,0.0013d0,0.0016d0,0.0020d0, &
                0.0025d0,0.0032d0,0.0040d0,0.0050d0,0.0063d0,0.0079d0, &
                0.0100d0,0.0126d0,0.0158d0,0.0200d0,0.0251d0,0.0316d0, &
                0.0398d0,0.0501d0,0.0631d0,0.0794d0,0.1000d0,0.1259d0, &
                0.1585d0,0.1995d0,0.2512d0,0.3162d0,0.3981d0,0.5012d0, &
                0.6310d0,0.7943d0,1.0000d0],'pt3456overm3456sq_trans1')
         histos(11) = plot_setup_custom([0.0010d0,0.0013d0,0.0016d0,0.0020d0, &
                0.0025d0,0.0032d0,0.0040d0,0.0050d0,0.0063d0,0.0079d0, &
                0.0100d0,0.0126d0,0.0158d0,0.0200d0,0.0251d0,0.0316d0, &
                0.0398d0,0.0501d0,0.0631d0,0.0794d0,0.1000d0,0.1259d0, &
                0.1585d0,0.1995d0,0.2512d0,0.3162d0,0.3981d0,0.5012d0, &
                0.6310d0,0.7943d0,1.0000d0],'pt3456overm3456sq_trans2')
! WW
         histos(12) = plot_setup_uniform(0._dp,600._dp,10._dp,'mTWW_transc')
         histos(13) = plot_setup_uniform(0._dp,600._dp,10._dp,'mTWW_trans1')
         histos(14) = plot_setup_uniform(0._dp,600._dp,10._dp,'mTWW_trans2')

! WW Fig. 7 (top left) from arXiv:2009.00119
          histos(15) = plot_setup_custom([ 55d0, 75d0 ,85d0 ,95d0 ,110d0,125d0,140d0,160d0, &
                        185d0,220d0,280d0,380d0,600d0,1500d0 ],'mll')
! WW Fig. 7 (top right) from arXiv:2009.00119
          histos(16) = plot_setup_custom([ 27d0,40d0,50d0,60d0,70d0,80d0,90d0,100d0, &
                        110d0,130d0,150d0,175d0,220d0,300d0,400d0 ],'ptlmax')
! WW Fig. 7 (bottom left) from arXiv:2009.00119
          histos(17) = plot_setup_custom([ 25d0,30d0,35d0,40d0,45d0,50d0,75d0,100d0,150d0],'ptlmin')
! WW Fig. 7 (bottom right) from arXiv:2009.00119
          histos(18) = plot_setup_custom([ 0d0,0.34907d0,0.69813d0,1.0472d0,1.3963d0, &
                        1.7453d0,2.0944d0,2.4435d0,2.7925d0,3.1416d0 ],'dphill')

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
          real(dp) :: pt3456,pt34,mll,dphill
          real(dp) :: ptfour,pttwo,pt,delphi
          real(dp) :: puremass, mWW,mWmW,m12,m34,trans,ptvals(mxpart),ptlmax,ptlmin
          real(dp) :: transm3456(9),pt36,mTWW
          integer :: id(mxpart)

          pt3456 = ptfour(3,4,5,6,p)
          m12 = puremass(p(3,:)+p(4,:))
          m34 = puremass(p(5,:)+p(6,:))
          mll = puremass(p(4,:)+p(5,:))
          mWmW = m12+m34
          mWW = puremass(p(3,:)+p(4,:)+p(5,:)+p(6,:))
          dphill=delphi(p(4,:),p(5,:))
          ptvals(1) = pt(3,p)
          ptvals(2) = pt(4,p)
          ptvals(3) = pt(5,p)
          ptvals(4) = pt(6,p)
          if (ptvals(2) > ptvals(3)) then
            ptlmax=ptvals(2)
            ptlmin=ptvals(3)
          else
            ptlmax=ptvals(3)
            ptlmin=ptvals(2)
          endif

! Definition of 1902.05759, Section 7, extended to Etmiss = pt36
          pt36=pttwo(3,6,p)
          mTWW = sqrt(max(0d0,(ptvals(2)+ptvals(3)+pt36)**2-pt3456**2))

          transm3456(:)=1._dp
! transition function is only applied when computing 'resLO', 'resNLO', 'resNNLO' (and +'p')
          if ((origKpart == kresummed) .and. (krespart == 0)) then
              if (abovecut .eqv. .false.) then
                  trans = transition((pt3456/mWW)**2d0,0.001d0,transitionSwitch,0.001d0)
                  transm3456(1) = transition((pt3456/mWW)**2d0,0.001d0,0.05d0,0.001d0)
                  transm3456(2) = transition((pt3456/mWW)**2d0,0.001d0,0.20d0,0.001d0)
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
          vals = [pt3456,mWW,mWmW,pt3456,pt3456, &
                  pt3456,pt3456,pt3456, &
                  (pt3456/mWW)**2,(pt3456/mWW)**2,(pt3456/mWW)**2, &
                  mTWW,mTWW,mTWW, &
                  mll,ptlmax,ptlmin,dphill ]
          wts = [wt*trans,wt*trans,wt*trans,wt*trans,wt*trans*pt3456, &
                 wt*trans,wt*transm3456(1),wt*transm3456(2), &
                 wt*trans,wt*transm3456(1),wt*transm3456(2), &
                 wt*trans,wt*transm3456(1),wt*transm3456(2), &
                 wt*trans,wt*trans,wt*trans,wt*trans ]

      end subroutine

end module
