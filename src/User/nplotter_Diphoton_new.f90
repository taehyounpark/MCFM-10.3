!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module nplotter_Diphoton
      use types
      use MCFMPlotting
      use qtResummation_params, only: transitionSwitch
      implicit none

      public :: setup, book
      private

      integer, save, allocatable :: histos(:)

      integer, save :: study

      contains

      function crosspr(a,b)
          use types
          implicit none
          real(dp) :: crosspr(3)
          real(dp), intent(in) :: a(3), b(3)

          crosspr = [a(2)*b(3) - a(3)*b(2), &
                     a(3)*b(1) - a(1)*b(3), &
                     a(1)*b(2) - a(2)*b(1)]
      end function


      subroutine setup()
          use types
          use parseinput
          implicit none
          include 'runstring.f'
          include 'mpicommon.f'

          allocate(histos(6))

          if (trim(runstring) == "1704.03839") then
              study = 170403839
          else
              study = 12111913
          endif

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
          histos(2) = plot_setup_custom([0.0010d0,0.0013d0,0.0016d0,0.0020d0, &
                0.0025d0,0.0032d0,0.0040d0,0.0050d0,0.0063d0,0.0079d0, &
                0.0100d0,0.0126d0,0.0158d0,0.0200d0,0.0251d0,0.0316d0, &
                0.0398d0,0.0501d0,0.0631d0,0.0794d0,0.1000d0,0.1259d0, &
                0.1585d0,0.1995d0,0.2512d0,0.3162d0,0.3981d0,0.5012d0, &
                0.6310d0,0.7943d0,1.0000d0,1.2589d0,1.5849d0,1.9953d0, &
                2.5119d0,3.1623d0,3.9811d0,5.0119d0,6.3096d0,7.9433d0, &
                10.0000d0,12.5893d0,15.8489d0,19.9526d0,25.1189d0, &
                31.6228d0,39.8107d0,50.1187d0,63.0957d0,79.4328d0,100.0000d0, &
                125.8925d0,158.4893d0,199.5262d0,251.1886d0,316.2278d0, &
                398.1072d0,501.1872d0,630.9573d0,794.3282d0,1000.0000d0, &
                1258.9254d0,1584.8932d0,1995.2623d0], 'm34_fine')
          if (study == 12111913) then
              histos(3) = plot_setup_custom([0d0,2d0,4d0,6d0,8d0,10d0,12d0,14d0, &
                  16d0,18d0,20d0,25d0,30d0,35d0,40d0,45d0,50d0,55d0,60d0,65d0, &
                  70d0,75d0,80d0,90d0,100d0,110d0,120d0,130d0,140d0,150d0,175d0, &
                  200d0,250d0,500d0], 'pt34_1211.1913')
              histos(4) = plot_setup_custom([0d0,0.5d0,1d0,1.5d0,1.75d0,2d0, &
                  2.25d0,2.35d0,2.45d0,2.55d0,2.65d0,2.7d0,2.75d0,2.8d0,2.85d0, &
                  2.9d0,2.95d0,3d0,3.05d0,3.1d0], 'delphi_1211.1913')

          endif

          if (study == 170403839) then
              histos(3) = plot_setup_custom([0d0,4d0,8d0,12d0,16d0,20d0,25d0, &
                  30d0,35d0,40d0,45d0,50d0,55d0,60d0,65d0,70d0,75d0,80d0,90d0, &
                  100d0,110d0,120d0,130d0,140d0,150d0,175d0,200d0,225d0,250d0, &
                  300d0,400d0,750d0], 'pt34_1704.03839')
              histos(4) = plot_setup_custom([0d0,0.25d0,0.5d0,0.75d0,1d0, &
                  1.25d0,1.5d0,1.625d0,1.75d0,1.875d0,2d0,2.125d0,2.25d0, &
                  2.3d0,2.35d0,2.4d0,2.45d0,2.5d0,2.55d0,2.6d0,2.65d0,2.675d0, &
                  2.7d0,2.725d0,2.75d0,2.775d0,2.8d0,2.825d0,2.85d0,2.875d0, &
                  2.9d0,2.925d0,2.95d0,2.975d0,3d0,3.025d0,3.050d0,3.075d0, &
                  3.1d0,3.12d0,3.142d0], 'delphi_1704.03839')
          endif

          histos(6) = plot_setup_custom([0d0,0.004d0,0.008d0,0.012d0,0.016d0, &
              0.02d0,0.024d0,0.029d0,0.034d0,0.039d0,0.045d0,0.051d0,0.057d0, &
              0.064d0,0.072d0,0.081d0,0.091d0,0.102d0,0.114d0,0.128d0,0.145d0, &
              0.165d0,0.189d0,0.219d0,0.258d0,0.312d0,0.391d0,0.524d0,0.695d0, &
              0.918d0,1.153d0,1.496d0,1.947d0,2.522d0,3.277d0,5d0,10d0,20d0,&
              50d0,100d0,50000d0], 'phistar_1704.03839')

      end subroutine

      subroutine book(p,wt,ids,vals,wts)
          use types
          use ResummationTransition
          use, intrinsic :: ieee_arithmetic
          implicit none
          include 'mxpart.f'
          include 'kpart.f'
          include 'taucut.f'! abovecut
          include 'frag.f'

          real(dp), intent(in) :: p(mxpart,4), wt
          integer, allocatable, intent(out) :: ids(:)
          real(dp), allocatable, intent(out) :: vals(:)
          real(dp), allocatable, intent(out) :: wts(:)
          real(dp) :: pttwo,massvec,delphi,etarap
          real(dp) :: pt34com
          common/pt34com/pt34com
!$omp threadprivate(/pt34com/)

          real(dp) :: pt34,m34
          real(dp) :: trans
          real(dp) :: phistar, phiacop, costhetastar, delphi34

          pt34 = pttwo(3,4,p)

          m34 = massvec(p(3,:)+p(4,:))
          if (m34 > 0) then
              m34 = sqrt(m34)
          else
              m34 = -1d0
          endif

          delphi34 = delphi(p(3,:),p(4,:))
          phiacop = 2._dp*atan(sqrt((1._dp+cos(delphi34))/(1._dp-cos(delphi34))))
          costhetastar = tanh((etarap(3,p)-etarap(4,p))/2._dp)
          phistar = tan(phiacop/2._dp)*sin(acos(costhetastar))

          ids = histos

          if (origKpart == kresummed) then
              if (abovecut .eqv. .false.) then
                  trans = transition((pt34/m34)**2d0,0.001d0,transitionSwitch,0.001d0)
              else
                  ! fo piece without transition
                  trans = 1._dp
              endif
          else
              trans = 1._dp
          endif

          if (ieee_is_nan(delphi34) .or. (.not. ieee_is_finite(delphi34))) then
              delphi34 = -1d0
          endif

          if (ieee_is_nan(pt34) .or. (.not. ieee_is_finite(pt34))) then
              pt34 = -1d0
          endif

          if (ieee_is_nan(m34) .or. (.not. ieee_is_finite(m34))) then
              m34 = -1d0
          endif

          if (ieee_is_nan(phistar) .or. (.not. ieee_is_finite(phistar))) then
              phistar = -1d0
          endif

          !DEBUG
          !trans = 1d0

          vals = [pt34,m34,pt34,delphi34,phistar]

          wts = [wt*trans,wt,wt*trans,wt*trans,wt*trans]

      end subroutine

end module
