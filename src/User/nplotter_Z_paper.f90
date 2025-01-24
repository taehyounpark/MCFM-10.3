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

          if ( index(rundir, 'above') > 0) then
              allocate(histos(4))
              resabove = .true.
          else
              allocate(histos(9))
              resabove = .false.
          endif

          if (rank == 0) then
              write (*,*) "RESUMMATION: Using transition with switch ", transitionSwitch
          endif

          qtvals = [1.0000d0,1.2589d0,1.5849d0,1.9953d0, &
                2.5119d0,3.1623d0,3.9811d0,5.0119d0,6.3096d0,7.9433d0, &
                10.0000d0,12.5893d0,15.8489d0,19.9526d0,25.1189d0, &
                31.6228d0,39.8107d0,50.1187d0,63.0957d0,79.4328d0,100.0000d0, &
                125.89d0, 158.49d0, 199.53d0, 251.19d0, 316.23d0, 398.11d0, &
                501.19d0, 630.96d0, 794.33d0]
          if (.not. resabove) then

          histos(1) = plot_setup_custom(qtvals,'pt34_fine_trans04')

          !histos(2) = plot_setup_uniform(20d0,80d0,1d0,'pt3_uni_trans04')

          histos(2) = plot_setup_custom([25d0,26d0,27d0,28d0,29d0,30d0,31d0, &
                    32d0,33d0,34d0,35d0,36d0,37d0,38d0,39d0,40d0,41d0,42d0, &
                    43d0,44d0,45d0,46d0,47d0,48d0,49d0,50d0,52d0,54d0,56d0,58d0,60d0,64d0,68d0, &
                    72d0,76d0,80d0], 'pt3_uninew_trans04')

          if ( index(runstring, 'CMS') > 0) then
              histos(3) = plot_setup_custom([0d0,1d0,2d0,3d0,4d0,5d0,6d0,7d0, &
                  8d0,9d0,10d0,11d0,12d0,13d0,14d0,16d0,18d0,20d0,22d0,25d0, &
                  28d0,32d0,37d0,43d0,52d0,65d0,85d0,120d0,160d0,190d0,220d0, &
                  250d0,300d0,400d0,500d0,800d0,1650d0],'pt34_cms13_trans04')            

              histos(4) = plot_setup_custom([0.001d0,0.002d0,0.003d0,0.004d0,0.005d0, &
                  0.006d0,0.007d0,0.008d0,0.009d0,0.01d0,0.02d0,0.03d0,0.04d0,0.05d0, &
                  0.06d0,0.07d0,0.08d0,0.09d0,0.1d0,0.2d0,0.3d0,0.4d0,0.5d0,0.6d0, &
                  0.7d0,0.8d0,0.9d0,1.0d0,3.0d0,5.0d0,7.0d0,10.0d0,20.0d0,30.0d0,50.0d0], 'phistar_cms13_trans04')
          elseif ( index(runstring, 'ATLAS') > 0 ) then
              histos(3) = plot_setup_custom([0.0d0,2.0d0,4.0d0,6.0d0,8.0d0,10.0d0,12.0d0, &
                  14.0d0,16.0d0,18.0d0,20.0d0,22.5d0,25.0d0,27.5d0,30.0d0,33.0d0,36.0d0, &
                  36.0d0,39.0d0,42.0d0,45.0d0,48.0d0,51.0d0,54.0d0,57.0d0,61.0d0,65.0d0, &
                  70.0d0,75.0d0,80.0d0,85.0d0,95.0d0,105.0d0,125.0d0,150.0d0,175.0d0,200.0d0, &
                  250.0d0,300.0d0,350.0d0,400.0d0,470.0d0,550.0d0,650.0d0,900.0d0,2500.0d0], 'pt34_atlas13_trans04')
              histos(4) = plot_setup_custom([0.0d0,0.004d0,0.008d0,0.012d0,0.016d0,0.02d0, &
                  0.024d0,0.029d0,0.034d0,0.039d0,0.045d0,0.051d0,0.057d0,0.064d0,0.072d0, &
                  0.081d0,0.091d0,0.102d0,0.114d0,0.128d0,0.145d0,0.165d0,0.189d0,0.219d0, &
                  0.258d0,0.312d0,0.391d0,0.524d0,0.695d0,0.918d0,1.153d0,1.496d0,1.947d0, &
                  2.522d0,3.277d0,5.000d0,10.0d0], 'phistar_atlas13_trans04')
          endif



          histos(5) = plot_setup_custom(qtvals,'pt34_fine_trans06')

          !histos(6) = plot_setup_uniform(20d0,80d0,1d0,'pt3_uni_trans06')

          histos(6) = plot_setup_custom([25d0,26d0,27d0,28d0,29d0,30d0,31d0, &
                    32d0,33d0,34d0,35d0,36d0,37d0,38d0,39d0,40d0,41d0,42d0, &
                    43d0,44d0,45d0,46d0,47d0,48d0,49d0,50d0,52d0,54d0,56d0,58d0,60d0,64d0,68d0, &
                    72d0,76d0,80d0], 'pt3_uninew_trans06')

          if ( index(runstring, 'CMS') > 0) then
              histos(7) = plot_setup_custom([0d0,1d0,2d0,3d0,4d0,5d0,6d0,7d0, &
                  8d0,9d0,10d0,11d0,12d0,13d0,14d0,16d0,18d0,20d0,22d0,25d0, &
                  28d0,32d0,37d0,43d0,52d0,65d0,85d0,120d0,160d0,190d0,220d0, &
                  250d0,300d0,400d0,500d0,800d0,1650d0],'pt34_cms13_trans06')            

              histos(8) = plot_setup_custom([0.001d0,0.002d0,0.003d0,0.004d0,0.005d0, &
                  0.006d0,0.007d0,0.008d0,0.009d0,0.01d0,0.02d0,0.03d0,0.04d0,0.05d0, &
                  0.06d0,0.07d0,0.08d0,0.09d0,0.1d0,0.2d0,0.3d0,0.4d0,0.5d0,0.6d0, &
                  0.7d0,0.8d0,0.9d0,1.0d0,3.0d0,5.0d0,7.0d0,10.0d0,20.0d0,30.0d0,50.0d0], 'phistar_cms13_trans06')
          elseif ( index(runstring, 'ATLAS') > 0 ) then
              histos(7) = plot_setup_custom([0.0d0,2.0d0,4.0d0,6.0d0,8.0d0,10.0d0,12.0d0, &
                  14.0d0,16.0d0,18.0d0,20.0d0,22.5d0,25.0d0,27.5d0,30.0d0,33.0d0,36.0d0, &
                  36.0d0,39.0d0,42.0d0,45.0d0,48.0d0,51.0d0,54.0d0,57.0d0,61.0d0,65.0d0, &
                  70.0d0,75.0d0,80.0d0,85.0d0,95.0d0,105.0d0,125.0d0,150.0d0,175.0d0,200.0d0, &
                  250.0d0,300.0d0,350.0d0,400.0d0,470.0d0,550.0d0,650.0d0,900.0d0,2500.0d0], 'pt34_atlas13_trans06')
              histos(8) = plot_setup_custom([0.0d0,0.004d0,0.008d0,0.012d0,0.016d0,0.02d0, &
                  0.024d0,0.029d0,0.034d0,0.039d0,0.045d0,0.051d0,0.057d0,0.064d0,0.072d0, &
                  0.081d0,0.091d0,0.102d0,0.114d0,0.128d0,0.145d0,0.165d0,0.189d0,0.219d0, &
                  0.258d0,0.312d0,0.391d0,0.524d0,0.695d0,0.918d0,1.153d0,1.496d0,1.947d0, &
                  2.522d0,3.277d0,5.000d0,10.0d0], 'phistar_atlas13_trans06')
          endif

              histos(9) = plot_setup_custom([0d0,1d0,2d0,3d0,4d0,5d0,6d0,7d0, &
                  8d0,9d0,10d0,11d0,12d0,13d0,14d0,16d0,18d0,20d0,22d0,25d0, &
                  28d0,32d0,37d0,43d0,52d0,65d0,85d0,120d0,160d0,190d0,220d0, &
                  250d0,300d0,400d0,500d0,800d0,1650d0],'pt34_cms13_notrans')

          else ! above

          histos(1) = plot_setup_custom(qtvals,'pt34_fine')

          histos(2) = plot_setup_uniform(20d0,80d0,1d0,'pt3_uni')

          if ( index(runstring, 'CMS') > 0) then
              histos(3) = plot_setup_custom([0d0,1d0,2d0,3d0,4d0,5d0,6d0,7d0, &
                  8d0,9d0,10d0,11d0,12d0,13d0,14d0,16d0,18d0,20d0,22d0,25d0, &
                  28d0,32d0,37d0,43d0,52d0,65d0,85d0,120d0,160d0,190d0,220d0, &
                  250d0,300d0,400d0,500d0,800d0,1650d0],'pt34_cms13')            

              histos(4) = plot_setup_custom([0.001d0,0.002d0,0.003d0,0.004d0,0.005d0, &
                  0.006d0,0.007d0,0.008d0,0.009d0,0.01d0,0.02d0,0.03d0,0.04d0,0.05d0, &
                  0.06d0,0.07d0,0.08d0,0.09d0,0.1d0,0.2d0,0.3d0,0.4d0,0.5d0,0.6d0, &
                  0.7d0,0.8d0,0.9d0,1.0d0,3.0d0,5.0d0,7.0d0,10.0d0,20.0d0,30.0d0,50.0d0], 'phistar_cms13')
          elseif ( index(runstring, 'ATLAS') > 0 ) then
              histos(3) = plot_setup_custom([0.0d0,2.0d0,4.0d0,6.0d0,8.0d0,10.0d0,12.0d0, &
                  14.0d0,16.0d0,18.0d0,20.0d0,22.5d0,25.0d0,27.5d0,30.0d0,33.0d0,36.0d0, &
                  36.0d0,39.0d0,42.0d0,45.0d0,48.0d0,51.0d0,54.0d0,57.0d0,61.0d0,65.0d0, &
                  70.0d0,75.0d0,80.0d0,85.0d0,95.0d0,105.0d0,125.0d0,150.0d0,175.0d0,200.0d0, &
                  250.0d0,300.0d0,350.0d0,400.0d0,470.0d0,550.0d0,650.0d0,900.0d0,2500.0d0], 'pt34_atlas13')
              histos(4) = plot_setup_custom([0.0d0,0.004d0,0.008d0,0.012d0,0.016d0,0.02d0, &
                  0.024d0,0.029d0,0.034d0,0.039d0,0.045d0,0.051d0,0.057d0,0.064d0,0.072d0, &
                  0.081d0,0.091d0,0.102d0,0.114d0,0.128d0,0.145d0,0.165d0,0.189d0,0.219d0, &
                  0.258d0,0.312d0,0.391d0,0.524d0,0.695d0,0.918d0,1.153d0,1.496d0,1.947d0, &
                  2.522d0,3.277d0,5.000d0,10.0d0], 'phistar_atlas13')
          endif



          endif

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

          pt3 = pt(3,p)
          pt34 = pttwo(3,4,p)

          delphi34 = delphi(p(3,:),p(4,:))
          phiacop = 2._dp*atan(sqrt((1._dp+cos(delphi34))/(1._dp-cos(delphi34))))
          costhetastar = tanh((etarap(3,p)-etarap(4,p))/2._dp)
          phistar = tan(phiacop/2._dp)*sin(acos(costhetastar))


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
          if (.not. resabove) then
              vals = [pt34,pt3,pt34,phistar, pt34, pt3, pt34,phistar, pt34]
              wts = [wt*trans04,wt*trans04,wt*trans04,wt*trans04, &
                     wt*trans06,wt*trans06,wt*trans06,wt*trans06,wt]
          else
              vals = [pt34,pt3,pt34,phistar]
              wts = [wt,wt,wt,wt]
          endif
                ! ,wt*trans*A0, &
                !wt*trans*A1,wt*trans*A2,wt*trans*A3, &
                !wt*trans*A4,wt*trans*A7]

      end subroutine
      
      
end module
