!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
     
    module ResummationIntegration
        use types
        use qtResummation_params
        implicit none

        public :: resint
        public :: resexpint
        public :: qtsubint

        private

        logical, parameter :: boostBorn = .true.
        ! do not set to .true.!
        logical, parameter :: phaseUseQt = .false.

        contains

        function qtsubint(r,wgt)
            use ieee_arithmetic
            use qtResummation
            use SCET
            use qqb_z_hardfun
            use Scalevar
            use MCFMSetupPlots
            implicit none
            include 'src/Inc/mxdim.f'
            include 'src/Inc/mxpart.f'
            include 'src/Inc/npart.f'
            include 'src/Inc/nproc.f'
            include 'src/Inc/kpart.f'
            include 'src/Inc/kprocess.f'
            include 'src/Inc/nf.f'
            include 'src/Inc/energy.f'
            include 'src/Inc/dynamicscale.f'
            include 'src/Inc/constants.f'
            include 'src/Inc/initialscales.f'
            include 'src/Inc/xmin.f'
            include 'src/Inc/scale.f'
            include 'src/Inc/qcdcouple.f'
            include 'src/Inc/ewcharge.f'
            include 'src/Inc/ewcouple.f'
            real(dp) :: qtsubint
            real(dp), intent(in) :: r(mxdim), wgt

            real(dp) :: p(mxpart,4), pjet(mxpart,4), pswt, xx(2)
            real(dp) :: val,flux,qsq
            real(dp) :: msqall(-nf:nf,-nf:nf,0:3)
            real(dp), allocatable :: results(:)
            real(dp), allocatable :: xmsq(:)
            real(dp) :: hard(0:3)
            integer :: j,k, order

            ! functions
            logical :: includedipole
            real(dp) :: twomass

            logical:: bin
            common/bin/bin
            real(dp) :: BrnRat
            common/BrnRat/BrnRat
            include 'beamtype.f'

            qtsubint = 0._dp
            p(:,:) = 0._dp
            pjet(:,:) = 0._dp

            if (any(nproc == [31,32])) then ! Z -> e e
                npart = 2
                ! generate phase space at zero pT
                if (gen_res(r,0._dp,p,pswt) .eqv. .false. ) then
                    qtsubint = 0._dp
                    return
                endif
            else
                write (*,*) __FILE__//": undefined nproc, line ", __LINE__
                error stop 
            endif


            if (any(ieee_is_nan(p(1:npart+2,:)) .eqv. .true.)) then
              !write(6,*) 'Discarding NaN or infinite phase space point'
              qtsubint = 0._dp
              return
            endif

! these are calculated before the boost
            xx(1)=-2._dp*p(1,4)/sqrts
            xx(2)=-2._dp*p(2,4)/sqrts

            if ( (xx(1) > one)  .or. (xx(2) > one) &
             .or.(xx(1) < xmin) .or. (xx(2) < xmin)) then
                qtsubint = 0._dp
                return
            endif

            includeTaucutgrid(0) = .TRUE. 
            if (includedipole(0,p) .eqv. .FALSE. ) then
                qtsubint = 0._dp
                return
            endif

            call storeptilde(0,p)
            call getptildejet(0,pjet)

            if (dynamicscale) then
                call scaleset(initscale, initfacscale, p)
            else
                call usescales(initscale, initfacscale)
            endif

            flux = fbGeV2/(two*xx(1)*xx(2)*sqrts**2)

            if (origkpart == ksnlo) then
                order = 1
            elseif (origkpart == knnlo) then
                order = 2
            elseif (origkpart == kn3lo) then
                order = 3
            else
                write (*,*) __FILE__//": undefined kpart, line ", __LINE__
                error stop
            endif
            
            if (kcase == kZ_only) then
                qsq = twomass(1,2,p)**2

                ! must not contain the axial singlet
                ! contribution since it has no counterpart in the Z+jet
                ! NNLO calculation
                call qqb_z_loop(p, musq, msqall, withaxsinglet=.false.)

                msqall(:,:,1) = msqall(:,:,1) * (as/4._dp/pi)
                msqall(:,:,2) = msqall(:,:,2) * (as/4._dp/pi)**2
                msqall(:,:,3) = msqall(:,:,3) * (as/4._dp/pi)**3
            else
                write (*,*) __FILE__//": undefined kcase, line ", __LINE__
                error stop 
            endif

            allocate(xmsq(size(tcutarray)+1), source=0._dp)

            do j=-5,5; do k =-5,5
                ! cycle if no contribution through NLO
                if ( all(msqall(j,k,[0,1]) == 0._dp) ) cycle

                ! coefficient in as/4/pi
                hard(0) = msqall(j,k,0) 
                hard(1) = msqall(j,k,1) / msqall(j,k,0) / (as/4._dp/pi)
                hard(2) = msqall(j,k,2) / msqall(j,k,0) / (as/4._dp/pi)**2
                hard(3) = msqall(j,k,2) / msqall(j,k,0) / (as/4._dp/pi)**2

                call qtsubtraction(p, qsq, xx(1), xx(2), j, k, hard, & 
                    sqrt(musq), order, 1._dp, results)

                xmsq = xmsq + results

            enddo; enddo


            if (doMultitaucut) then
                scetreweight(:) = 0._dp
                if (xmsq(1) /= 0._dp) then
                    scetreweight(:) = xmsq(2:size(tcutarray)+1) / xmsq(1)
                endif
            endif

            qtsubint = flux*pswt*xmsq(1)/BrnRat

            val=qtsubint*wgt

            if (ieee_is_nan(val) .OR. ( .NOT. ieee_is_finite(val))) then
                !write(6,*) 'Discarded NaN, val=',val
                qtsubint = 0._dp
                return
            endif
                  
            if (bin) then
                includeTaucutgrid(0) = .true. 
                call nplotter_new(pjet,val)
            endif

        end function qtsubint

        function resint(r,wgt)
            use ieee_arithmetic
            use types
            use Scalevar
            use PDFerrors
            use SCET
            use MCFMStorage
            use qtResummation
            use LHAPDF
            use MCFMSetupPlots
            use SafetyCuts, only : passed_smallnew
            use ptveto, only: usept, jetptveto
            implicit none
            include 'constants.f'
            include 'mxpart.f'
            include 'vegas_common.f'
            include 'kprocess.f'
            include 'kpart.f'
            include 'energy.f'
            include 'npart.f'
            include 'dynamicscale.f'
            include 'initialscales.f'
            include 'scalevar.f'
            include 'scale.f'
            include 'facscale.f'
            include 'qcdcouple.f'
            include 'couple.f'
            include 'nlooprun.f'
            include 'x1x2.f'
            include 'xmin.f'
            include 'nf.f'
            include 'nflav.f'
            include 'nwz.f'
            include 'nproc.f'
            include 'ipsgen.f'
            include 'masses.f'
            include 'beamtype.f'
            logical:: bin
            common/bin/bin
            real(dp) :: BrnRat
            common/BrnRat/BrnRat

            real(dp) :: resint
            real(dp), intent(in) :: r(mxdim), wgt

            real(dp):: p(mxpart,4),pjet(mxpart,4),xmsq,xmsqMu,xmsqPdf
            real(dp) :: val,pswt,flux

            real(dp) :: qsq
            real(dp) :: p_born(mxpart,4)

            ! functions
            logical :: includedipole
            real(dp) :: ptpure, twomass

            real(dp) :: xjac
            real(dp) :: qt, phi, pttwo
            real(dp) :: hard(2)
            real(dp) :: res, resexp
            real(dp) :: msq(-nf:nf, -nf:nf)
            real(dp) :: phaseqt
            real(dp) :: safety

            integer :: j,k

            resint = 0._dp

            currentPDF = 0
            currentNd = 0

            if (doScalevar .and. bin) then
                scalereweight(:) = 1._dp
            endif

            p(:,:) = 0._dp
            pjet(:,:) = 0._dp

            if (qtmaxRes .eq. qtminRes) then
                qt = qtminRes
                xjac = 2._dp*qt
            else
                qt = r(ndim)*(qtmaxRes-qtminRes) + qtminRes
                xjac = 2._dp*qt*(qtmaxRes-qtminRes)
            endif
            
            phi = 2._dp*pi*r(ndim-1)

    !       ! for debugging purposes we can fix qt here
            !qt = 8.0_dp
            !xjac = 2._dp*qt

            if (qt > qtmaxRes) then
                resint = 0._dp
                return
            endif

            if (qt < qtminRes) then
                resint = 0._dp
                return
            endif

            if (usept) then
              qt = jetptveto
              xjac=1._dp
            endif

            p(:,:) = 0._dp

            if (phaseUseQt) then
                phaseqt = qt
            else
                phaseqt = 0._dp
            endif

            safety = 1d-9

            if (nproc == 1 .or. nproc == 6) then ! W+, W-
                npart = 2
                if (gen_res(r,phaseqt,p,pswt) .eqv. .false. ) then
                    resint = 0._dp
                endif
            elseif (any(nproc == [31,32])) then ! Z -> e e
                npart = 2
                ! generate phase space at zero pT
                if (gen_res(r,phaseqt,p,pswt) .eqv. .false. ) then
                    resint = 0._dp
                endif
            elseif (kcase == kggfus0) then ! H -> ga ga
                npart = 2
                !if (gen_res(r,phaseqt,p,pswt) .eqv. .false.) then
                    !resint = 0._dp
                    !return
                !endif
                call gen2m(r,p,pswt,*999)
            elseif (kcase == kHi_Zga) then
                npart=3
                call gen3h(r,p,pswt,*999)
            elseif (kcase == kWHbbar .or. kcase == kWHgaga) then
                npart=4
                call genVH(r,p,pswt,*999)
            elseif (kcase == kZHbbar .or. kcase == kZHgaga) then
                npart=4
                call genVH(r,p,pswt,*999)
            elseif (kcase == kZH__WW) then
                npart=6
                call gen6(r,p,pswt,*999)
            elseif (any(nproc == [290,295])) then ! W -> e nu
                npart = 3
                if (ipsgen == 1) then
                    call gen_Vphotons_jets(r,1,0,p,pswt,*999)
                elseif (ipsgen == 2) then
                    call gen_Vphotons_jets_dkrad(r,1,0,p,pswt,*999)
                else
                  error stop
                endif
                wgaips: block
                    real(dp) :: s34, s345, wt34, wt345, wtips(2)
                    real(dp) :: dot
                  s34 = 2._dp*dot(p,3,4)
                  s345 = s34+2._dp*dot(p,3,5)+2._dp*dot(p,4,5)
                  wt34 = (s34-wmass**2)**2+(wmass*wwidth)**2
                  wt345 = (s345-wmass**2)**2+(wmass*wwidth)**2
                  wtips(1) = wt345
                  wtips(2) = wt34
                  pswt = pswt*wtips(ipsgen)/(wtips(1)+wtips(2))
                end block wgaips
            elseif (any(nproc == [300,305])) then ! Z -> e+ e-
                npart = 3
                if (ipsgen == 1) then
                    call gen_Vphotons_jets(r,1,0,p,pswt,*999)
                elseif (ipsgen == 2) then
                    call genVphoton(r,p,pswt,*999)
                else
                  error stop
                endif
                zgaips: block
                    real(dp) :: s34, s345, wt34, wt345, wtips(2)
                    real(dp) :: dot
                  s34 = 2._dp*dot(p,3,4)
                  s345 = s34+2._dp*dot(p,3,5)+2._dp*dot(p,4,5)
                  wt34 = s34*(s34-zmass**2)**2+(zmass*zwidth)**2
                  wt345 = s345*(s345-zmass**2)**2+(zmass*zwidth)**2
                  wtips(1) = wt345
                  wtips(2) = wt34
                  pswt = pswt*wtips(ipsgen)/(wtips(1)+wtips(2))
                end block zgaips
            elseif (nproc == 285 .or. nproc == 2851 .or. nproc == 2852) then ! gamma gamma
                npart = 2
                call gen_photons_jets_res(r,2,0,p,pswt,*999)
                safety = 1d-7
            elseif ( ((kcase == kWWqqbr) .or. (kcase == kZZlept)) &
                .and. (ipsgen == 2)) then
                npart=4
                call gen4handc(r,p,pswt,*999)
            elseif ((kcase == kWWqqbr) .or. (kcase == kWZbbar) .or. (kcase == kZZlept)) then
                npart = 4
                call gen4(r,p,pswt,*999)
            else
                write (*,*) __FILE__//": undefined nproc, line ", __LINE__
                error stop 
            endif

! back-to-back check
            if ((kcase == kWWqqbr) .or. (kcase == kWZbbar) .or. (kcase == kZZlept)) then
                if (pttwo(3,4,p) < 1.e-3_dp) goto 999
            endif

            if (.not. passed_smallnew(p,npart,safety)) then
     999        continue
                resint = 0._dp
                return
            endif


            if (any(ieee_is_nan(p(1:npart,:)) .eqv. .true.)) then
              !write(6,*) 'Discarding NaN or infinite phase space point'
              resint = 0._dp
              return
            endif


    ! these are calculated before the boost
            xx(1)=-2._dp*p(1,4)/sqrts
            xx(2)=-2._dp*p(2,4)/sqrts

            if ( (xx(1) > one)  .or. (xx(2) > one) &
             .or.(xx(1) < xmin) .or. (xx(2) < xmin)) then
                resint = 0._dp
                return
            endif

    ! save before boost for matrix element evaluation
            p_born(1:npart+2,:) = p(1:npart+2,:)
    ! mostly for plotting purposes
            call recoilBoost(qt, phi, p)

! Removed this: it only seems to spoil p_born
!            if (boostBorn) then
!                p_born = p
!            endif

    ! check cuts, fills ptildejet
            includeTaucutgrid(0) = .TRUE. 
            if (includedipole(0,p) .eqv. .FALSE. ) then
                resint = 0._dp
                return
            endif

    ! fill pjet array
            call getptildejet(0,pjet)


            flux = fbGeV2/(two*xx(1)*xx(2)*sqrts**2)

            ! central value
            call computeResint(r, p, p_born, qt, kresorder, xx, resout=xmsq)

            if (doScalevar .and. bin) then
                do j=1,maxscalevar
                    call computeResint(r,p,p_born,qt,kresorder,xx, resout=xmsqMu, &
                        muMultHard_in=scalevarmult(j), muMultRes_in=facscalevarmult(j))
                    scalereweight(j) = xmsqMu / xmsq
                enddo

                if (scalevar_rapidity) then
                    do j=1,2
                        scalevar_rapidity_i = j
                        call computeResint(r,p,p_born,qt,kresorder,xx,resout=xmsqMu)
                        scalereweight(maxscalevar+j) = xmsqMu / xmsq
                        if (ieee_is_nan(scalereweight(maxscalevar+j))) then
                            write (*,*) "NaN in scalevar_rapidity", xmsqMu, xmsq
                        endif
                    enddo
                    scalevar_rapidity_i = 0
                endif
            endif

            if ((maxPDFsets > 0) .and. bin) then
                pdfreweight(:) = 0._dp
                do j=1,maxPDFsets
                    currentPDF = j
                    call computeResint(r, p, p_born, qt, kresorder, xx, resout=xmsqPdf)
                    pdfreweight(currentPDF) = flux*xjac*pswt*(xmsq-xmsqPdf)/BrnRat*wgt
                enddo
            endif

            resint = flux*xjac*pswt*xmsq/BrnRat

            !write (*,*) "BrnRat = ", BrnRat
            !write (*,*) "resint = ", resint
            !pause

            val=resint*wgt

            if (ieee_is_nan(val) .OR. ( .NOT. ieee_is_finite(val))) then
                write(6,*) 'Discarded NaN, val=',val
                resint = 0._dp
                return
            endif

            if (bin) then
                includeTaucutgrid(0) = .true. 
                call nplotter_new(pjet,val)
            endif

        end function resint

        subroutine computeResint(r, p, p_born, qt, order, xx, resout, resexpout, muMultHard_in, muMultRes_in)
            use types
            use qtResummation
            use Scalevar
            use GammaGamma3L
            use qqb_z_hardfun
            use ptveto
            implicit none
            include 'nf.f'
            include 'mxpart.f'
            include 'kprocess.f'
            include 'constants.f'
            include 'initialscales.f'
            include 'dynamicscale.f'
            include 'scale.f'
            include 'qcdcouple.f'
            include 'ewcharge.f'
            include 'ewcouple.f'
            include 'ipsgen.f'
            include 'newpoint.f'
            include 'Rcut.f'
            include 'taucut.f'
            include 'vegas_common.f'
            include 'nflav.f'

            real(dp), intent(in) :: r(mxdim), p(mxpart,4), p_born(mxpart,4)
            real(dp), intent(in) :: qt
            integer, intent(in) :: order
            real(dp), intent(in) :: xx(2)
            real(dp), intent(out), optional :: resout, resexpout
            real(dp), intent(in), optional :: muMultHard_in, muMultRes_in

            real(dp) :: twomass, threemass, fourmass, puremass
            real(dp) :: gg_2gam_msbar

            integer :: j,k
            real(dp) :: muMultHard, muMultRes, hard(2), hard3(3), Ct2(2), f0sq
            real(dp) :: xmsq, qsq, msq(-nf:nf,-nf:nf), res, resexp, xmsq1, xmsq2

            real(dp) :: msqall(-nf:nf,-nf:nf,0:3), ggcontrib
            real(dp) :: mtcorr, NfV

            real(dp) :: msq0(-nf:nf,-nf:nf),msq1(-nf:nf,-nf:nf),msq2(-nf:nf,-nf:nf)
            real(dp) :: msqall_tmp(-nf:nf,-nf:nf,0:3)

            integer :: hardorder

            if (present(muMultHard_in)) then
                muMultHard = muMultHard_in
            else
                muMultHard = 1._dp
            endif

            if (present(muMultRes_in)) then
                muMultRes = muMultRes_in
            else
                muMultRes = 1._dp
            endif

! For the pure-resummed pt-veto calculation set special central scales
            if (usept .and. present(resout)) then
! note that, for this suite of routines
!  mu = facscale, central value = ptveto (qt)
! muh = scale,  central values = sqrt(qsq)
              qsq = twomass(1,2,p)**2
              call usescales(sqrt(qsq)*initscale*muMultHard, qt*initfacscale*muMultres)
            else
              if (dynamicscale) then
                  call scaleset(initscale*muMultHard, initfacscale*muMultRes, p)
              else
                  call usescales(initscale*muMultHard, initfacscale*muMultRes)
              endif
            endif

            msqall = 0._dp

            if (order < 3) then
               hardorder = 0
            elseif (order < 5) then
               hardorder = 1
            elseif (order < 7) then
                hardorder = 2
            else
               hardorder = 3
            endif

            if (kcase == kW_only) then
                qsq = twomass(1,2,p)**2
                call hardqq(qsq, musq, hard)
                if (boostBorn) then
                    call qqb_w(p, msq)
                else
                    call qqb_w(p_born, msq)
                endif

                where(msq /= 0._dp) msqall(:,:,0) = msq(:,:)
                where(msq /= 0._dp) msqall(:,:,1) = hard(1)*msq(:,:)*(as/2._dp/pi)
                where(msq /= 0._dp) msqall(:,:,2) = hard(2)*msq(:,:)*(as/2._dp/pi)**2

            elseif (.false. .and. kcase == kZ_only) then
                qsq = twomass(1,2,p)**2
                call hardqq(qsq, musq, hard)
                if (boostBorn) then
                    call qqb_z(p, msq)
                else
                    call qqb_z(p_born, msq)
                endif

                where(msq /= 0._dp) msqall(:,:,0) = msq(:,:)
                where(msq /= 0._dp) msqall(:,:,1) = hard(1)*msq(:,:)*(as/2._dp/pi)
                where(msq /= 0._dp) msqall(:,:,2) = hard(2)*msq(:,:)*(as/2._dp/pi)**2
                call qqb_z_loop(p, musq, msqall_tmp)

                if (order >= 7) then
                    msqall(:,:,3) = 0._dp
                    call hardqq3(Qsq,musq,NfV,hard3)

                    do j=-nf,nf
                        if (j==0) cycle
                        k=-j
                        msqall(j,k,3) = hard3(3)*msq(j,k)*(as/2._dp/pi)**3 
                    enddo
                endif

                 !for debugging
                do j=-nf,nf
                write (*,*) "tree", j,-j, msqall(j,-j,0) / msqall_tmp(j,-j,0)
                write (*,*) "1loop", j,-j, msqall(j,-j,1) / (msqall_tmp(j,-j,1) * (as/4._dp/pi))
                write (*,*) "2loop", j,-j, msqall(j,-j,2) / (msqall_tmp(j,-j,2) * (as/4._dp/pi)**2)
                write (*,*) "3loop", j,-j, msqall(j,-j,3) / (msqall_tmp(j,-j,3) * (as/4._dp/pi)**3)
                write (*,*) "3loop fraction", j,-j, msqall_tmp(j,-j,3)*(as/4._dp/pi)**3 / sum(msqall(j,-j,1:2))
                enddo
                pause

            elseif (kcase == kZ_only) then
                qsq = twomass(1,2,p)**2

                if (present(resexpout)) then
                    ! for resexp hard(2) must not contain the axial singlet
                    ! contribution since it has no counterpart in the Z+jet
                    ! NNLO calculation
                    call qqb_z_loop(p, musq, msqall, withaxsinglet=.false.)
                else
                    ! in the resummation we include axial singlet contributions
                    call qqb_z_loop(p, musq, msqall, withaxsinglet=.true.)
                endif

                msqall(:,:,1) = msqall(:,:,1) * (as/4._dp/pi)
                msqall(:,:,2) = msqall(:,:,2) * (as/4._dp/pi)**2
                msqall(:,:,3) = msqall(:,:,3) * (as/4._dp/pi)**3

            elseif (kcase == kggfus0) then
                if (boostBorn) then
                    call gg_h(p,msq)
                else
                    call gg_h(p_born, msq)
                endif
                qsq = twomass(1,2,p)**2
                if (gghsinglestep) then
                  call CggHsq_onestep(qsq,musq,f0sq,hard)
! apply finite-mt correction to LO (HEFT) result
                  msq(:,:) = msq(:,:)*f0sq
                else
                  call CSsq(hardorder, qsq, musq, hard)
                  call Ctsq(hardorder, musq, Ct2)
! apply either strict fixed-order expansion  ...
                  hard(2)=hard(2)+Ct2(2)+hard(1)*Ct2(1)
                  hard(1)=hard(1)+Ct2(1)
! ... or just multiply, including some higher-order termss
!                  hard(2)=(1d0+hard(1)+hard(2))*(1d0+Ct2(1)+Ct2(2))-(1d0+hard(1))*(1d0+Ct2(1))
!                  hard(1)=(1d0+hard(1))*(1d0+Ct2(1))-1d0
                endif

                where(msq /= 0._dp) msqall(:,:,0) = msq(:,:)
                where(msq /= 0._dp) msqall(:,:,1) = hard(1)*msq(:,:)
                where(msq /= 0._dp) msqall(:,:,2) = hard(2)*msq(:,:)

!!!!---- new two-step/one-step comparison
                if (1 == 2) then
                call Gammafill(nflav)
! fix scales to be the same to demonstrate equivalence
!                call usescales(sqrt(qsq)*initscale*muMultHard, sqrt(qsq)*initfacscale*muMultres)
!                write(6,*) 'as,ason4pi',as,ason4pi

                gghsinglestep=.false.
                call gg_h(p,msq)
                if (gghsinglestep) then
                  call CggHsq_onestep(qsq,musq,f0sq,hard)
! apply finite-mt correction to LO (HEFT) result
                  msq(:,:) = msq(:,:)*f0sq
                else
                  call CSsq(hardorder, qsq, musq, hard)
                  call Ctsq(hardorder, musq, Ct2)
! apply either strict fixed-order expansion  ...
                  hard(2)=hard(2)+Ct2(2)+hard(1)*Ct2(1)
                  hard(1)=hard(1)+Ct2(1)
! ... or just multiply, including some higher-order termss
!                  hard(2)=(1d0+hard(1)+hard(2))*(1d0+Ct2(1)+Ct2(2))-(1d0+hard(1))*(1d0+Ct2(1))
!                  hard(1)=(1d0+hard(1))*(1d0+Ct2(1))-1d0
                endif
                where(msq /= 0._dp) msqall(:,:,0) = msq(:,:)
                where(msq /= 0._dp) msqall(:,:,1) = hard(1)*msq(:,:)
                where(msq /= 0._dp) msqall(:,:,2) = hard(2)*msq(:,:)
                msq(0,0) = msqall(0,0,0)
                call resummation_BNRptlumxmsq(p,xx,r(ndim-3),r(ndim-2),qt,sqrt(qsq),Rcut,msq,1,xmsq)
                write(6,*) 'two-step order=0 xmsq',xmsq
                msq(0,0) = sum(msqall(0,0,[0,1]))
                call resummation_BNRptlumxmsq(p,xx,r(ndim-3),r(ndim-2),qt,sqrt(qsq),Rcut,msq,1,xmsq)
                write(6,*) 'two-step order=1 xmsq',xmsq
                xmsq1=xmsq
                msq(0,0) = sum(msqall(0,0,[0,1,2]))
                call resummation_BNRptlumxmsq(p,xx,r(ndim-3),r(ndim-2),qt,sqrt(qsq),Rcut,msq,2,xmsq)
                write(6,*) 'two-step order=2 xmsq',xmsq
                write(6,*) 'two-step hard(1), hard(2)', hard(1),hard(2)
                xmsq2=xmsq

                gghsinglestep=.true.
                call gg_h(p,msq)
                if (gghsinglestep) then
                  call CggHsq_onestep(qsq,musq,f0sq,hard)
! apply finite-mt correction to LO (HEFT) result
                  msq(:,:) = msq(:,:)*f0sq
                  write(6,*) 'f0sq=',f0sq
                else
                  call CSsq(hardorder, qsq, musq, hard)
                  call Ctsq(hardorder, musq, Ct2)
! apply either strict fixed-order expansion  ...
                  hard(2)=hard(2)+Ct2(2)+hard(1)*Ct2(1)
                  hard(1)=hard(1)+Ct2(1)
! ... or just multiply, including some higher-order termss
!                  hard(2)=(1d0+hard(1)+hard(2))*(1d0+Ct2(1)+Ct2(2))-(1d0+hard(1))*(1d0+Ct2(1))
!                  hard(1)=(1d0+hard(1))*(1d0+Ct2(1))-1d0
                endif
                where(msq /= 0._dp) msqall(:,:,0) = msq(:,:)
                where(msq /= 0._dp) msqall(:,:,1) = hard(1)*msq(:,:)
                where(msq /= 0._dp) msqall(:,:,2) = hard(2)*msq(:,:)
                msq(0,0) = msqall(0,0,0)
                call resummation_BNRptlumxmsq(p,xx,r(ndim-3),r(ndim-2),qt,sqrt(qsq),Rcut,msq,1,xmsq)
                write(6,*) 'one-step order=0 xmsq',xmsq
                msq(0,0) = sum(msqall(0,0,[0,1]))
                call resummation_BNRptlumxmsq(p,xx,r(ndim-3),r(ndim-2),qt,sqrt(qsq),Rcut,msq,1,xmsq)
                write(6,*) 'one-step order=1 xmsq',xmsq,xmsq/xmsq1
                msq(0,0) = sum(msqall(0,0,[0,1,2]))
                call resummation_BNRptlumxmsq(p,xx,r(ndim-3),r(ndim-2),qt,sqrt(qsq),Rcut,msq,2,xmsq)
                write(6,*) 'one-step order=2 xmsq',xmsq,xmsq/xmsq2
                write(6,*) 'one-step hard(1), hard(2)', hard(1),hard(2)
                pause

                endif

!--- debug
!                if (hardorder == 0) then
!                    msq(0,0) = msqall(0,0,0)
!                elseif (hardorder == 1) then
!                    msq(0,0) = sum(msqall(0,0,[0,1]))
!                elseif (hardorder == 2) then
!                    msq(0,0) = sum(msqall(0,0,[0,1,2]))
!                elseif (hardorder == 3) then
!                    msq(0,0) = sum(msqall(0,0,[0,1,2,3]))
!                endif
!                msq(:,:)=msq(:,:)/msqall(:,:,0)
!                call resummation_BNRptlumxmsq(p,xx,r(ndim-3),r(ndim-2),qt,sqrt(qsq),Rcut,msq,hardorder,xmsq)
!--- debug

            elseif (kcase == kHi_Zga) then
                qsq = threemass(3,4,5,p)**2
                call hardgg(qsq, musq, hard)
                if (boostBorn) then
                    call gg_hzgam(p, msq)
                else
                    call gg_hzgam(p_born, msq)
                endif

                where(msq /= 0._dp) msqall(:,:,0) = msq(:,:)
                where(msq /= 0._dp) msqall(:,:,1) = hard(1)*msq(:,:)*(as/2._dp/pi)
                where(msq /= 0._dp) msqall(:,:,2) = hard(2)*msq(:,:)*(as/2._dp/pi)**2
            elseif (kcase == kWHbbar) then
                qsq = fourmass(3,4,5,6,p)**2
                call hardqq(qsq, musq, hard)
                if (boostBorn) then
                    call qqb_wh(p, msq)
                else
                    call qqb_wh(p_born, msq)
                endif

                where(msq /= 0._dp) msqall(:,:,0) = msq(:,:)
                where(msq /= 0._dp) msqall(:,:,1) = hard(1)*msq(:,:)*(as/2._dp/pi)
                where(msq /= 0._dp) msqall(:,:,2) = hard(2)*msq(:,:)*(as/2._dp/pi)**2
            elseif (kcase == kWHgaga) then
                qsq = fourmass(3,4,5,6,p)**2
                call hardqq(qsq, musq, hard)
                if (boostBorn) then
                    call qqb_wh_gaga(p, msq)
                else
                    call qqb_wh_gaga(p_born, msq)
                endif

                where(msq /= 0._dp) msqall(:,:,0) = msq(:,:)
                where(msq /= 0._dp) msqall(:,:,1) = hard(1)*msq(:,:)*(as/2._dp/pi)
                where(msq /= 0._dp) msqall(:,:,2) = hard(2)*msq(:,:)*(as/2._dp/pi)**2
            elseif (kcase == kZHbbar) then
                qsq = fourmass(3,4,5,6,p)**2
                call hardqq(qsq, musq, hard)
                if (boostBorn) then
                    call qqb_zh(p, msq)
                else
                    call qqb_zh(p_born, msq)
                endif

                where(msq /= 0._dp) msqall(:,:,0) = msq(:,:)
                where(msq /= 0._dp) msqall(:,:,1) = hard(1)*msq(:,:)*(as/2._dp/pi)
                where(msq /= 0._dp) msqall(:,:,2) = hard(2)*msq(:,:)*(as/2._dp/pi)**2
            elseif (kcase == kZH__WW) then
                qsq = puremass(p(3,:)+p(4,:)+p(5,:)+p(6,:)+p(7,:)+p(8,:))**2
                call hardqq(qsq, musq, hard)
                if (boostBorn) then
                    call qqb_zh_ww(p, msq)
                else
                    call qqb_zh_ww(p_born, msq)
                endif

                where(msq /= 0._dp) msqall(:,:,0) = msq(:,:)
                where(msq /= 0._dp) msqall(:,:,1) = hard(1)*msq(:,:)*(as/2._dp/pi)
                where(msq /= 0._dp) msqall(:,:,2) = hard(2)*msq(:,:)*(as/2._dp/pi)**2
            elseif (kcase == kZHgaga) then
                qsq = fourmass(3,4,5,6,p)**2
                call hardqq(qsq, musq, hard)
                if (boostBorn) then
                    call qqb_zh_gaga(p, msq)
                else
                    call qqb_zh_gaga(p_born, msq)
                endif

                where(msq /= 0._dp) msqall(:,:,0) = msq(:,:)
                where(msq /= 0._dp) msqall(:,:,1) = hard(1)*msq(:,:)*(as/2._dp/pi)
                where(msq /= 0._dp) msqall(:,:,2) = hard(2)*msq(:,:)*(as/2._dp/pi)**2
            elseif (kcase == kWgamma) then
                qsq = threemass(3,4,5,p)**2
                if (boostBorn) then
                    call wgam_mat(p,msqall)
                else
                    call wgam_mat(p_born,msqall)
                endif
                msq(:,:) = msqall(:,:,0)
            elseif (kcase == kZgamma) then
                qsq = threemass(3,4,5,p)**2
                if (boostBorn) then
                    call zgam_mat(p,msqall)
                    call gg_zgam(p,ggcontrib)
                else
                    call zgam_mat(p_born,msqall)
                    call gg_zgam(p_born,ggcontrib)
                endif

                msqall(0,0, 2) = ggcontrib
                msq(:,:) = msqall(:,:,0)
            elseif (kcase == kgg2gam) then
                qsq = twomass(3,4,p)**2
                if (boostBorn) then
                    call gg_gamgam(p,ggcontrib)
                    !call gggaga_mt(p,ggcontrib) ! with mt effect
                    msqall(0,0, 0) = ggcontrib
                    if (order >= 3) then
                        msqall(0,0, 1) = gg_2gam_msbar(p)
                    endif
                else
                    call gg_gamgam(p_born,ggcontrib)
                    !call gggaga_mt(p_born,ggcontrib) ! with mt effect
                    msqall(0,0, 0) = ggcontrib
                    if (order >= 3) then
                        msqall(0,0, 1) = gg_2gam_msbar(p_born)
                    endif
                endif

            elseif (kcase == kgamgam) then
                qsq = twomass(3,4,p)**2
                if (boostBorn) then
                    call qqb_gamgam(p,msq)
                else
                    call qqb_gamgam(p_born,msq)
                endif

                where(msq /= 0._dp) msqall(:,:,0) = msq(:,:)

                gamgam: block
                    real(dp) :: qqb
                    real(dp) :: qsum, q2sum, fac
                    integer :: jj

                    real(dp) :: twoloop_nfv2, threeloop_nfv, threeloop_nfv2
                    real(dp) :: hard3(3), tree

                    hard3(:) = 0._dp
                    
                    if (boostBorn) then
                        if (order >= 7) then
                            call gamgamhard(p,tree,hard3,twoloop_nfv2, &
                                            threeloop_nfv, threeloop_nfv2)
                            ! normalize as old routines
                            hard3(:) = hard3(:) / tree

                            ! for debugging:
                            !call gamgamampsq_new(2,p,1,2,3,4,qqb,hard,tlrem)
                        else
                            call gamgamampsq_new(2,p,1,2,3,4,qqb,hard3(1:2),twoloop_nfv2)
                        endif
                    else
                        if (order >= 7) then
                            call gamgamhard(p_born,tree,hard3,twoloop_nfv2, &
                                            threeloop_nfv, threeloop_nfv2)
                            ! normalize as old routines
                            hard3(:) = hard3(:) / tree
                        else
                            call gamgamampsq_new(2,p_born,1,2,3,4,qqb,hard3(1:2),twoloop_nfv2)
                        endif
                    endif

                    !write (*,*) "Comparison ", order
                    !write (*,*) "ONELOOP", hard3(1) / (hard(1))
                    !write (*,*) "DIFF", (hard3(1) - hard(1))
                    !write (*,*) "TWOLOOP", hard3(2) / (hard(2))
                    !write (*,*) "TWOLOOP NFV2", twoloop_nfv2 / tlrem
                    !pause

                    fac = xn*esq**2*0.5_dp
                    q2sum = sum(Q(1:5)**2)
                    qsum = sum(Q(1:5))

                    where(msq /= 0._dp) msqall(:,:,1) = hard3(1)*msq(:,:)*(as/2._dp/pi)

                    if (order >= 5) then
                        where(msq /= 0._dp) msqall(:,:,2) = hard3(2)*msq(:,:)*(as/2._dp/pi)**2

                        do jj=-nf,nf
                            if (jj == 0) cycle
                            msqall(jj,-jj,2) = msqall(jj,-jj,2) + &
                                fac*aveqq*twoloop_nfv2 * &
                                q2sum*Q(abs(jj))**2 * (as/2._dp/pi)**2
                        enddo
                    endif

                    if (order >= 7) then
                        where(msq /= 0._dp) msqall(:,:,3) = hard3(3)*msq(:,:)*(as/2._dp/pi)**3

                        do jj=-nf,nf
                            if (jj == 0) cycle

                            msqall(jj,-jj,3) = msqall(jj,-jj,3) + &
                                fac*aveqq*threeloop_nfv2 * &
                                q2sum*Q(abs(jj))**2 * (as/2._dp/pi)**3

                            msqall(jj,-jj,3) = msqall(jj,-jj,3) + &
                                fac*aveqq*threeloop_nfv * &
                                qsum*Q(abs(jj))**2 * (as/2._dp/pi)**3
                        enddo
                endif

            end block gamgam

            !! gg contributions

            !ggcontrib = 0._dp
            !if (boostBorn) then
                !if (order >= 5) then
                    !call gg_gamgam(p,ggcontrib)
                    !!call gggaga_mt(p,ggcontrib) ! with mt effect
                    !msqall(0,0, 2) = ggcontrib
                !endif
                !if (order >= 7) then
                    !msqall(0,0, 3) = gg_2gam_msbar(p)
                !endif
            !else
                !if (order >= 5) then
                    !call gg_gamgam(p_born,ggcontrib)
                    !!call gggaga_mt(p_born,ggcontrib) ! with mt effect
                !endif
                !if (order >= 7) then
                    !msqall(0,0, 3) = gg_2gam_msbar(p_born)
                !endif
            !endif
            
        elseif ((kcase == kWWqqbr) .or. (kcase == kWZbbar) .or. (kcase == kZZlept)) then
            newpoint = .true.
            qsq = fourmass(3,4,5,6,p)**2
!            if (boostBorn) then
!              call hard_VV(p, 0, msq, msq1, msq2)
!            else
!              call hard_VV(p_born, 0, msq, msq1, msq2)
!            endif
! Should use p_born in call below since hard_VV checks pt(34)
#ifdef WITH_VVAMP
            if (present(resexpout) .and. (ipsgen == 1)) then
! hard(2) is never needed but expensive to compute
              call hard_VV(p_born, 1, msq0, msq1, msq2)
            else
              call hard_VV(p_born, hardorder, msq0, msq1, msq2)
            endif
#else
            error stop "Please recompile with VVamp support to run WW/WZ/ZZ resummation"
#endif

            msqall(:,:,0) = msq0(:,:)
            msqall(:,:,1) = msq1(:,:)*(as/4._dp/pi)
            msqall(:,:,2) = msq2(:,:)*(as/4._dp/pi)**2

!--- debug
!                do j=1,2
!                if (hardorder == 0) then
!                    msq(j,-j) = msqall(j,-j,0)
!                elseif (hardorder == 1) then
!                    msq(j,-j) = sum(msqall(j,-j,[0,1]))
!                elseif (hardorder == 2) then
!                    msq(j,-j) = sum(msqall(j,-j,[0,1,2]))
!                elseif (hardorder == 3) then
!                    msq(j,-j) = sum(msqall(j,-j,[0,1,2,3]))
!                endif
!                enddo
!                msq(:,:)=msq(:,:)/msqall(:,:,0)
!                call resummation_BNRptlumxmsq(p,xx,r(ndim-3),r(ndim-2),qt,sqrt(qsq),Rcut,msq,hardorder,xmsq)
!--- debug

!            where(msq /= 0._dp) msqall(:,:,0) = msq(:,:)
!            where(msq /= 0._dp) msqall(:,:,1) = msq(:,:)*msq1(:,:)/msq0(:,:)*(as/4._dp/pi)
!            where(msq /= 0._dp) msqall(:,:,2) = msq(:,:)*msq2(:,:)/msq0(:,:)*(as/4._dp/pi)**2

!            if (ipsgen == 2) msqall(0,0, 2) = msq2(0,0)     ! gg contribution
        else
            write (*,*) __FILE__//": undefined kcase, line", __LINE__
            error stop
        endif



        if (present(resout)) then
! jet-veto resummation
            if (usept) then
                do j=-nflav,nflav; do k =-nflav,nflav
                    if ( all(msqall(j,k,0:hardorder) == 0._dp) ) cycle
                    if (hardorder == 0) then
                        msq(j,k) = msqall(j,k,0)
                    elseif (hardorder == 1) then
                        msq(j,k) = sum(msqall(j,k,[0,1]))
                    elseif (hardorder == 2) then
                        msq(j,k) = sum(msqall(j,k,[0,1,2]))
                    elseif (hardorder == 3) then
                        msq(j,k) = sum(msqall(j,k,[0,1,2,3]))
                    endif
                enddo; enddo
                call resummation_BNRptlumxmsq(p,xx,r(ndim-3),r(ndim-2),qt,sqrt(qsq),Rcut,msq,hardorder,xmsq)
            else
! qt resummation
            xmsq = 0._dp
                do j=-nflav,nflav; do k =-nflav,nflav
                if ( all(msqall(j,k,0:hardorder) == 0._dp) ) cycle
                res = 0._dp

                ! hard function set to 0._dp, not required for resummed piece
                call resummation(qsq, qt, xx(1), xx(2), j, k, 0._dp, sqrt(musq), &
                    order, muMult=muMultRes, res=res)

                if (hardorder == 0) then
                   xmsq = xmsq + res * msqall(j,k,0)
                elseif (hardorder == 1) then
                    xmsq = xmsq + res * sum(msqall(j,k,[0,1]))
                elseif (hardorder == 2) then
                    xmsq = xmsq + res * sum(msqall(j,k,[0,1,2]))
                elseif (hardorder == 3) then
                    xmsq = xmsq + res * sum(msqall(j,k,[0,1,2,3]))
                    !write (*,*) j,k, msqall(j,k,3)/sum(msqall(j,k,[0,1,2]))
                endif
            enddo; enddo
            endif

            resout = xmsq
            return
        endif

        if (present(resexpout)) then
! jet-veto resummation
            if (usept) then
                msq0(:,:)=msqall(:,:,0)
                msq1(:,:)=msqall(:,:,1)/(as/4._dp/pi)
!                msq2(:,:)=msqall(:,:,2)/(as/4._dp/pi)**2
! msq2 is not required since it is also removed in the above-cut piece,
! c.f. calls to hard_routine in qtlumxmsq.f, GLYlumxmsq.f, ptlumxmsq.f, BNRptlumxmsq.f
                msq2(:,:)=0._dp
                call resexp_BNRptlumxmsq(p,xx,r(ndim-3),r(ndim-2),qt,sqrt(qsq),Rcut,hardorder,msq0,msq1,msq2,xmsq)
                xmsq = -xmsq  ! flip sign in order to form difference for matching corrections
            else
            xmsq = 0._dp
            do j=-5,5; do k =-5,5
                ! cycle if no contribution through NLO
                if ( all(msqall(j,k,[0,1]) == 0._dp) ) cycle

                resexp = 0._dp

                ! coefficient in as/4/pi
                hard(1) = msqall(j,k,1) / msqall(j,k,0) / (as/4._dp/pi)

                if (order >= 7) then
                    hard(2) = msqall(j,k,2) / msqall(j,k,0) / (as/4._dp/pi)**2
                endif

                call resummation(qsq, qt, xx(1), xx(2), j, k, hard(1), sqrt(musq), order, resexp=resexp, hard2=hard(2))

                ! matching correction
                xmsq = xmsq - msqall(j,k,0)*resexp

                !write (*,*) "prefactor", j,k, msqall(j,k,0)*fbgev2
            enddo; enddo
            endif

            resexpout = xmsq
            return
        endif


    end subroutine computeResint

    function resexpint(r,wgt)
        use ieee_arithmetic
        use types
        use Scalevar
        use PDFerrors
        use SCET
        use MCFMStorage
        use qtResummation
        use LHAPDF
        use MCFMSetupPlots
        use scaleset_m
        use MCFMSettings, only: resexp_linPC
        use SafetyCuts, only : passed_smallnew
        use ptveto
        implicit none
        include 'constants.f'
        include 'mxpart.f'
        include 'vegas_common.f'
        include 'kprocess.f'
        include 'kpart.f'
        include 'energy.f'
        include 'npart.f'
        include 'dynamicscale.f'
        include 'initialscales.f'
        include 'scalevar.f'
        include 'scale.f'
        include 'facscale.f'
        include 'qcdcouple.f'
        include 'couple.f'
        include 'nlooprun.f'
        include 'x1x2.f'
        include 'xmin.f'
        include 'nf.f'
        include 'nflav.f'
        include 'nwz.f'
        include 'nproc.f'
        include 'ipsgen.f'
        include 'masses.f'
        include 'beamtype.f'
        logical:: bin
        common/bin/bin
        real(dp) :: BrnRat
        common/BrnRat/BrnRat

        real(dp) :: resexpint
        real(dp), intent(in) :: r(mxdim), wgt

        real(dp):: p(mxpart,4),pjet(mxpart,4),xmsq
        real(dp) :: val,pswt,flux

        real(dp) :: qsq
        real(dp) :: p_born(mxpart,4)

        ! functions
        logical :: includedipole
        real(dp) :: ptpure, twomass, threemass, pttwo
        real(dp) :: xjac
        real(dp) :: qt, phi
        real(dp) :: hard(2)
        real(dp) :: res, resexp, xmsqMu, xmsqPdf
        real(dp) :: msq(-nf:nf, -nf:nf)
        real(dp) :: phaseqt

        real(dp) :: pt34com
        common/pt34com/pt34com
!$omp threadprivate(/pt34com/)

        integer :: j,k

        logical :: withRecoil, withoutRecoil
        real(dp) :: withRecoiln, withoutRecoiln

        resexpint = 0._dp

        currentPDF = 0
        currentNd = 0

        if (doScalevar .and. bin) then
            scalereweight(:) = 1._dp
        endif

        p(:,:) = 0._dp
        pjet(:,:) = 0._dp

        if (qtmaxResexp .eq. qtminResexp) then
            qt = qtminResexp
            xjac = 2._dp*qt
        else
            qt = r(ndim)*(qtmaxResexp-qtminResexp) + qtminResexp
            xjac = 2._dp*qt*(qtmaxResexp-qtminResexp)
        endif
        phi = 2._dp*pi*r(ndim-1)

!       ! for debugging purposes we can fix qt here
        !qt = 1.5_dp
        !xjac = 2._dp*qt

        pt34com = qt

        if (qt > qtmaxResexp) then
            resexpint = 0._dp
            return
        endif

        if (qt < qtminResexp) then
            resexpint = 0._dp
            return
        endif

        if (usept) then
          qt = jetptveto
          xjac=1._dp
        endif

        p(:,:) = 0._dp

        if (phaseUseQt) then
            phaseqt = qt
        else
            phaseqt = 0._dp
        endif

        if (nproc == 1 .or. nproc == 6) then
            npart = 2
            if (gen_res(r,phaseqt,p,pswt) .eqv. .false. ) then
                resexpint = 0._dp
                return
            endif
        elseif (any(nproc == [31,32])) then ! Z -> e e or Z -> 3*(nu nubar)
            npart = 2
            ! generate phase space at zero pT
            if (gen_res(r,phaseqt,p,pswt) .eqv. .false. ) then
                resexpint = 0._dp
                return
            endif
        elseif (kcase == kggfus0) then ! H -> ga ga
            npart = 2
            if (gen_res(r,phaseqt,p,pswt) .eqv. .false.) then
                resexpint = 0._dp
                return
            endif
        elseif (kcase == kHi_Zga) then ! H -> Z ga
            npart=3
            call gen3h(r,p,pswt,*888)
        elseif (kcase == kWHbbar .or. kcase == kWHgaga) then
            npart=4
            call genVH(r,p,pswt,*888)
        elseif (kcase == kZHbbar .or. kcase == kZHgaga) then
            npart=4
            call genVH(r,p,pswt,*888)
        elseif (kcase == kZH__WW) then
            npart=6
            call gen6(r,p,pswt,*888)
        elseif (any(nproc == [290,295])) then ! W -> e nu
            npart = 3
            if (ipsgen == 1) then
                call gen_Vphotons_jets(r,1,0,p,pswt,*888)
            elseif (ipsgen == 2) then
                call gen_Vphotons_jets_dkrad(r,1,0,p,pswt,*888)
            else
              error stop
            endif
            wgaips: block
                real(dp) :: s34, s345, wt34, wt345, wtips(2)
                real(dp) :: dot
              s34 = 2._dp*dot(p,3,4)
              s345 = s34+2._dp*dot(p,3,5)+2._dp*dot(p,4,5)
              wt34 = (s34-wmass**2)**2+(wmass*wwidth)**2
              wt345 = (s345-wmass**2)**2+(wmass*wwidth)**2
              wtips(1) = wt345
              wtips(2) = wt34
              pswt = pswt*wtips(ipsgen)/(wtips(1)+wtips(2))
            end block wgaips
        elseif (any(nproc == [300,305])) then ! Z -> e+ e- + photon or Z -> 3*(nu nubar) + photon
            npart = 3
            if (ipsgen == 1) then
                call gen_Vphotons_jets(r,1,0,p,pswt,*888)
            elseif (ipsgen == 2) then
                call genVphoton(r,p,pswt,*888)
            else
              error stop
            endif
            zgaips: block
                real(dp) :: s34, s345, wt34, wt345, wtips(2)
                real(dp) :: dot
              s34 = 2._dp*dot(p,3,4)
              s345 = s34+2._dp*dot(p,3,5)+2._dp*dot(p,4,5)
              wt34 = s34*(s34-zmass**2)**2+(zmass*zwidth)**2
              wt345 = s345*(s345-zmass**2)**2+(zmass*zwidth)**2
              wtips(1) = wt345
              wtips(2) = wt34
              pswt = pswt*wtips(ipsgen)/(wtips(1)+wtips(2))
            end block zgaips
        elseif (nproc == 285 &      ! q qbar -> gam gam
           .or. nproc == 2851 &     ! g g -> gam gam (massless quarks)
           .or. nproc == 2852) then ! g g -> gam gam (including top loops)
            npart = 2
            call gen_photons_jets_res(r,2,0,p,pswt,*888)
        elseif ( ((kcase == kWWqqbr) .or. (kcase == kZZlept)) &   ! WW and ZZ gg contribution (ipsgen = 2)
            .and. (ipsgen == 2)) then
            npart=4
            call gen4handc(r,p,pswt,*888)
        elseif ((kcase == kWWqqbr) .or. (kcase == kWZbbar) .or. (kcase == kZZlept)) then  ! WW, ZZ, WZ
            npart = 4
            call gen4(r,p,pswt,*888)
        else
            write (*,*) __FILE__//": undefined nproc, line ", __LINE__
            error stop 
        endif


        if (.not. passed_smallnew(p,npart,1e-9_dp)) then
            goto 888
        endif

! back-to-back check
        if ((kcase == kWWqqbr) .or. (kcase == kWZbbar) .or. (kcase == kZZlept)) then
            if (pttwo(3,4,p) < 1.e-3_dp) goto 888
        endif



        !call smallnew(p,npart,*999)
        if (.false.) then
 888        continue
            resexpint = 0._dp
            return
        endif


        if (any(ieee_is_nan(p(1:npart+2,:)) .eqv. .true.)) then
          !write(6,*) 'Discarding NaN or infinite phase space point'
          resexpint = 0._dp
          return
        endif

! these are calculated before the boost
        xx(1)=-2._dp*p(1,4)/sqrts
        xx(2)=-2._dp*p(2,4)/sqrts

        if ( (xx(1) > one)  .or. (xx(2) > one) &
         .or.(xx(1) < xmin) .or. (xx(2) < xmin)) then
            resexpint = 0._dp
            return
        endif

! save before boost for matrix element evaluation
        p_born(1:npart+2,:) = p(1:npart+2,:)
! mostly for plotting purposes
        call recoilBoost(qt, phi, p)
        if (resexp_linPC) then
            withoutRecoil = includedipole(0,p_born)
            withRecoil = includedipole(0,p)

            if (boostBorn) then
                p_born = p
            endif

            if (withRecoil .eqv. withoutRecoil) then
                resexpint = 0._dp
                return
            endif

            if (withRecoil .eqv. .true.) then
                withRecoiln = 1._dp
            else
                withRecoiln = 0._dp
            endif

            if (withoutRecoil .eqv. .true.) then
                withoutRecoiln = 1._dp
            else
                withoutRecoiln = 0._dp
            endif

            pswt = pswt*(withRecoiln - withoutRecoiln)

        else ! normal resexp operation

        if (boostBorn) then
            p_born = p
        endif

! check cuts, fills ptildejet
        includeTaucutgrid(0) = .TRUE. 
        if (includedipole(0,p) .eqv. .FALSE. ) then
            resexpint = 0._dp
            return
          endif
        endif

        call storeptilde(0,p)

! fill pjet array
        call getptildejet(0,pjet)

!! In this part of the calculation just use the normal scale choice!
!        if (usept) then
!! note that, for this suite of routines
!!  mu = facscale, central value = ptveto (qt)
!! muh = scale,  central values = sqrt(qsq)
!          qsq = twomass(1,2,p)**2
!!          call usescales(sqrt(qsq)*initscale, qt*initfacscale)
!          call usescales(qt*initscale, qt*initfacscale)
!        else
        if (dynamicscale) then
            if (boostBorn) then
                use_resummation_recoil = .true.
                resummation_recoil = qt
            else
                use_resummation_recoil = .false.
                resummation_recoil = 0._dp
            endif
            call scaleset(initscale, initfacscale, p)
        else
            call usescales(initscale, initfacscale)
        endif
!        endif

        flux = fbGeV2/(two*xx(1)*xx(2)*sqrts**2)

        ! central value
        xmsq = 0._dp
        call computeResint(r, p, p_born, qt, kresorder, xx, resexpout=xmsq)

        if (doScalevar .and. bin) then
            do j=1,maxscalevar
                call computeResint(r,p,p_born,qt,kresorder,xx, resexpout=xmsqMu, &
                    muMultHard_in=scalevarmult(j), muMultRes_in=facscalevarmult(j))
                scalereweight(j) = xmsqMu / xmsq
            enddo
        endif

        if ((maxPDFsets > 0) .and. bin) then
            pdfreweight(:) = 0._dp
            do j=1,maxPDFsets
                currentPDF = j
                call computeResint(r,p, p_born, qt, kresorder, xx, resexpout=xmsqPdf)
                pdfreweight(currentPDF) = flux*xjac*pswt*(xmsq-xmsqPdf)/BrnRat*wgt
            enddo
        endif

        resexpint = flux*xjac*pswt*xmsq/BrnRat

        val=resexpint*wgt

        if (ieee_is_nan(val) .OR. ( .NOT. ieee_is_finite(val))) then
            !write(6,*) 'Discarded NaN, val=',val
            resexpint = 0._dp
            return
        endif
              
        if (bin) then
            includeTaucutgrid(0) = .true. 
            call nplotter_new(pjet,val)
        endif

    end function resexpint

end module
