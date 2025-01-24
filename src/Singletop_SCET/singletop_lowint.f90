!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module singletop_lowint
      use singletop_int
      implicit none

      private

      public :: lowint

      contains

      function lowint(r,wgt)
          use Superhisto, only : shtmpreset, shtmpcommit
          use MCFMSetupPlots, only: nplotter_new
          implicit none
          include 'nf.f'
          include 'mxdim.f'
          include 'mxpart.f'
          include 'npart.f'
          include 'x1x2.f'
          include 'kprocess.f'
          include 'nproc.f'
          include 'constants.f'
          include 'energy.f'
          include 'debug.f'
          include 'xmin.f'
          include 'maxwt.f'
          include 'taucut.f'! usescet
          include 'beamtype.f'

          logical :: bin
          common/bin/bin

          real(dp) :: lowint
          real(dp), intent(in) :: r(mxdim), wgt

          ! external functions
          logical :: includedipole

          real(dp) :: p(mxpart,4), pjet(mxpart,4)
          real(dp) :: pswt

          real(dp) :: lowint_b(max_bcontrib, max_corr_on_beam), &
              tmp_b(max_bcontrib, max_corr_on_beam)
          real(dp), allocatable :: scalereweight_b(:,:,:)
          real(dp), allocatable :: pdfreweight_b(:,:,:)
          real(dp), allocatable :: scetreweight_b(:,:)

          logical :: passed_taucut(max_corr_on_beam)

          real(dp) :: msqall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)
          real(dp) :: msqall_tmp(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

          logical :: bcontribs(max_bcontrib)

          integer :: truescalevar

          integer :: j,m,ipdf

          lowint = 0._dp
          lowint_b(:,:) = 0._dp
          p(:,:) = 0._dp
          pjet(:,:) = 0._dp
          currentPDF = 0

          if (.not. gen_singletop(r,p,pswt)) then
              lowint = 0._dp
              return
          endif

          if (any(ieee_is_nan(p(1:npart,:)))) then
              if (debug) then
                   write(6,*) 'Discarding NaN or infinite phase space point'
              endif
              lowint = 0._dp
              return
          endif

          !!!
          !!! note that we do not have a smallnew/smalltau here (yet)
          !!!

          xx(1)=-2._dp*p(1,4)/sqrts
          xx(2)=-2._dp*p(2,4)/sqrts
          if ((xx(1) >  1._dp) .or. (xx(2) >  1._dp) &
              .or.(xx(1) < xmin)   .or. (xx(2) < xmin)) then
              lowint = 0._dp
              return
          endif

          ! bcontribs are
          ! 1: main b
          ! 2: additional pair b~ b~
          ! 3: additional b~
          ! 4: additional pair b b~

          ! first disable all, next we determine which ones are needed
          ! dependent on the process and on cuts in includedipole
          bcontribs(:) = .false.

          ! set up all possible b contributions
          if (kcase == kbq_tpq) then
              bcontribs(1) = .true.
          elseif (kcase == kbq_tpq_jet) then
              bcontribs(1) = .true.

              if (currentContrib == 2) then
                  bcontribs(3) = .true.
              endif
          endif

          ! first determine b contributions that do not vanish
          ! and the taucuts on the two different beam contributions
          passed_taucut = .true.
          if (doMultitaucut) then
              allocate(scetreweight_b(size(scetreweight),max_corr_on_beam))

              if (.not. usescet) then ! this is not a scet piece, so fill with 1
                  scetreweight_b(:,:) = 1d0 
              endif
          endif

          do b_contrib=1,max_bcontrib
              if (bcontribs(b_contrib) .eqv. .false.) cycle

              if (usescet) then

                  do m=1,maxbeams
                      corr_on_beam = beams_enabled(m)
                      if (currentContrib == 1) then
                          passed_taucut(corr_on_beam) = &
                              passed_taucut_light(p,scetreweight_b(:,corr_on_beam))
                      elseif (currentContrib == 2) then
                          passed_taucut(corr_on_beam) = &
                              passed_taucut_heavyprod(p,scetreweight_b(:,corr_on_beam))
                      elseif (currentContrib == 3) then
                          passed_taucut(corr_on_beam) = &
                              passed_taucut_decay(p,scetreweight_b(:,corr_on_beam))
                      endif
                  enddo
              endif

              currentNd = 0
              bcontribs(b_contrib) = includedipole(0,p)
          enddo

          if (all(bcontribs .eqv. .false.)) then
              lowint = 0._dp
              return
          endif

!!! central value

          ! set all DDIS singletop scale variables based on p
          ! or as constant or dynamic if use_DDIS is .false.
          call singletop2_scale_setup(p)
          ! fill pdf ingredients
          call calc_singletop_pdfs_lord(xx)
          ! fill msq ingredients
          call lowint_fillmsq(p, bcontribs, msqall)
          lowint_b = lowint_assemble(msqall,bcontribs,pswt,xx)

           if (any(ieee_is_nan(lowint_b)) .or. (.not. all(ieee_is_finite(lowint_b)))) then
               !write(6,*) 'Discarding NaN matrix element!'
               lowint = 0._dp
               return
           endif

!!! pdf uncertainties
! reuses msqall, have it before scale variation
          if (maxPDFsets > 0 .and. bin) then
              allocate(pdfreweight_b(maxPDFsets,max_bcontrib,max_corr_on_beam), source=0._dp)

              do ipdf=1,maxPDFsets
                  currentPDF = ipdf
                  call calc_singletop_pdfs_lord(xx)
                  if (doPDFAlphas) then
                      ! update values of alpha_s
                      call singletop2_scale_setup(p)
                      ! recompute with new alphas
                      call lowint_fillmsq(p, bcontribs, msqall)
                  endif
                  tmp_b = lowint_assemble(msqall,bcontribs,pswt,xx)
                  do b_contrib=1,max_bcontrib
                      if (bcontribs(b_contrib) .eqv. .false.) cycle
                      do m=1,maxbeams
                          corr_on_beam = beams_enabled(m)
                          pdfreweight_b(currentPDF, b_contrib, corr_on_beam) = &
                              (lowint_b(b_contrib,corr_on_beam) - tmp_b(b_contrib,corr_on_beam))*wgt
                      enddo
                  enddo

              enddo

              currentPDF = 0
          endif

!!! scale variation

          if (doScalevar .and. bin) then
              allocate(scalereweight_b(maxscalevar, max_bcontrib, max_corr_on_beam), source=0._dp)

              ! let's not be clever and rather be general, for now.
              ! * kbq_tpq  does not depend on alphas and we can just vary factscale
              ! * kbq_tpq_jet does depend on alphas, and for each corr_on_beam
              !   we could rescale with the appropriate alphas.

              ! but in the following way it's completely general

              if (use_DDIS) then
                  truescalevar = (maxscalevar-1)/2
              else
                  truescalevar = maxscalevar
              endif

              do j=1,maxscalevar
                  if (j <= truescalevar) then
                      call singletop2_scale_setup(p, mult_in_ren=scalevarmult(j), mult_in_fac=facscalevarmult(j))
                  elseif (j == (truescalevar + 1)) then
                      call singletop2_scale_setup(p, mult_in_ren=1d0, mult_in_fac=1d0, forcemt=.true.)
                  else
                      call singletop2_scale_setup(p, mult_in_ren=scalevarmult(j-truescalevar-1), &
                          mult_in_fac=facscalevarmult(j-truescalevar-1), forcemt=.true.)
                  endif
                  call calc_singletop_pdfs_lord(xx)
                  call lowint_fillmsq(p, bcontribs, msqall_tmp)
                  tmp_b = lowint_assemble(msqall_tmp,bcontribs,pswt,xx)

                  do b_contrib=1,max_bcontrib
                      if (bcontribs(b_contrib) .eqv. .false.) cycle
                      do m=1,maxbeams
                          corr_on_beam = beams_enabled(m)
                          if (lowint_b(b_contrib,corr_on_beam) == 0._dp) cycle
                          scalereweight_b(j, b_contrib, corr_on_beam) = &
                              tmp_b(b_contrib,corr_on_beam) / lowint_b(b_contrib,corr_on_beam)
                      enddo
                  enddo
              enddo

              ! reset to unmodified scales
              call singletop2_scale_setup(p)

              ! pdfs can stay modified, they are at most used by pdf variation
              ! and will be recalculated
          endif


          if (bin) then
              call getptildejet(0,pjet)
              currentNd = 0

              do m=1,maxbeams
                  corr_on_beam = beams_enabled(m)
                  includeTaucutgrid(currentNd) = passed_taucut(corr_on_beam)

                  do b_contrib=1,max_bcontrib
                      if (bcontribs(b_contrib) .eqv. .false.) cycle

                      if (doMultitaucut) then
                          scetreweight(:) = scetreweight_b(:,corr_on_beam)
                      endif

                      if (doScalevar) then
                          scalereweight(1:maxscalevar) = scalereweight_b(:,b_contrib,corr_on_beam)
                      endif
                      if (maxPDFsets > 0) then
                          pdfreweight(:) = pdfreweight_b(:,b_contrib,corr_on_beam)
                      endif

                      call nplotter_new(pjet,lowint_b(b_contrib,corr_on_beam)*wgt)
                  enddo
              enddo
              call threadStorageOp(shtmpcommit)
          endif

          lowint = 0._dp
          do m=1,maxbeams
              corr_on_beam = beams_enabled(m)

              if (passed_taucut(corr_on_beam)) then
                  lowint = lowint + sum(lowint_b(:,corr_on_beam))
              endif
          enddo

          if (abs(lowint*wgt) > wtmax) then
              wtmax = abs(lowint*wgt)
          endif

      end function

end module
