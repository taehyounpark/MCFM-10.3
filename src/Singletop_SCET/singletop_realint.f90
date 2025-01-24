!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module singletop_realint
      use singletop_int
      implicit none

      private

      public :: realint

      contains

      function realint(r,wgt)
          use Superhisto, only : shtmpreset, shtmpcommit
          use SafetyCuts, only : passed_smallnew
          use MCFMSetupPlots, only: nplotter_new
          use singletop_rvint, only: rvint
          implicit none
          include 'nf.f'
          include 'mxdim.f'
          include 'mxpart.f'
          include 'npart.f'
          include 'x1x2.f'
          include 'kprocess.f'
          include 'nproc.f'
          include 'ipsgen.f'
          include 'constants.f'
          include 'energy.f'
          include 'debug.f'
          include 'xmin.f'
          include 'maxwt.f'
          include 'taucut.f'! usescet
          include 'maxd.f'
          include 'incldip.f'
          include 'beamtype.f'

          logical :: bin
          common/bin/bin

          real(dp) :: realint 
          real(dp), intent(in) :: r(mxdim), wgt

          ! external functions
          logical :: includedipole

          real(dp) :: p(mxpart,4), pjet(mxpart,4)
          real(dp) :: ptilde(mxpart,4)
          real(dp) :: pswt

          real(dp), allocatable :: realint_b(:,:,:), tmp_b(:,:,:)
          real(dp), allocatable :: scalereweight_b(:,:,:,:)
          real(dp), allocatable :: pdfreweight_b(:,:,:,:)
          real(dp), allocatable :: scetreweight_b(:,:,:)


          real(dp) :: msqall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)
          real(dp) :: msqall_tmp(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

          integer :: ndmax
          logical, allocatable :: passed_taucut(:,:)
          logical, allocatable :: bcontribs(:,:)

          real(dp), allocatable :: msqcall(:,:,:,:,:)
          real(dp), allocatable :: msqcall_tmp(:,:,:,:,:)

          real(dp), allocatable :: singletop2_dipscale_save(:,:)

          integer :: truescalevar

          integer :: j,ipdf,nd,m

          if ((currentContrib == 4 .and. ipsgen == 5) .or. &
              (currentContrib == 5 .and. ipsgen == 7) .or. &
              (currentContrib == 6 .and. ipsgen == 9)) then
              realint = rvint(r,wgt)
              return
          endif

           p(:,:) = 0._dp
           pjet(:,:) = 0._dp
           currentPDF = 0

           ! SCET for all interference pieces
           if (currentContrib > 3) then
               usescet = .true.
           endif

           if (.not. gen_singletop(r,p,pswt)) then
               realint= 0._dp
               return
           endif
 
           if (any(ieee_is_nan(p(1:npart+2,:)))) then
               write(6,*) 'Discarding NaN phase space point'
               realint = 0._dp
               return
           endif
 
           if (.not. passed_smallnew(p,npart,1d-8)) then
               realint = 0._dp
               return
           endif


           ! Unfortunately the small taucut alone is not sufficient
           ! so we added smallnew above. 
!          if (currentContrib == 1) then
!              do m=1,maxbeams
!                  corr_on_beam = beams_enabled(m)
!                  if (.not. passed_taucut_light(p, taucut_in=1d-6)) then
!                      realint = 0._dp
!                      return
!                  endif
!              enddo
!          elseif (currentContrib == 2) then
!              do m=1,maxbeams
!                  corr_on_beam = beams_enabled(m)
!                  if (.not. passed_taucut_heavyprod(p, taucut_in=1d-6)) then
!                      realint = 0._dp
!                      return
!                  endif
!              enddo
!          elseif (currentContrib == 3) then
!              !error stop "implement safety taucut in singletop_realint for currentContrib=3"
!          endif
 
           xx(1)=-2._dp*p(1,4)/sqrts
           xx(2)=-2._dp*p(2,4)/sqrts
           if ((xx(1) >  1._dp) .or. (xx(2) >  1._dp) &
               .or.(xx(1) < xmin)   .or. (xx(2) < xmin)) then
               realint = 0._dp
               return
           endif
 
           ! bcontribs are
           ! 1: main b
           ! 2: additional b
           ! 3: additional b~
           ! 4: additional pair b b~
           ! 5: additional pair b~ b~
 
           ! first disable all, next we determine which ones are needed
           ! dependent on the process and on cuts in includedipole
           ! and for all dipole contributions!
           ndmax = 0
           if (kcase == kbq_tpq) then
               if (currentContrib == 1) then
                   ndmax = 8
               elseif (currentContrib == 2) then
                   ndmax = 6
               elseif (currentContrib == 3) then
                   ndmax = 2
               elseif (any(currentContrib == [4,5])) then
                   ndmax = 8
               elseif (currentContrib == 6) then
                   ndmax = 6
               endif
           elseif (kcase == kbq_tpq_jet) then
               if (currentContrib == 1) then
                   ndmax = 32
               elseif (currentContrib == 2) then
                   ndmax = 20
               elseif (currentContrib == 3) then
                   ndmax = 20
               endif
           endif

           allocate(passed_taucut(0:ndmax,max_corr_on_beam), source=.true.)
           allocate(bcontribs(0:ndmax,max_bcontrib), source=.false.)
 
           ! set up all possible b contributions
           if (kcase == kbq_tpq) then
               bcontribs(:,1) = .true.

               if (currentContrib == 2) then
                   ! initial state gluon splits into b b~, b~ goes extra
                   bcontribs(:,3) = .true.
               elseif (currentContrib == 4 .and. ipsgen == 4) then
                   bcontribs(:,3) = .true.
               elseif (currentContrib == 6 .and. ipsgen == 8) then
                   bcontribs(:,3) = .true.
               endif
           elseif (kcase == kbq_tpq_jet) then
               if (currentContrib == 1) then
                   bcontribs(:,[1,2,3,4]) = .true.
               elseif (currentContrib == 2) then
                   ! just an additional b (case 2) can't appear
                   bcontribs(:,[1,3,4,5]) = .true.
               elseif (currentContrib == 3) then
                   bcontribs(:,[1,4]) = .true.
               endif
           endif

           call singletop2_scale_setup(p)

           ! this array is huge, so filling it is expensive
           ! so we don't initialize with zero here, but in the routines
           allocate(msqcall(ndmax,-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam))

           ! fill msq ingredients 
           ! this implicitly sets all dipole scales which we need to
           ! evaluate maketaucut and includedipole 
           call realint_fillmsq(p, ndmax, bcontribs, msqall, msqcall)
           ! reset scales, which have been modified by ptrans
           call singletop2_scale_setup(p)

           ! now we can determine contributions that vanish
           passed_taucut = .true.
           call storeptilde(0,p)

           if (doMultitaucut) then
               allocate(scetreweight_b(size(scetreweight),0:ndmax,max_corr_on_beam))

               if (.not. usescet) then ! this is not a scet piece, so fill with 1
                   scetreweight_b(:,:,:) = 1d0 
               endif
           endif

           do nd=0,ndmax
              ! otherwise passed_taucut will be called with zero ptilde
              if ((nd > 0) .and. (incldip(nd) .eqv. .false.)) cycle

              call getptilde(nd,ptilde)

              do b_contrib=1,max_bcontrib
                  if (bcontribs(nd,b_contrib) .eqv. .false.) cycle

                  currentNd = nd ! to correctly fill jetcontent
                  bcontribs(nd,b_contrib) = includedipole(nd,ptilde)

                  if (bcontribs(nd,b_contrib) .eqv. .false.) cycle

                  if (usescet) then
                      do m=1,maxbeams
                          corr_on_beam = beams_enabled(m)
                           if (currentContrib == 1) then
                               passed_taucut(nd,corr_on_beam) = &
                                   passed_taucut_light(ptilde,scetreweight_b(:,nd,corr_on_beam))
                           elseif (currentContrib == 2) then
                               passed_taucut(nd,corr_on_beam) = &
                                   passed_taucut_heavyprod(ptilde,scetreweight_b(:,nd,corr_on_beam))
                           elseif (currentContrib == 3) then
                               passed_taucut(nd,corr_on_beam) = &
                                   passed_taucut_decay(ptilde,scetreweight_b(:,nd,corr_on_beam))
                           elseif (currentContrib == 4 .and. ipsgen == 4) then
                               passed_taucut(nd,corr_on_beam) = &
                                   passed_taucut_lxh(ptilde,scetreweight_b(:,nd,corr_on_beam))
                           elseif (currentContrib == 5 .or. currentContrib == 6) then
                               passed_taucut(nd,corr_on_beam) = &
                                   passed_taucut_decay(ptilde,scetreweight_b(:,nd,corr_on_beam), &
                                        ptilde(5,:) + ptilde(7,:))
                           endif
                      enddo
                  endif

              enddo
           enddo

           if (all(bcontribs .eqv. .false.)) then
               realint = 0._dp
               return
           endif

           ! only after dipole scales are set
           call calc_singletop_pdfs_real(xx,ndmax)

           allocate(realint_b(0:ndmax,max_bcontrib,max_corr_on_beam), source=0._dp)
           realint_b = realint_assemble(ndmax,msqall,msqcall,bcontribs,pswt,xx)

           if (any(ieee_is_nan(realint_b)) .or. (.not. all(ieee_is_finite(realint_b)))) then
               !write(6,*) 'Discarding NaN matrix element!'
               realint = 0._dp
               return
           endif


           if (doScalevar .and. bin) then
               allocate(scalereweight_b(0:ndmax,maxscalevar,max_bcontrib,max_corr_on_beam), source=0._dp)
               if (.not. allocated(msqcall_tmp)) then
                   allocate(msqcall_tmp(ndmax,-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam), source=0._dp)
               endif
               if (.not. allocated(tmp_b)) then
                   allocate(tmp_b(0:ndmax,max_bcontrib,max_corr_on_beam), source=0._dp)
               endif

               ! we might need this later for the pdf uncertainties
               ! but let's restore the environment as before in any case
               allocate(singletop2_dipscale_save(ndmax,2),source=0._dp)
               singletop2_dipscale_save(1:ndmax,:) = singletop2_dipscale(1:ndmax,:)

               ! in principle we can reuse msqall and msqcall
               ! since these are no longer accessed later on. let's do so.

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
                   call realint_fillmsq(p, ndmax, bcontribs, msqall_tmp, msqcall_tmp)

                   ! reset scales, which have been modified by ptrans
                   if (j <= truescalevar) then
                       call singletop2_scale_setup(p, mult_in_ren=scalevarmult(j), mult_in_fac=facscalevarmult(j))
                   elseif (j == (truescalevar + 1)) then
                       call singletop2_scale_setup(p, mult_in_ren=1d0, mult_in_fac=1d0, forcemt=.true.)
                   else
                       call singletop2_scale_setup(p, mult_in_ren=scalevarmult(j-truescalevar-1), &
                           mult_in_fac=facscalevarmult(j-truescalevar-1), forcemt=.true.)
                   endif

                   ! pdfs after dipole scales are set
                   call calc_singletop_pdfs_real(xx,ndmax)
                   tmp_b = realint_assemble(ndmax,msqall_tmp,msqcall_tmp,bcontribs,pswt,xx)

                   do b_contrib=1,max_bcontrib
                       if (all(bcontribs(:,b_contrib) .eqv. .false.)) cycle
                       do m=1,maxbeams
                           corr_on_beam = beams_enabled(m)
                           do nd=0,ndmax
                               if (bcontribs(nd,b_contrib) .eqv. .false.) cycle
                               if (realint_b(nd,b_contrib,corr_on_beam) == 0._dp) cycle

                               scalereweight_b(nd, j, b_contrib, corr_on_beam) = &
                                   tmp_b(nd,b_contrib,corr_on_beam) / &
                                      realint_b(nd,b_contrib,corr_on_beam)
                           enddo
                       enddo
                   enddo
               enddo

               ! reset old scales
               call singletop2_scale_setup(p)
               ! restore old dipscales
               singletop2_dipscale(1:ndmax,:) = singletop2_dipscale_save(1:ndmax,:)
           endif

           if (maxPDFsets > 0 .and. bin) then
               allocate(pdfreweight_b(0:ndmax,maxPDFsets,max_bcontrib,max_corr_on_beam), source=0._dp)

               if (.not. allocated(msqcall_tmp)) then
                   allocate(msqcall_tmp(ndmax,-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam), source=0._dp)
               endif
               if(.not. allocated(tmp_b)) then
                   allocate(tmp_b(0:ndmax,max_bcontrib,max_corr_on_beam), source=0._dp)
               endif

               do ipdf=1,maxPDFsets
                   currentPDF = ipdf
                   if (doPDFalphas) then
                       call singletop2_scale_setup(p)
                       call realint_fillmsq(p, ndmax, bcontribs, msqall_tmp, msqcall_tmp)
                       ! reset scales, which have been modified by ptrans
                       call singletop2_scale_setup(p)
                       ! pdfs after dipole scales are set
                       call calc_singletop_pdfs_real(xx,ndmax)
                       tmp_b = realint_assemble(ndmax,msqall_tmp,msqcall_tmp,bcontribs,pswt,xx)
                   else
                       call singletop2_scale_setup(p)
                       call calc_singletop_pdfs_real(xx,ndmax)
                       tmp_b = realint_assemble(ndmax,msqall,msqcall,bcontribs,pswt,xx)
                   endif

                   do b_contrib=1,max_bcontrib
                       if (all(bcontribs(:,b_contrib) .eqv. .false.)) cycle
                       do m=1,maxbeams
                           corr_on_beam = beams_enabled(m)
                           do nd=0,ndmax
                               if (bcontribs(nd,b_contrib) .eqv. .false.) cycle

                               pdfreweight_b(nd,ipdf,b_contrib,corr_on_beam) = &
                                   (realint_b(nd,b_contrib,corr_on_beam) - tmp_b(nd,b_contrib,corr_on_beam))*wgt
                           enddo
                       enddo
                   enddo
               enddo

               ! reset environment (otherwise bookplot will fail)
               currentPDF = 0
           endif

!!! BINNING
           if (bin) then
               do m=1,maxbeams
                   corr_on_beam = beams_enabled(m)
                   do b_contrib=1,max_bcontrib
                       if (all(bcontribs(:,b_contrib) .eqv. .false.)) cycle

                       do nd=0,ndmax
                           if (bcontribs(nd,b_contrib) .eqv. .false.) cycle

                           ! Each nd corresponds to a unique
                           ! corr_on_beam,b_contrib combination, so we should be
                           ! save to have this here.

                           ! incldip gets overwritten after scale variation and
                           ! pdf uncertainties, but *should* stay the same.
                           if ((nd > 0) .and. (incldip(nd) .eqv. .false.)) cycle

                           currentNd = nd
                           call getptildejet(nd,pjet)
                           includeTaucutgrid(currentNd) = passed_taucut(nd,corr_on_beam)

                           if (doMultitaucut) then
                               scetreweight(:) = scetreweight_b(:,nd,corr_on_beam)
                           endif

                           if (doScalevar) then
                               scalereweight(1:maxscalevar) = scalereweight_b(nd,1:maxscalevar,b_contrib,corr_on_beam)
                           endif

                           if (maxPDFsets > 0) then
                               pdfreweight(:) = pdfreweight_b(nd,:,b_contrib,corr_on_beam)
                           endif

                           call nplotter_new(pjet,realint_b(nd,b_contrib,corr_on_beam)*wgt)

                       enddo

                   enddo
               enddo
               call threadStorageOp(shtmpcommit)
           endif

!! finalization for vegas

           realint = 0._dp
           do nd=0,ndmax
               do m=1,maxbeams
                   corr_on_beam = beams_enabled(m)
                   if (passed_taucut(nd,corr_on_beam)) then
                       do b_contrib=1,max_bcontrib
                           if (bcontribs(nd,b_contrib)) then
                               realint = realint + realint_b(nd,b_contrib,corr_on_beam)
                           endif
                       enddo
                   endif
               enddo
           enddo

           if (abs(realint*wgt) > wtmax) then
               wtmax = abs(realint*wgt)
           endif

      end function

end module

