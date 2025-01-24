!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
module singletop_virtint
      use singletop_int
      implicit none

      private

      public :: virtint

#define PRINT_EPINV 0

      contains

      ! mostly copied from lowint, with added dipole pieces
      function virtint(r,wgt)
          use Superhisto, only : shtmpreset, shtmpcommit
          use SafetyCuts, only : passed_smallnew
          use MCFMSetupPlots, only: nplotter_new
          use singletop_vvint, only: vvint
          implicit none
          include 'nf.f'
          include 'vegas_common.f'! ndim, mxdim
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
          include 'epinv.f'
          include 'epinv2.f'
          include 'taucut.f'! usescet
          include 'beamtype.f'

          logical :: bin
          common/bin/bin
          real(dp) :: BrnRat
          common/BrnRat/BrnRat

          real(dp) :: virtint 
          real(dp), intent(in) :: r(mxdim), wgt

          ! external functions
          logical :: includedipole

          real(dp) :: p(mxpart,4), pjet(mxpart,4)
          real(dp) :: pswt

          integer :: kcasemap
          ! ipsgen=1,2,3, kcase
          logical, save :: epinv_checked(6,2) = .false.
!$omp threadprivate(epinv_checked)

          real(dp) :: virtint_b(max_bcontrib,max_corr_on_beam), &
              tmp, tmp_b(max_bcontrib,max_corr_on_beam)

          real(dp), allocatable :: scalereweight_b(:,:,:)
          real(dp), allocatable :: pdfreweight_b(:,:,:)
          real(dp), allocatable :: scetreweight_b(:,:)

          logical :: passed_taucut(max_corr_on_beam)

          real(dp) :: msqall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)
          real(dp) :: msqvall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)
          real(dp) :: msqvall_tmp(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

          integer :: truescalevar


          real(dp) :: z

          logical :: bcontribs(max_bcontrib)

          integer :: j,m,ipdf

          ! shortcut for interference
          if ((currentContrib == 4 .and. ipsgen == 5) .or. &
              (currentContrib == 5 .and. ipsgen == 7) .or. &
              (currentContrib == 6 .and. ipsgen == 9)) then
              virtint = vvint(r,wgt)
              return
          endif

          virtint = 0._dp
          virtint_b = 0._dp
          p(:,:) = 0._dp
          pjet(:,:) = 0._dp

          ! SCET for all interference pieces
          if (currentContrib > 3) then
              usescet = .true.
          endif

          if (.not. gen_singletop(r,p,pswt)) then
              virtint = 0._dp
              return
          endif

          if (all(.not. ieee_is_nan(p(1:npart,:))) .eqv. .false.) then
              if (debug) then
                   write(6,*) 'Discarding NaN or infinite phase space point'
              endif
              virtint = 0._dp
              return
          endif

          if (.not. passed_smallnew(p,npart,1d-12)) then
              virtint = 0._dp
              return
          endif

          ! additional tau cutoff for one additional radiation
          if (kcase == kbq_tpq_jet) then
              if (currentContrib == 1) then
                  do corr_on_beam=1,2
                      if (.not. passed_taucut_light(p, taucut_in=1d-6)) then
                          virtint= 0._dp
                          return
                      endif
                  enddo
              elseif (currentContrib == 2) then
                  do corr_on_beam=1,2
                      if (.not. passed_taucut_heavyprod(p, taucut_in=1d-6)) then
                          virtint= 0._dp
                          return
                      endif
                  enddo
              elseif (currentContrib == 3) then
                  ! ...
              endif
          endif


          xx(1)=-2._dp*p(1,4)/sqrts
          xx(2)=-2._dp*p(2,4)/sqrts
          if ((xx(1) >  1._dp) .or. (xx(2) >  1._dp) &
              .or.(xx(1) < xmin)   .or. (xx(2) < xmin)) then
              virtint = 0._dp
              return
          endif

          z=r(ndim)**2
          ! taken into account in assembly routine
          !xjac=two*sqrt(z)
      
          ! things that only need to be set up once for all pieces

          bcontribs = .false.

          ! bcontribs are
          ! 1: main b
          ! 2: additional pair b~ b~
          ! 3: additional b~
          ! 4: additional pair b b~

          ! set up all possible b contributions
          if (kcase == kbq_tpq) then
              bcontribs(1) = .true.
              if (currentContrib == 2) then
                  bcontribs(3) = .true.
              elseif (currentContrib == 4) then
                  bcontribs(3) = .true.
              endif
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
                      elseif (currentContrib == 4 .and. ipsgen == 4) then
                          passed_taucut(corr_on_beam) = &
                              passed_taucut_lxh(p,scetreweight_b(:,corr_on_beam))
                      elseif (currentContrib == 5 .and. ipsgen == 6) then
                          passed_taucut(corr_on_beam) = &
                              passed_taucut_decay(p,scetreweight_b(:,corr_on_beam))
                      elseif (currentContrib == 6 .and. ipsgen == 8) then
                          passed_taucut(corr_on_beam) = &
                              passed_taucut_decay(p,scetreweight_b(:,corr_on_beam))
                      endif
                  enddo
              endif

              currentNd = 0
              bcontribs(b_contrib) = includedipole(0,p)
          enddo

          if (all(bcontribs .eqv. .false.)) then
              virtint = 0._dp
              return
          endif

!!! central value
          ! set all DDIS singletop scale variables based on p
          ! or as constant or dynamic if use_DDIS is .false.
          call singletop2_scale_setup(p)
          ! fill pdf ingredients
          call calc_singletop_pdfs_virt(xx,z)
          ! fill msq ingredients
          call lowint_fillmsq(p, bcontribs, msqall)

          if (kcase == kbq_tpq) then 
              kcasemap = 1
          else
              kcasemap = 2
          endif

          if (epinv_checked(currentContrib,kcasemap) .eqv. .false.) then
              ! gives 1/eps piece
              epinv = 10._dp
              epinv2 = 10._dp

              call virtint_fillmsq(p, z, bcontribs, msqvall_tmp)
              tmp_b = virtint_assemble(msqall,msqvall_tmp,z,bcontribs,pswt,xx)
              tmp = sum(tmp_b)
          endif

          epinv = 0._dp
          epinv2 = 0._dp
          call virtint_fillmsq(p, z, bcontribs, msqvall)
          virtint_b = virtint_assemble(msqall,msqvall,z,bcontribs,pswt,xx)

          virtint = sum(virtint_b)

          if (epinv_checked(currentContrib,kcasemap) .eqv. .false.) then
#if PRINT_EPINV == 1
!$omp critical(epinv_print)
              if ( virtint == 0._dp .and. tmp == 0._dp ) then
                  write (*,*) "INFO: virtint zero for this phase space point"
                  write (*,*) "Delaying epinv check"
              elseif ( abs(virtint/tmp - 1d0) > 1d-8 .or. &
                      ieee_is_nan(virtint/tmp) ) then
                  write (*,'(A,I1,A,I1)') "ERROR: poles do not cancel"
                  write (*,*) "results: ", virtint, tmp
                  write (*,*) "ratio = ", virtint/tmp
                  write (*,*) "difference 25 to born", (msqvall(2,5,1,2) - &
                  msqvall_tmp(2,5,1,2)) / msqall(2,5,1,2)
              else
#endif
                  epinv_checked(currentContrib,kcasemap) = .true.
#if PRINT_EPINV == 1
                  write (*,'(A,E7.1)') "Cancelation of poles checked to precision = ", &
                      abs(virtint/tmp- 1d0)
              endif
!$omp end critical(epinv_print)
#endif
          endif

           if (any(ieee_is_nan(virtint_b)) .or. (.not. all(ieee_is_finite(virtint_b)))) then
               !write(6,*) 'Discarding NaN matrix element!'
               virtint = 0._dp
               return
           endif


!!! pdf uncertainties
! reuses msqall, msqvall, have it before scale variation
          if (maxPDFsets > 0 .and. bin) then
              allocate(pdfreweight_b(maxPDFsets,max_bcontrib,max_corr_on_beam), source=0._dp)

              do ipdf=1,maxPDFsets
                  currentPDF = ipdf
                  call calc_singletop_pdfs_virt(xx, z)
                  if (doPDFAlphas) then
                      call singletop2_scale_setup(p)
                      call lowint_fillmsq(p, bcontribs, msqall)
                      call virtint_fillmsq(p, z, bcontribs, msqvall)
                  endif
                  tmp_b = virtint_assemble(msqall,msqvall,z,bcontribs,pswt,xx)
                  do b_contrib=1,max_bcontrib
                      if (bcontribs(b_contrib) .eqv. .false.) cycle
                      do m=1,maxbeams
                          corr_on_beam = beams_enabled(m)
                          pdfreweight_b(currentPDF, b_contrib, corr_on_beam) = &
                              (virtint_b(b_contrib,corr_on_beam) - tmp_b(b_contrib,corr_on_beam))*wgt
                      enddo
                  enddo
              enddo

              currentPDF = 0
          endif

!!! scale variation

          if (doScalevar .and. bin) then
              allocate(scalereweight_b(maxscalevar, max_bcontrib, max_corr_on_beam), source=0._dp)

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
                  call calc_singletop_pdfs_virt(xx,z)
                  call lowint_fillmsq(p, bcontribs, msqall)
                  call virtint_fillmsq(p, z, bcontribs, msqvall)
                  tmp_b = virtint_assemble(msqall,msqvall,z,bcontribs,pswt,xx)

                  do b_contrib=1,max_bcontrib
                      if (bcontribs(b_contrib) .eqv. .false.) cycle
                      do m=1,maxbeams
                          corr_on_beam = beams_enabled(m)
                          if (virtint_b(b_contrib,corr_on_beam) == 0._dp) cycle
                          scalereweight_b(j, b_contrib, corr_on_beam) = &
                              tmp_b(b_contrib,corr_on_beam) / virtint_b(b_contrib,corr_on_beam)
                      enddo
                  enddo
              enddo

              ! reset to unmodified scales
              call singletop2_scale_setup(p)
          endif


!!! final plotting and result for vegas

!         !where (ieee_is_nan(virtint_b)) virtint_b = 0._dp

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
                          scalereweight(1:maxscalevar) = &
                              scalereweight_b(1:maxscalevar,b_contrib,corr_on_beam)
                      endif

                      if (maxPDFsets > 0) then
                          pdfreweight(:) = pdfreweight_b(:,b_contrib,corr_on_beam)
                      endif

                      call nplotter_new(pjet,virtint_b(b_contrib,corr_on_beam)*wgt)
                  enddo
              enddo
              call threadStorageOp(shtmpcommit)
          endif

          virtint = 0._dp
          do m=1,maxbeams
              corr_on_beam = beams_enabled(m)
              if (passed_taucut(corr_on_beam)) then
                  do b_contrib=1,max_bcontrib
                      if (bcontribs(b_contrib) .eqv. .false.) cycle
                      virtint = virtint + virtint_b(b_contrib,corr_on_beam)
                  enddo
              endif
          enddo

          if (abs(virtint*wgt) > wtmax) then
              wtmax = abs(virtint*wgt)
          endif

          return
      end function

end module
