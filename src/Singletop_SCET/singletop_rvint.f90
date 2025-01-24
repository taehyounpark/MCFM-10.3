!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module singletop_rvint
      use singletop_int
      implicit none

      private

      public :: rvint

      contains

#define ONLY_MAIN_B 1
#define WITH_BBAR 3

      function rvint_assemble_hxd(ndmax,bcontribs,msq,msqc,pswt,xx)
          implicit none
          include 'nf.f'
          include 'kpart.f'! coeffonly
          include 'constants.f'
          include 'energy.f'
          include 'ipsgen.f'
          include 'maxd.f'
          include 'incldip.f'
          include 'beamtype.f'
          real(dp) :: BrnRat
          common/BrnRat/BrnRat

          real(dp) :: rvint_assemble_hxd(0:ndmax,max_bcontrib)
          integer, intent(in) :: ndmax
          logical, intent(in) :: bcontribs(0:ndmax,max_bcontrib)
          real(dp), intent(in) :: msq(-nf:nf,-nf:nf,max_bcontrib)
          real(dp), intent(in) :: msqc(ndmax,-nf:nf,-nf:nf,max_bcontrib)
          real(dp), intent(in) :: pswt, xx(2)

          integer :: j,k,nd
          real(dp) :: flux

          rvint_assemble_hxd(:,:) = 0._dp

          do b_contrib=1,max_bcontrib
              do nd=0,ndmax
                  if ((bcontribs(nd,b_contrib) .eqv. .false.) .or. &
                      ((nd > 0) .and. (incldip(nd) .eqv. .false.)) ) cycle

                  do j=-nf,nf; do k=-nf,nf
                      if ( any(j == [2,4,-1,-3]) .and. any(k == [5,0]) ) then
                          if (nd == 0) then
                              rvint_assemble_hxd(nd,b_contrib) = rvint_assemble_hxd(nd,b_contrib) + &
                                   msq(j,k,b_contrib) * &
                                       fdist_one(ih1, xx(1), facscale_beam1_islight_onheavy, j, 1) * &
                                       fdist_one(ih2, xx(2), facscale_beam2_isheavy_onheavy, k, 2)
                          else
                              rvint_assemble_hxd(nd,b_contrib) = rvint_assemble_hxd(nd,b_contrib) - &
                                    msqc(nd,j,k,b_contrib) * &
                                        fdist_one(ih1, xx(1), facscale_beam1_islight_onheavy, j, 1) * &
                                        fdist_one(ih2, xx(2), singletop2_dipscale(nd,2), k, 2)
                          endif
                      elseif ( any(j == [5,0]) .and. any(k == [2,4,-1,-3]) ) then
                          if (nd == 0) then
                              rvint_assemble_hxd(nd,b_contrib) = rvint_assemble_hxd(nd,b_contrib) + &
                                    msq(j,k,b_contrib) * &
                                        fdist_one(ih1, xx(1), facscale_beam1_isheavy_onheavy, j, 1) * &
                                        fdist_one(ih2, xx(2), facscale_beam2_islight_onheavy, k, 2)
                          else
                              rvint_assemble_hxd(nd,b_contrib) = rvint_assemble_hxd(nd,b_contrib) - &
                                    msqc(nd,j,k,b_contrib) * &
                                        fdist_one(ih1, xx(1), singletop2_dipscale(nd,1), j, 1) * &
                                        fdist_one(ih2, xx(2), facscale_beam2_islight_onheavy, k, 2)
                          endif
                      endif
                  enddo; enddo
              enddo
          enddo

          flux = fbGev2/(2._dp*xx(1)*xx(2)*sqrts**2)
          rvint_assemble_hxd = rvint_assemble_hxd*flux*pswt/BrnRat

      end function rvint_assemble_hxd

      function rvint_assemble_lxd(ndmax,bcontribs,msq,msqc,pswt,xx)
          implicit none
          include 'nf.f'
          include 'kpart.f'! coeffonly
          include 'constants.f'
          include 'energy.f'
          include 'ipsgen.f'
          include 'maxd.f'
          include 'incldip.f'
          include 'beamtype.f'
          real(dp) :: BrnRat
          common/BrnRat/BrnRat

          real(dp) :: rvint_assemble_lxd(0:ndmax,max_bcontrib)
          integer, intent(in) :: ndmax
          logical, intent(in) :: bcontribs(0:ndmax,max_bcontrib)
          real(dp), intent(in) :: msq(-nf:nf,-nf:nf,max_bcontrib)
          real(dp), intent(in) :: msqc(ndmax,-nf:nf,-nf:nf,max_bcontrib)
          real(dp), intent(in) :: pswt, xx(2)

          integer :: j,k,nd
          real(dp) :: flux

          rvint_assemble_lxd(:,:) = 0._dp

          do b_contrib=1,max_bcontrib
              do nd=0,ndmax
                  if ((bcontribs(nd,b_contrib) .eqv. .false.) .or. &
                      ((nd > 0) .and. (incldip(nd) .eqv. .false.)) ) cycle

                  do j=-nf,nf; do k=-nf,nf
                      if ( any(j == [2,4,-1,-3,0]) .and. k == 5 ) then
                          if (nd == 0) then
                              rvint_assemble_lxd(nd,b_contrib) = rvint_assemble_lxd(nd,b_contrib) + &
                                   msq(j,k,b_contrib) * &
                                       fdist_one(ih1, xx(1), facscale_beam1_islight_onlight, j, 1) * &
                                       fdist_one(ih2, xx(2), facscale_beam2_isheavy_onlight, k, 2)
                          else
                              rvint_assemble_lxd(nd,b_contrib) = rvint_assemble_lxd(nd,b_contrib) - &
                                    msqc(nd,j,k,b_contrib) * &
                                        fdist_one(ih1, xx(1), singletop2_dipscale(nd,1), j, 1) * &
                                        fdist_one(ih2, xx(2), facscale_beam2_isheavy_onlight, k, 2)
                          endif
                      elseif (j == 5 .and. any(k == [2,4,-1,-3,0])) then
                          if (nd == 0) then
                              rvint_assemble_lxd(nd,b_contrib) = rvint_assemble_lxd(nd,b_contrib) + &
                                    msq(j,k,b_contrib) * &
                                        fdist_one(ih2, xx(2), facscale_beam2_islight_onlight, k, 2) * &
                                        fdist_one(ih1, xx(1), facscale_beam1_isheavy_onlight, j, 1)
                          else
                              rvint_assemble_lxd(nd,b_contrib) = rvint_assemble_lxd(nd,b_contrib) - &
                                    msqc(nd,j,k,b_contrib) * &
                                        fdist_one(ih2, xx(2), singletop2_dipscale(nd,2), k, 2) * &
                                        fdist_one(ih1, xx(1), facscale_beam1_isheavy_onlight, j, 1)
                          endif
                      endif
                  enddo; enddo
              enddo
          enddo

          flux = fbGev2/(2._dp*xx(1)*xx(2)*sqrts**2)
          rvint_assemble_lxd = rvint_assemble_lxd*flux*pswt/BrnRat

      end function rvint_assemble_lxd


      function rvint_assemble_lxh(ndmax,bcontribs,msq,msqc,z1,z2,pswt,xx)
          implicit none
          include 'nf.f'
          include 'kpart.f'! coeffonly
          include 'constants.f'
          include 'energy.f'
          include 'ipsgen.f'
          include 'maxd.f'
          include 'incldip.f'
          include 'beamtype.f'
          real(dp) :: BrnRat
          common/BrnRat/BrnRat

          real(dp) :: rvint_assemble_lxh(0:ndmax,max_bcontrib)
          integer, intent(in) :: ndmax
          logical, intent(in) :: bcontribs(0:ndmax,max_bcontrib)
          real(dp), intent(in) :: msq(-nf:nf,-nf:nf,max_bcontrib)
          real(dp), intent(in) :: msqc(ndmax,-nf:nf,-nf:nf,max_bcontrib)
          real(dp), intent(in) :: z1, z2, pswt, xx(2)

          integer :: j,k,nd
          real(dp) :: flux

          rvint_assemble_lxh(:,:) = 0._dp

          do b_contrib=1,max_bcontrib
              do nd=0,ndmax
                  if ((bcontribs(nd,b_contrib) .eqv. .false.) .or. &
                      ((nd > 0) .and. (incldip(nd) .eqv. .false.)) ) cycle

                  do j=-nf,nf; do k=-nf,nf
                      if ( any(j == [2,4,-1,-3]) .and. k == 5 ) then
                          if (nd == 0) then
                              rvint_assemble_lxh(nd,b_contrib) = rvint_assemble_lxh(nd,b_contrib) + &
                                   msq(j,k,b_contrib) * fdist_one(ih1, xx(1), facscale_beam1_islight_onlight, j, 1)
                          else
                              rvint_assemble_lxh(nd,b_contrib) = rvint_assemble_lxh(nd,b_contrib) - &
                                    msqc(nd,j,k,b_contrib) * fdist_one(ih1, xx(1), singletop2_dipscale(nd,1), j, 1)
                          endif
                      elseif (j == 5 .and. any(k == [2,4,-1,-3])) then
                          if (nd == 0) then
                              rvint_assemble_lxh(nd,b_contrib) = rvint_assemble_lxh(nd,b_contrib) + &
                                    msq(j,k,b_contrib) * fdist_one(ih2, xx(2), facscale_beam2_islight_onlight, k, 2)
                          else
                              rvint_assemble_lxh(nd,b_contrib) = rvint_assemble_lxh(nd,b_contrib) - &
                                    msqc(nd,j,k,b_contrib) * fdist_one(ih2, xx(2), singletop2_dipscale(nd,2), k, 2)
                          endif
                      elseif (j==0 .and. k==5) then
                          if (nd == 0) then
                              rvint_assemble_lxh(nd,b_contrib) = rvint_assemble_lxh(nd,b_contrib) + &
                                    msq(j,k,b_contrib) * fdist_one(ih1, xx(1), facscale_beam1_islight_onlight, j, 1)
                          else
                              rvint_assemble_lxh(nd,b_contrib) = rvint_assemble_lxh(nd,b_contrib) - &
                                    msqc(nd,j,k,b_contrib) * fdist_one(ih1, xx(1), singletop2_dipscale(nd,1), j, 1)
                          endif
                      elseif (j==5 .and. k==0) then
                          if (nd == 0) then
                              rvint_assemble_lxh(nd,b_contrib) = rvint_assemble_lxh(nd,b_contrib) + &
                                    msq(j,k,b_contrib) * fdist_one(ih2, xx(2), facscale_beam2_islight_onlight, k, 2)
                          else
                              rvint_assemble_lxh(nd,b_contrib) = rvint_assemble_lxh(nd,b_contrib) - &
                                    msqc(nd,j,k,b_contrib) * fdist_one(ih2, xx(2), singletop2_dipscale(nd,2), k, 2)
                          endif
                      endif
                  enddo; enddo
              enddo
          enddo

          flux = fbGev2/(2._dp*xx(1)*xx(2)*sqrts**2)
          rvint_assemble_lxh = rvint_assemble_lxh*(4*sqrt(z1*z2))*flux*pswt/BrnRat

      end function rvint_assemble_lxh

      function rvint(r,wgt)
          use Superhisto, only : shtmpreset, shtmpcommit
          use SafetyCuts, only : passed_smallnew
          use MCFMSetupPlots, only: nplotter_new
          use singletop_interf_lxh
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
          include 'taucut.f'! usescet
          include 'maxd.f'
          include 'incldip.f'
          include 'epinv.f'
          include 'epinv2.f'
          include 'beamtype.f'

          logical :: bin
          common/bin/bin

          real(dp) :: rvint
          real(dp), intent(in) :: r(mxdim), wgt

          ! external functions
          logical :: includedipole

          real(dp) :: p(mxpart,4), pjet(mxpart,4)
          real(dp) :: ptilde(mxpart,4)
          real(dp) :: pswt

          real(dp), allocatable :: realint_b(:,:), tmp_b(:,:)
          real(dp), allocatable :: scalereweight_b(:,:,:)
          real(dp), allocatable :: pdfreweight_b(:,:,:)
          real(dp), allocatable :: scetreweight_b(:,:,:)
          !real(dp), allocatable :: singletop2_dipscale_save(:,:)


          real(dp) :: msq(-nf:nf,-nf:nf,max_bcontrib)
          real(dp) :: msq_tmp(-nf:nf,-nf:nf,max_bcontrib)

          integer :: ndmax
          logical, allocatable :: bcontribs(:,:)

          real(dp), allocatable :: msqc(:,:,:,:)
          real(dp), allocatable :: msqc_tmp(:,:,:,:)

          integer :: j,k,nd,ipdf
          real(dp) :: taucut_save

          integer :: truescalevar

           p(:,:) = 0._dp
           pjet(:,:) = 0._dp
           currentPDF = 0

           epinv = 0._dp
           epinv2 = 0._dp
            
           b_contrib = 1


           ! SCET for all interference pieces
           if (currentContrib > 3) then
               usescet = .true.
           endif

           if (.not. gen_singletop(r,p,pswt)) then
               rvint = 0._dp
               return
           endif
 
           if (any(ieee_is_nan(p(1:npart+2,:)))) then
               write(6,*) 'Discarding NaN phase space point'
               rvint = 0._dp
               return
           endif
 
           if (.not. passed_smallnew(p,npart,1d-8)) then
               rvint = 0._dp
               return
           endif

 
           xx(1)=-2._dp*p(1,4)/sqrts
           xx(2)=-2._dp*p(2,4)/sqrts
           if ((xx(1) >  1._dp) .or. (xx(2) >  1._dp) &
               .or.(xx(1) < xmin)   .or. (xx(2) < xmin)) then
               rvint = 0._dp
               return
           endif

 
           if (currentContrib == 4) then
               ndmax = 8
               z1int = r(ndim-1)**2
               z2int = r(ndim)**2
           elseif (currentContrib == 5) then
               ndmax = 8
           elseif (currentContrib == 6) then
               ndmax = 6
           endif

           call singletop2_scale_setup(p)

           allocate(msqc(ndmax,-nf:nf,-nf:nf,max_bcontrib))
           allocate(msqc_tmp(ndmax,-nf:nf,-nf:nf,max_bcontrib))

           if (currentContrib == 4) then ! lxh
               call singletop_jet_light_heavy_rv(p,msq)
               call singletop_jet_light_heavy_rv_gs(p,ndmax,msqc)
           elseif (currentContrib == 5) then ! lxd
               call singletop_light_decay_rv(p,msq)
               call singletop_light_decay_rv_gs(p,ndmax,msqc)
           elseif (currentContrib == 6) then ! hxd
               call singletop_heavy_decay_rv(p,msq)
               call singletop_heavy_decay_rv_gs(p,ndmax,msqc)
           else
               error stop "unknown currentContrib in singletop_rvint"
           endif

           ! reset scales, which have been modified by ptrans
           call singletop2_scale_setup(p)

           ! now we can determine contributions that vanish
           allocate(bcontribs(0:ndmax,max_bcontrib), source=.false.)


           bcontribs(:,ONLY_MAIN_B) = .true.

           ! hxd only piece with bbar
           if (currentContrib == 6) then
               bcontribs(:,WITH_BBAR) = .true.
           endif

           call storeptilde(0,p)
           do nd=0,ndmax
              if ((nd > 0) .and. (incldip(nd) .eqv. .false.)) cycle
              call getptilde(nd,ptilde)
              do b_contrib=1,max_bcontrib
                  if (bcontribs(nd,b_contrib) .eqv. .false.) cycle

                  currentNd = nd ! to correctly fill jetcontent
                  ! cuts may depend on b_contrib, that's why we need to call includedipole for all pieces
                  bcontribs(nd,b_contrib) = includedipole(nd,ptilde)
              enddo
           enddo

           if (all(bcontribs .eqv. .false.)) then
               rvint = 0._dp
               return
           endif

           allocate(realint_b(0:ndmax,max_bcontrib), source = 0._dp)
           allocate(tmp_b(0:ndmax,max_bcontrib))
            
           ! WARNING!!
           ! need to reset scales before assembly, since dipole routines have
           ! overwritten scales that are accessed in assembly routine
           call singletop2_scale_setup(p)

           if (currentContrib == 4) then
               realint_b = rvint_assemble_lxh(ndmax,bcontribs,msq,msqc,z1int,z2int,pswt,xx)
           elseif (currentContrib == 5) then
               realint_b = rvint_assemble_lxd(ndmax,bcontribs,msq,msqc,pswt,xx)
           elseif (currentContrib == 6) then
               realint_b = rvint_assemble_hxd(ndmax,bcontribs,msq,msqc,pswt,xx)
           endif

           !call singcheck(singletop_jet_light_heavy_rv, &
               !singletop_jet_light_heavy_rv_gs, p)

           !call singcheck(singletop_light_decay_rv, &
               !singletop_light_decay_rv_gs, p)

           !call singcheck(singletop_heavy_decay_rv, &
               !singletop_heavy_decay_rv_gs, p)

           if (any(ieee_is_nan(realint_b)) .or. (.not. all(ieee_is_finite(realint_b)))) then
               !write(6,*) 'Discarding NaN matrix element!'
               rvint = 0._dp
               return
           endif

           if (bin .and. doMultitaucut) then
               allocate(scetreweight_b(size(tcutarray),0:ndmax,max_bcontrib))
               taucut_save = taucut
               do j=1,size(tcutarray)
                   taucut = tcutarray(j)
                   call singletop2_scale_setup(p)
                   if (currentContrib == 4) then
                       call singletop_jet_light_heavy_rv(p,msq_tmp)
                       call singletop_jet_light_heavy_rv_gs(p,ndmax,msqc_tmp)
                   elseif (currentContrib == 5) then
                       call singletop_light_decay_rv(p,msq_tmp)
                       call singletop_light_decay_rv_gs(p,ndmax,msqc_tmp)
                   elseif (currentContrib == 6) then
                       call singletop_heavy_decay_rv(p,msq_tmp)
                       call singletop_heavy_decay_rv_gs(p,ndmax,msqc_tmp)
                   endif
                   call singletop2_scale_setup(p)
                   if (currentContrib == 4) then
                       tmp_b = rvint_assemble_lxh(ndmax,bcontribs,msq_tmp,msqc_tmp,z1int,z2int,pswt,xx)
                   elseif (currentContrib == 5) then
                       tmp_b = rvint_assemble_lxd(ndmax,bcontribs,msq_tmp,msqc_tmp,pswt,xx)
                   elseif (currentContrib == 6) then
                       tmp_b = rvint_assemble_hxd(ndmax,bcontribs,msq_tmp,msqc_tmp,pswt,xx)
                   endif
    
                   do b_contrib=1,max_bcontrib
                       do nd=0,ndmax
                           if (bcontribs(nd,b_contrib) .eqv. .false.) cycle

                           if (realint_b(nd,b_contrib) /= 0._dp) then
                               scetreweight_b(j,nd,b_contrib) = tmp_b(nd,b_contrib) / realint_b(nd,b_contrib)
                           else
                               scetreweight_b(:,nd,b_contrib) = 1._dp
                           endif
                       enddo
                   enddo
               enddo

               taucut = taucut_save
           endif

           if (doScalevar .and. bin) then
               allocate(scalereweight_b(0:ndmax,maxscalevar,max_bcontrib), source=0._dp)

               !allocate(singletop2_dipscale_save(ndmax,2),source=0._dp)
               !singletop2_dipscale_save(1:ndmax,:) = singletop2_dipscale(1:ndmax,:)

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

                   if (currentContrib == 4) then
                       call singletop_jet_light_heavy_rv(p,msq_tmp)
                       call singletop_jet_light_heavy_rv_gs(p,ndmax,msqc_tmp)
                   elseif (currentContrib == 5) then
                       call singletop_light_decay_rv(p,msq_tmp)
                       call singletop_light_decay_rv_gs(p,ndmax,msqc_tmp)
                   elseif (currentContrib == 6) then
                       call singletop_heavy_decay_rv(p,msq_tmp)
                       call singletop_heavy_decay_rv_gs(p,ndmax,msqc_tmp)
                   endif

                   if (j <= truescalevar) then
                       call singletop2_scale_setup(p, mult_in_ren=scalevarmult(j), mult_in_fac=facscalevarmult(j))
                   elseif (j == (truescalevar + 1)) then
                       call singletop2_scale_setup(p, mult_in_ren=1d0, mult_in_fac=1d0, forcemt=.true.)
                   else
                       call singletop2_scale_setup(p, mult_in_ren=scalevarmult(j-truescalevar-1), &
                           mult_in_fac=facscalevarmult(j-truescalevar-1), forcemt=.true.)
                   endif

                   if (currentContrib == 4) then
                       tmp_b = rvint_assemble_lxh(ndmax,bcontribs,msq_tmp,msqc_tmp,z1int,z2int,pswt,xx)
                   elseif (currentContrib == 5) then
                       tmp_b = rvint_assemble_lxd(ndmax,bcontribs,msq_tmp,msqc_tmp,pswt,xx)
                   elseif (currentContrib == 6) then
                       tmp_b = rvint_assemble_hxd(ndmax,bcontribs,msq_tmp,msqc_tmp,pswt,xx)
                   endif

                   do b_contrib=1,max_bcontrib
                       do k=0,ndmax
                           if (realint_b(k,b_contrib) == 0._dp) cycle
                           if (bcontribs(k,b_contrib) .eqv. .false.) cycle
                           if ((k > 0) .and. (incldip(k) .eqv. .false.)) cycle
                           scalereweight_b(k, j, b_contrib) = tmp_b(k,b_contrib) / realint_b(k,b_contrib)
                       enddo
                   enddo

               enddo

               ! reset old scales
               call singletop2_scale_setup(p)
               !restore old dipscales
               !singletop2_dipscale(1:ndmax,:) = singletop2_dipscale_save(1:ndmax,:)
           endif


           if (maxPDFsets > 0 .and. bin) then
               allocate(pdfreweight_b(0:ndmax,maxPDFsets,max_bcontrib), source=0._dp)

               do ipdf=1,maxPDFsets
                   currentPDF = ipdf

                   call singletop2_scale_setup(p)
                   if (currentContrib == 4) then
                       call singletop_jet_light_heavy_rv(p,msq_tmp)
                       call singletop_jet_light_heavy_rv_gs(p,ndmax,msqc_tmp)
                   elseif (currentContrib == 5) then
                       call singletop_light_decay_rv(p,msq_tmp)
                       call singletop_light_decay_rv_gs(p,ndmax,msqc_tmp)
                   elseif (currentContrib == 6) then
                       call singletop_heavy_decay_rv(p,msq_tmp)
                       call singletop_heavy_decay_rv_gs(p,ndmax,msqc_tmp)
                   endif
                   call singletop2_scale_setup(p)

                   if (currentContrib == 4) then
                       tmp_b = rvint_assemble_lxh(ndmax,bcontribs,msq_tmp,msqc_tmp,z1int,z2int,pswt,xx)
                   elseif (currentContrib == 5) then
                       tmp_b = rvint_assemble_lxd(ndmax,bcontribs,msq_tmp,msqc_tmp,pswt,xx)
                   elseif (currentContrib == 6) then
                       tmp_b = rvint_assemble_hxd(ndmax,bcontribs,msq_tmp,msqc_tmp,pswt,xx)
                   endif

                   do b_contrib=1,max_bcontrib
                       do nd=0,ndmax
                           if (bcontribs(nd,b_contrib) .eqv. .false.) cycle
                           if ((nd > 0) .and. (incldip(nd) .eqv. .false.)) cycle
                           
                           pdfreweight_b(nd,ipdf,b_contrib) = &
                               (realint_b(nd,b_contrib) - tmp_b(nd,b_contrib))*wgt
                       enddo
                   enddo
               enddo

               currentPDF = 0
           endif

!!! BINNING
           if (bin) then
               do b_contrib=1,max_bcontrib
                   if (all(bcontribs(:,b_contrib) .eqv. .false.)) cycle

                   do nd=0,ndmax
                       if (bcontribs(nd,b_contrib) .eqv. .false.) cycle

                       ! incldip gets overwritten after scale variation and
                       ! pdf uncertainties, but *should* stay the same.
                       if ((nd > 0) .and. (incldip(nd) .eqv. .false.)) cycle

                       currentNd = nd
                       call getptildejet(nd,pjet)
                       includeTaucutgrid(currentNd) = .true.
           
                       if (doMultitaucut) then
                          scetreweight(:) = scetreweight_b(:,nd,b_contrib)
                       endif

                       if (doScalevar) then
                           scalereweight(1:maxscalevar) = scalereweight_b(nd,1:maxscalevar,b_contrib)
                       endif

                       if (maxPDFsets > 0) then
                           pdfreweight(:) = pdfreweight_b(nd,:,b_contrib)
                       endif

                       call nplotter_new(pjet,realint_b(nd,b_contrib)*wgt)

                   enddo
               enddo

               call threadStorageOp(shtmpcommit)
           endif

!! finalization for vegas

           rvint = sum(realint_b)

           if (abs(rvint*wgt) > wtmax) then
               wtmax = abs(rvint*wgt)
           endif

      end function

end module

