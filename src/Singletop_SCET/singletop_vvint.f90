!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module singletop_vvint
      use singletop_int
      implicit none

      private

      public :: vvint 

#define PRINT_EPINV 0

      contains

       function vvint_assemble_hxd(msq,msqv,z,pswt,xx)
          implicit none
          include 'nf.f'
          include 'kpart.f'! coeffonly
          include 'constants.f'
          include 'energy.f'
          include 'PR_new.f'
          include 'PR_stop.f'
          include 'agq.f'
          include 'epinv.f'
          include 'ipsgen.f'
          include 'beamtype.f'
          real(dp) :: BrnRat
          common/BrnRat/BrnRat

          real(dp) :: vvint_assemble_hxd
          real(dp), intent(in) :: msq(-nf:nf,-nf:nf)
          real(dp), intent(in) :: msqv(-nf:nf,-nf:nf)
          real(dp), intent(in) :: z, pswt, xx(2)

          real(dp) :: AP(-1:1, -1:1, 3)
          integer :: j,k
          real(dp) :: flux
          real(dp) :: xmsq, pdf

          xmsq = 0._dp

          do j=-nf,nf; do k=-nf,nf
              if (any(j == [2,4,-1,-3]) .and. k == 5) then
                  xmsq = xmsq + msqv(j,k) ! comes with PDfs

                  call singletop2_fillAP(z, i_beam2_heavy, AP)
                  xmsq = xmsq + msq(j,k)*(AP(q,q,1)-AP(q,q,3)+B2(b,b,q,1)-B2(b,b,q,3)) * &
                      fdist_one(ih1, xx(1), facscale_beam1_islight_onheavy, j, 1) * &
                      fdist_one(ih2, xx(2), facscale_beam2_isheavy_onheavy, k, 2)
                  xmsq = xmsq + msq(j,k)*(AP(q,q,2)+AP(q,q,3)+B2(b,b,q,2)+B2(b,b,q,3)) * &
                      fdist_one(ih1, xx(1), facscale_beam1_islight_onheavy, j, 1) * &
                      fdist_one(ih2, xx(2)/z, facscale_beam2_isheavy_onheavy, k, 2) / z

              elseif (j == 5 .and. any(k == [2,4,-1,-3])) then
                  xmsq = xmsq + msqv(j,k)

                  call singletop2_fillAP(z, i_beam1_heavy, AP)
                  xmsq = xmsq + msq(j,k)*(AP(q,q,1)-AP(q,q,3)+B1(b,b,q,1)-B1(b,b,q,3)) * &
                      fdist_one(ih1, xx(1), facscale_beam1_isheavy_onheavy, j, 1) * &
                      fdist_one(ih2, xx(2), facscale_beam2_islight_onheavy, k, 2)
                  xmsq = xmsq + msq(j,k)*(AP(q,q,2)+AP(q,q,3)+B1(b,b,q,2)+B1(b,b,q,3)) * &
                      fdist_one(ih1, xx(1)/z, facscale_beam1_isheavy_onheavy, j, 1) / z * &
                      fdist_one(ih2, xx(2), facscale_beam2_islight_onheavy, k, 2)
              elseif (any(j == [2,4,-1,-3]) .and. k == 0) then
                  call singletop2_fillAP(z, i_beam2_heavy, AP)
                  xmsq = xmsq + msq(j,+5)*(AP(q,g,2) + B2(q,g,q,2)) * &
                      fdist_one(ih1, xx(1), facscale_beam1_islight_onheavy, j, 1) * &
                      fdist_one(ih2, xx(2)/z, facscale_beam2_isheavy_onheavy, k, 2) / z
              elseif (j == 0 .and. any(k == [2,4,-1,-3])) then
                  call singletop2_fillAP(z, i_beam1_heavy, AP)
                  xmsq = xmsq + msq(+5,k)*(AP(q,g,2) + B1(q,g,q,2)) * &
                      fdist_one(ih1, xx(1)/z, facscale_beam1_isheavy_onheavy, j, 1) / z * &
                      fdist_one(ih2, xx(2), facscale_beam2_islight_onheavy, k, 2)
              endif
          enddo; enddo
                  
          flux = fbGev2/(2._dp*xx(1)*xx(2)*sqrts**2)
          vvint_assemble_hxd = xmsq*(2*sqrt(z))*flux*pswt/BrnRat

       end function vvint_assemble_hxd

       function vvint_assemble_lxd(msq,msqv,z,pswt,xx)
          implicit none
          include 'nf.f'
          include 'kpart.f'! coeffonly
          include 'constants.f'
          include 'energy.f'
          include 'PR_new.f'
          include 'PR_stop.f'
          include 'agq.f'
          include 'epinv.f'
          include 'ipsgen.f'
          include 'beamtype.f'
          real(dp) :: BrnRat
          common/BrnRat/BrnRat

          real(dp) :: vvint_assemble_lxd
          real(dp), intent(in) :: msq(-nf:nf,-nf:nf)
          real(dp), intent(in) :: msqv(-nf:nf,-nf:nf)
          real(dp), intent(in) :: z, pswt, xx(2)

          real(dp) :: AP(-1:1, -1:1, 3)
          integer :: j,k
          real(dp) :: flux
          real(dp) :: xmsq, pdf

          xmsq = 0._dp

                            !fdist_one(ih1, xx(1), facscale_beam1_islight_onlight, j, 1)

          do j=-nf,nf; do k=-nf,nf
              if ( any(j == [2,4,-1,-3]) .and. k == 5 ) then
                  ! pdfs included
                  xmsq = xmsq + msqv(j,k)

                  call singletop2_fillAP(z, i_beam1_light, AP)
                  xmsq = xmsq + msq(j,k)*(AP(q,q,1)-AP(q,q,3)+B1(q,q,b,1)-B1(q,q,b,3)) * &
                          fdist_one(ih1, xx(1), facscale_beam1_islight_onlight, j, 1) * &
                          fdist_one(ih2, xx(2), facscale_beam2_isheavy_onlight, k, 2)
                  xmsq = xmsq + msq(j,k)*(AP(q,q,2)+AP(q,q,3)+B1(q,q,b,2)+B1(q,q,b,3)) * &
                          fdist_one(ih1, xx(1)/z, facscale_beam1_islight_onlight, j, 1)/z * &
                          fdist_one(ih2, xx(2), facscale_beam2_isheavy_onlight, k, 2)
              elseif (j == 5 .and. any(k == [2,4,-1,-3])) then
                  ! pdfs included
                  xmsq = xmsq + msqv(j,k)  !* singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light,k)

                  call singletop2_fillAP(z, i_beam2_light, AP)
                  xmsq = xmsq + msq(j,k)*(AP(q,q,1)-AP(q,q,3)+B2(q,q,b,1)-B2(q,q,b,3)) * &
                          fdist_one(ih1, xx(1), facscale_beam1_isheavy_onlight, j, 1) * &
                          fdist_one(ih2, xx(2), facscale_beam2_islight_onlight, k, 2)
                  xmsq = xmsq + msq(j,k)*(AP(q,q,2)+AP(q,q,3)+B2(q,q,b,2)+B2(q,q,b,3)) * &
                          fdist_one(ih1, xx(1), facscale_beam1_isheavy_onlight, j, 1) * &
                          fdist_one(ih2, xx(2)/z, facscale_beam2_islight_onlight, k, 2)/z
              elseif (j==5 .and. k==0) then
                  call singletop2_fillAP(z, i_beam2_light, AP)
                  xmsq = xmsq + sum(msq(j,[2,4,-1,-3])) * (AP(q,g,2) + Q2(q,g,q,2)) * &
                          fdist_one(ih1, xx(1), facscale_beam1_isheavy_onlight, j, 1) * &
                          fdist_one(ih2, xx(2)/z, facscale_beam2_islight_onlight, k, 2)/z
              elseif (j==0 .and. k==5) then
                  call singletop2_fillAP(z, i_beam1_light, AP)
                  xmsq = xmsq + sum(msq([2,4,-1,-3],k)) * (AP(q,g,2) + Q1(q,g,q,2)) * &
                      fdist_one(ih1, xx(1)/z, facscale_beam1_islight_onlight, j, 1)/z * &
                      fdist_one(ih2, xx(2), facscale_beam2_isheavy_onlight, k, 2)
              endif
          enddo; enddo
                  
          flux = fbGev2/(2._dp*xx(1)*xx(2)*sqrts**2)
          vvint_assemble_lxd = xmsq*(2*sqrt(z))*flux*pswt/BrnRat

       end function vvint_assemble_lxd

      function vvint_assemble_lxh(msq,msqv,z,z1,z2,pswt,xx)
          implicit none
          include 'nf.f'
          include 'kpart.f'! coeffonly
          include 'constants.f'
          include 'energy.f'
          include 'PR_new.f'
          include 'PR_stop.f'
          include 'agq.f'
          include 'epinv.f'
          include 'ipsgen.f'
          include 'beamtype.f'
          real(dp) :: BrnRat
          common/BrnRat/BrnRat

          real(dp) :: vvint_assemble_lxh
          real(dp), intent(in) :: msq(-nf:nf,-nf:nf)
          real(dp), intent(in) :: msqv(-nf:nf,-nf:nf)
          real(dp), intent(in) :: z, z1, z2, pswt, xx(2)

          real(dp) :: AP(-1:1, -1:1, 3)
          integer :: j,k
          real(dp) :: flux
          real(dp) :: xmsq, pdf

          xmsq = 0._dp

          do j=-nf,nf; do k=-nf,nf
              if ( any(j == [2,4,-1,-3]) .and. k == 5 ) then
                  ! comes pre-dressed with PDFs
                  xmsq = xmsq + msqv(j,k)
                  !* &
                  !     fdist_one(ih1, xx(1), facscale_beam1_islight_onlight, j, 1)

                  ! splitting term matrix elements come pre-dressed with heavy-line beam functions
                  call singletop2_fillAP(z, i_beam1_light, AP)

                  xmsq = xmsq + msq(j,k)*(AP(q,q,1)-AP(q,q,3)+B1(q,q,b,1)-B1(q,q,b,3)) * &
                            fdist_one(ih1, xx(1), facscale_beam1_islight_onlight, j, 1)

                  xmsq = xmsq + msq(j,k)*(AP(q,q,2)+AP(q,q,3)+B1(q,q,b,2)+B1(q,q,b,3)) * &
                            fdist_one(ih1, xx(1)/z, facscale_beam1_islight_onlight, j, 1) / z

              elseif (j == 5 .and. any(k == [2,4,-1,-3])) then
                  xmsq = xmsq + msqv(j,k)

                  call singletop2_fillAP(z, i_beam2_light, AP)
                  xmsq = xmsq + msq(j,k)*(AP(q,q,1)-AP(q,q,3)+B2(q,q,b,1)-B2(q,q,b,3)) * &
                      fdist_one(ih2,xx(2),facscale_beam2_islight_onlight,k,2)

                  xmsq = xmsq + msq(j,k)*(AP(q,q,2)+AP(q,q,3)+B2(q,q,b,2)+B2(q,q,b,3)) * &
                      fdist_one(ih2,xx(2)/z,facscale_beam2_islight_onlight,k,2) / z
              elseif (j==5 .and. k==0) then
                  call singletop2_fillAP(z, i_beam2_light, AP)
                  xmsq = xmsq + sum(msq(j,[2,4,-1,-3])) * (AP(q,g,2) + Q2(q,g,q,2)) * &
                         fdist_one(ih2,xx(2)/z,facscale_beam2_islight_onlight,k,2) / z
              elseif (j==0 .and. k==5) then
                  call singletop2_fillAP(z, i_beam1_light, AP)
                  xmsq = xmsq + sum(msq([2,4,-1,-3],k)) * (AP(q,g,2) + Q1(q,g,q,2)) * &
                         fdist_one(ih1,xx(1)/z,facscale_beam1_islight_onlight,j,1) / z
              endif
            
          enddo; enddo
                  
          flux = fbGev2/(2._dp*xx(1)*xx(2)*sqrts**2)
          vvint_assemble_lxh = xmsq*(4*sqrt(z1*z2))*(2*sqrt(z))*flux*pswt/BrnRat

      end function

      ! mostly copied from lowint, with added dipole pieces
      function vvint(r,wgt)
          use Superhisto, only : shtmpreset, shtmpcommit
          use SafetyCuts, only : passed_smallnew
          use MCFMSetupPlots, only: nplotter_new
          use singletop_interf_lxh ! z1int, z2int
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

          real(dp) :: vvint, tmp
          real(dp), intent(in) :: r(mxdim), wgt

          ! external functions
          logical :: includedipole

          real(dp) :: p(mxpart,4), pjet(mxpart,4)
          real(dp) :: msq(-nf:nf,-nf:nf), msqv(-nf:nf,-nf:nf), pswt
          real(dp) :: msq_tmp(-nf:nf,-nf:nf), msqv_tmp(-nf:nf,-nf:nf)

          ! ipsgen=1,2,3, kcase
          logical, save :: epinv_checked(6) = .false.
!$omp threadprivate(epinv_checked)

          real(dp) :: z

          real(dp) :: taucut_save

          integer :: truescalevar

          !logical :: bcontribs(max_bcontrib)

          integer :: j,ipdf

          vvint = 0._dp
          p(:,:) = 0._dp
          pjet(:,:) = 0._dp

          currentPDF = 0

          ! SCET for all interference pieces
          if (currentContrib > 3) then
              usescet = .true.
          else
              error stop "routine only for VxV interference"
          endif

          if (.not. gen_singletop(r,p,pswt)) then
              vvint = 0._dp
              return
          endif

          if (all(.not. ieee_is_nan(p(1:npart,:))) .eqv. .false.) then
          ! hard function without as/4/pi
              if (debug) then
                   write(6,*) 'Discarding NaN or infinite phase space point'
              endif
              vvint = 0._dp
              return
          endif

          if (.not. passed_smallnew(p,npart,1d-10)) then
              vvint = 0._dp
              return
          endif

          xx(1)=-2._dp*p(1,4)/sqrts
          xx(2)=-2._dp*p(2,4)/sqrts
          if ((xx(1) >  1._dp) .or. (xx(2) >  1._dp) &
              .or.(xx(1) < xmin)   .or. (xx(2) < xmin)) then
              vvint = 0._dp
              return
          endif

          if (currentContrib == 4) then ! lxh
              z=r(ndim-2)**2
              ! set as module level variable in module singletop_interf_lxh
              ! for use in below cut routine
              z1int = r(ndim-1)**2
              z2int = r(ndim)**2
          elseif (currentContrib == 5) then ! lxd
              z = r(ndim)**2
          elseif (currentContrib == 6) then ! hxd
              z = r(ndim)**2
          else
              error stop "undefined currentContrib in singletop_vvint"
          endif

          ! only main B for all contributions
          b_contrib = 1

          call singletop2_scale_setup(p)

          if (epinv_checked(currentContrib) .eqv. .false.) then
              ! gives 1/eps piece
              epinv = 10._dp
              epinv2 = 10._dp

              if (currentContrib == 4) then
                  call singletop_jet_light_heavy_vv_tree(p,msq_tmp)
                  call singletop_jet_light_heavy_vv(p,msqv_tmp)
                  call singletop_jet_light_heavy_vr_z(p,z)
                  tmp = vvint_assemble_lxh(msq_tmp,msqv_tmp,z,z1int,z2int,pswt,xx)
              elseif (currentContrib == 5) then
                  call singletop_light_decay_vv_tree(p,msq_tmp)
                  call singletop_light_decay_vv(p,msqv_tmp)
                  call singletop_light_decay_vr_z(p,z)
                  tmp = vvint_assemble_lxd(msq_tmp,msqv_tmp,z,pswt,xx)
              elseif (currentContrib == 6) then
                  call singletop_heavy_decay_vv_tree(p,msq_tmp)
                  call singletop_heavy_decay_vv(p,msqv_tmp)
                  call singletop_heavy_decay_vr_z(p,z)
                  tmp = vvint_assemble_hxd(msq_tmp,msqv_tmp,z,pswt,xx)
              else
                  error stop "undefined currentContrib in singletop_vvint"
              endif
          endif

          epinv = 0._dp
          epinv2 = 0._dp


          if (currentContrib == 4) then
              call singletop_jet_light_heavy_vv_tree(p,msq)
              call singletop_jet_light_heavy_vv(p,msqv)
              call singletop_jet_light_heavy_vr_z(p,z)
              vvint = vvint_assemble_lxh(msq,msqv,z,z1int,z2int,pswt,xx)
          elseif (currentContrib == 5) then
              call singletop_light_decay_vv_tree(p,msq)
              call singletop_light_decay_vv(p,msqv)
              call singletop_light_decay_vr_z(p,z)
              vvint = vvint_assemble_lxd(msq,msqv,z,pswt,xx)
          elseif (currentContrib == 6) then
              call singletop_heavy_decay_vv_tree(p,msq)
              call singletop_heavy_decay_vv(p,msqv)
              call singletop_heavy_decay_vr_z(p,z)
              vvint = vvint_assemble_hxd(msq,msqv,z,pswt,xx)
          else
              error stop "undefined currentContrib in singletop_vvint"
          endif

!         QB(1) = -2*p(1,4)
!         QB(2) = -2*p(2,4)
!         call lumxmsq_singletop_prod(p,xx,z1int,z2int,QB,1,xmsq,.true.)
!         write (*,*) msq(2,5) * fdist_one(ih1, xx(1), facscale_beam1_islight_onlight, 2, 1)/ vvint, &
!               vvint / xmsq

          if (epinv_checked(currentContrib) .eqv. .false.) then
#if PRINT_EPINV == 1
!$omp critical(epinv_print)
              if ( vvint == 0._dp .and. tmp == 0._dp ) then
                  write (*,*) "INFO: virtint zero for this phase space point"
                  write (*,*) "Delaying epinv check"
              elseif ( abs(vvint/tmp - 1d0) > 1d-8 .or. &
                      ieee_is_nan(vvint/tmp) ) then
                  write (*,'(A,I1,A,I1)') "ERROR: poles do not cancel"
                  write (*,*) "results: ", vvint, tmp
                  write (*,*) "ratio = ", vvint/tmp
              else
#endif
                  epinv_checked(currentContrib) = .true.
#if PRINT_EPINV == 1
                  write (*,'(A,E7.1)') "Cancelation of poles checked to precision = ", &
                      abs(vvint/tmp - 1d0)
              endif
!$omp end critical(epinv_print)
#endif
          endif

          if (ieee_is_nan(vvint) .or. (.not. ieee_is_finite(vvint))) then
              write (*,*) "Discarding NaN/Inf matrix element!"
              vvint = 0._dp
              return
          endif

          if (doMultitaucut .and. bin) then
              if (vvint /= 0._dp) then
                  taucut_save = taucut
                  do j=1,size(tcutarray)
                      taucut = tcutarray(j)
                      if (currentContrib == 4) then
                          call singletop_jet_light_heavy_vv_tree(p,msq_tmp)
                          call singletop_jet_light_heavy_vv(p,msqv_tmp)
                          call singletop_jet_light_heavy_vr_z(p,z)
                          scetreweight(j) = vvint_assemble_lxh(msq_tmp,msqv_tmp,z,z1int,z2int,pswt,xx)
                      elseif (currentContrib == 5) then
                          call singletop_light_decay_vv_tree(p,msq_tmp)
                          call singletop_light_decay_vv(p,msqv_tmp)
                          call singletop_light_decay_vr_z(p,z)
                          scetreweight(j) = vvint_assemble_lxd(msq_tmp,msqv_tmp,z,pswt,xx)
                      elseif (currentContrib == 6) then
                          call singletop_heavy_decay_vv_tree(p,msq_tmp)
                          call singletop_heavy_decay_vv(p,msqv_tmp)
                          call singletop_heavy_decay_vr_z(p,z)
                          scetreweight(j) = vvint_assemble_hxd(msq_tmp,msqv_tmp,z,pswt,xx)
                      else
                          error stop "undefined currentContrib in singletop_vvint"
                      endif
                  enddo
                  taucut = taucut_save

                  scetreweight(:) = scetreweight(:) / vvint
              else
                  scetreweight(:) = 1._dp
              endif
          endif

          if (vvint == 0._dp) then
              return
          endif

          if(.not. includedipole(0,p)) then
              vvint = 0._dp
              return
          endif
          call storeptilde(0,p)


          if (doScalevar .and. bin) then
              if (use_DDIS) then
                  truescalevar = (maxscalevar-1)/2
              else
                  truescalevar = maxscalevar
              endif

              if (vvint /= 0._dp) then
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
                          call singletop_jet_light_heavy_vv_tree(p,msq_tmp)
                          call singletop_jet_light_heavy_vv(p,msqv_tmp)
                          call singletop_jet_light_heavy_vr_z(p,z)
                          scalereweight(j) = vvint_assemble_lxh(msq_tmp,msqv_tmp,z,z1int,z2int,pswt,xx)
                      elseif (currentContrib == 5) then
                          call singletop_light_decay_vv_tree(p,msq_tmp)
                          call singletop_light_decay_vv(p,msqv_tmp)
                          call singletop_light_decay_vr_z(p,z)
                          scalereweight(j) = vvint_assemble_lxd(msq_tmp,msqv_tmp,z,pswt,xx)
                      elseif (currentContrib == 6) then
                          call singletop_heavy_decay_vv_tree(p,msq_tmp)
                          call singletop_heavy_decay_vv(p,msqv_tmp)
                          call singletop_heavy_decay_vr_z(p,z)
                          scalereweight(j) = vvint_assemble_hxd(msq_tmp,msqv_tmp,z,pswt,xx)
                      else
                          error stop "undefined currentContrib in singletop_vvint"
                      endif
                  enddo

                  scalereweight(1:maxscalevar) = scalereweight(1:maxscalevar) / vvint
                  call singletop2_scale_setup(p)
              else
                  scalereweight(1:maxscalevar) = 1._dp
              endif
          endif


          if (maxPDFsets > 0 .and. bin) then
              do ipdf=1,maxPDFsets
                  currentPDF = ipdf

                  if (doPDFAlphas) then
                      call singletop2_scale_setup(p)
                  endif

                  if (currentContrib == 4) then
                      call singletop_jet_light_heavy_vv_tree(p,msq_tmp)
                      call singletop_jet_light_heavy_vv(p,msqv_tmp)
                      call singletop_jet_light_heavy_vr_z(p,z)
                      pdfreweight(ipdf) = (vvint - vvint_assemble_lxh(msq,msqv,z,z1int,z2int,pswt,xx))*wgt
                  elseif (currentContrib == 5) then
                      call singletop_light_decay_vv_tree(p,msq_tmp)
                      call singletop_light_decay_vv(p,msqv_tmp)
                      call singletop_light_decay_vr_z(p,z)
                      pdfreweight(ipdf) = (vvint - vvint_assemble_lxd(msq_tmp,msqv_tmp,z,pswt,xx))*wgt
                  elseif (currentContrib == 6) then
                      call singletop_heavy_decay_vv_tree(p,msq_tmp)
                      call singletop_heavy_decay_vv(p,msqv_tmp)
                      call singletop_heavy_decay_vr_z(p,z)
                      pdfreweight(ipdf) = (vvint - vvint_assemble_hxd(msq_tmp,msqv_tmp,z,pswt,xx))*wgt
                  else
                      error stop "undefined currentContrib in singletop_vvint"
                  endif
              enddo

              currentPDF = 0
          endif

          if (bin) then
              call getptildejet(0,pjet)
              currentNd = 0
              includeTaucutgrid(currentNd) = .true.

              call nplotter_new(pjet,vvint*wgt)
              call threadStorageOp(shtmpcommit)
          endif


          if (abs(vvint*wgt) > wtmax) then
              wtmax = abs(vvint*wgt)
          endif

          return
      end function

end module
