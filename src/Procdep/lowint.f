!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function lowint(r,wgt)
        use ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
        use MCFMStorage
        use Scalevar
        use PDFerrors
        use LHAPDF
        use SCET
        use singletop2_m
        use singletop2_scale_m
        use singletop2_realamps_m
        use singletop_lowint, only: singletop_lowintf => lowint
        use VVconfig_m
        use bbfrac_m, only : bbfrac
        use m_gencuts, only : reweight_user, enable_reweight_user
        use Multichannel
        use SafetyCuts, only : passed_smallnew
        use MCFMSetupPlots, only: nplotter_new
        use MCFMSettings
        use nnlo_z1jet, only: genps_z1jet_r
        use ptveto, only: jetptveto
      implicit none
      real(dp):: lowint
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'limits.f'
      include 'npart.f'
      include 'debug.f'
      include 'new_pspace.f'
      include 'vegas_common.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'noglue.f'
      include 'kprocess.f'
      include 'maxwt.f'
      include 'phasemin.f'
      include 'xmin.f'
      include 'wts_bypart.f'
      include 'stopscales.f'
      include 'frag.f'
      include 'ipsgen.f'
      include 'outputoptions.f'
      include 'runstring.f'
      include 'energy.f'
      include 'VVstrong.f'
      include 'dm_params.f'
      include 'initialscales.f'
      include 'toploopgaga.f'
      include 'qcdcouple.f'
      include 'couple.f'
      include 'nlooprun.f'
      include 'badpoint.f'
      include 'ewcorr.f'
      include 'scalevar.f'
      include 'nflav.f'
      include 'x1x2.f'
      include 'bypart.f'
      include 'taucut.f'! for usescet
      include 'nproc.f'
      include 'kpart.f'
      include 'beamtype.f'
      include 'bsm_higgs.f'
      integer:: i,j,k,l,nvec
      integer:: itrial
      real(dp):: alphas,msqtrial,xmsqvar(2),
     & fx1up(-nf:nf),fx2up(-nf:nf),fx1dn(-nf:nf),fx2dn(-nf:nf),
     & r(mxdim),W,xmsq,val,val2,
     & fx1(-nf:nf),fx2(-nf:nf),p(mxpart,4),pjet(mxpart,4),pswt,
     & fx1_H(-nf:nf),fx2_H(-nf:nf),fx1_L(-nf:nf),fx2_L(-nf:nf),
     & fxb1(-nf:nf),fxb2(-nf:nf),
     & wgt,msq(-nf:nf,-nf:nf),xmsqjk,
     & msq_qqb(-nf:nf,-nf:nf), msq_sig_sm(-nf:nf,-nf:nf),msq_bkg_sm(-nf:nf,-nf:nf),msq_int_sm(-nf:nf,-nf:nf),msq_sbi_sm(-nf:nf,-nf:nf),msq_sig_bsm(c6_nval*ct_nval*cg_nval,-nf:nf,-nf:nf),msq_int_bsm(c6_nval*ct_nval*cg_nval,-nf:nf,-nf:nf),msq_sbi_bsm(c6_nval*ct_nval*cg_nval,-nf:nf,-nf:nf),
     & xmsqjk_noew,msq_noew(-nf:nf,-nf:nf),xmsq_noew,lowint_noew,
     & flux,vol,vol_mass,vol3_mass,vol_wt,BrnRat,scaleup,scaledn 
      logical:: bin,includedipole,nopdf
      external qg_tbq,BSYqqb_QQbdk_gvec,qqb_QQbdk,qg_tbqdk,qg_tbqdk_gvec,
     & qqb_Waa,qqb_Waa_mad,qqb_Zbbmas,qqb_totttZ,qqb_totttZ_mad
      common/bin/bin
      common/BrnRat/BrnRat
      include 'bqscale.f'
      external qq_tchan_ztq,qq_tchan_ztq_mad
      external qq_tchan_htq,qq_tchan_htq_mad,qq_tchan_htq_amp
      external qqb_gamgam_g,qqb_gmgmjt_gvec,gg_hzgamg,gg_hg_zgam_gvec

      if (nproc == 1610 .or. nproc == 1650) then
          lowint = singletop_lowintf(r,wgt)
          return
      endif

      lowint=0._dp
c--- ensure isolation code does not think this is fragmentation piece
      z_frag=0._dp

      W=sqrts**2
      p(:,:)=0._dp
      pjet(:,:)=0._dp

      currentNd = 0
      currentPDF=0

      if (doScalevar .and. bin) then
          scalereweight(:) = 1._dp
      endif

      if (maxPDFsets > 0) pdfreweight(:) = 0._dp

      if (usescet .and. kcase == kZ_2jet) then
          if ( genps_z1jet_r(r,p,pswt) .eqv. .false. ) then
              goto 999
          endif
      else
          call gen_lops(r,p,pswt,*999)
      endif

      if (all(.not. ieee_is_nan(p(1:npart+2,:))) .eqv. .false.) then
          if (debug) then
              write(6,*) 'Discarding NaN or infinite phase space point'
          endif
          goto 999
      endif

      nvec=npart+2
      call dotem(nvec,p,s)

c small safety cuts
      if (origkpart == kresummed) then
          if (.not. passed_smallnew(p,npart,1d-9)) then
              goto 999
          endif
      else
          if (.not. passed_smallnew(p,npart,1d-9)) then
              goto 999
          endif
      endif

      xx(1)=-2._dp*p(1,4)/sqrts
      xx(2)=-2._dp*p(2,4)/sqrts
      if ((xx(1) >  1._dp) .or. (xx(2) >  1._dp)
     &.or.(xx(1) < xmin)   .or. (xx(2) < xmin)) then
         goto 999
      endif

  47  continue ! 47 is the label to loop over single top beam contributions
      lowint=0._dp

c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead bail out
      includeTaucutgrid(0) = .true.
      if (includedipole(0,p) .eqv. .false.) then
        goto 888
      endif

  771 continue
      if (dynamicscale) then
          call scaleset(initscale,initfacscale,p)
      else
          call usescales(initscale,initfacscale)
      endif

      if (doPDFAlphas) then
          call updateAlphas(scale)
      endif

c      if (debug) write(*,*) 'Reconstructed x1,x2 ',xx(1),xx(2)

      if ((doScalevar .and. currentPDF == 0 .and. bin) .and. (foundpow(currentPart) .eqv. .false.)) then
        itrial=1
      endif
   66 continue

c--- flux factor
      flux=fbGeV2/(2._dp*xx(1)*xx(2)*W)

      nopdf=.false.

c--- Calculate the required matrix elements
      if     ((kcase==kW_only) .or. (kcase==kWln_ew)) then
        call qqb_w(p,msq)
      elseif (kcase==kW_1jet) then
        call qqb_w_g(p,msq)
c        call qqb_w_gbis(p,msq1)
c        do j=-nf,nf
c        do k=-nf,nf
c        if (msq(j,k)  /=  0._dp) write(6,*) msq(j,k)/msq1(j,k)
c        enddo
c        enddo
c        stop
      elseif (kcase==kWgamma) then
        call qqb_wgam(p,msq)
      elseif (kcase==kWgajet) then
         call qqb_wgamg(p,msq)
      elseif (kcase==kWga2jt) then
         call qqb_waj_g(p,msq)
      elseif ((kcase==kWgaj_a) .or. (kcase==kWga_ew) .or. (kcase==kWln_aq)) then
        write(6,*) 'Not a leading order process'
        stop
      elseif (kcase==kWgajew) then
        call qqb_wgamg(p,msq)
      elseif (kcase==kWbfrmc) then
        call qqb_wbfromc(p,msq)
      elseif (kcase==kW_cjet) then
        call qqb_w_cjet(p,msq)
      elseif (kcase==kWcjet0) then
        call qqb_w_cjet_massless(p,msq)
      elseif (kcase==kWbbmas) then
        call qqb_wbbm(p,msq)
      elseif (kcase==kWbbjem) then
        call qqb_wbbm_g(p,msq)
      elseif (kcase==kWttmas) then
        call qqb_wbbm(p,msq)
      elseif (kcase==kWbbbar) then
        call qqb_wbb(p,msq)
      elseif (kcase==kW_2jet) then
        call qqb_w2jet(p,msq)
      elseif (kcase==kW_3jet) then
        call qqb_w2jet_g(p,msq)
      elseif (kcase==kWbbjet) then
        call qqb_wbb_g(p,msq)
      elseif (kcase==kZ_only) then
        call qqb_z(p,msq)
        if (kewcorr /= knone) then
          msq_noew=msq
          if     (kewcorr == ksudakov) then
            call qqb_z_ew_sudakov(p,msq)
          elseif (kewcorr == kexact) then
            call qqb_z_ew_exact(p,msq)
          endif
        endif
      elseif (kcase==kgg2lep) then
        call ggdilep(p,msq)
      elseif (kcase==kZ_1jet) then
        call qqb_z1jet(p,msq)
c        call qqb_z1jetbis(p,msq1)
c        do j=-nf,nf
c        do k=-nf,nf
c        if (msq(j,k)  /=  0._dp) write(6,*) msq(j,k)/msq1(j,k)
c        enddo
c        enddo
c        pause
      elseif (kcase==kZ_2jet) then
        call qqb_z2jet(p,msq)
      elseif (kcase==kZ_3jet) then
        call qqb_z2jet_g(p,msq)
      elseif (kcase==kZgamma) then
        call set_anomcoup(p)
        ! distinction between decays is inside the subroutine
        call qqb_zgam_new(p,msq)
      elseif (kcase==kZ_2gam) then
        call qqb_zaa(p,msq)
      elseif (kcase==kW_2gam) then
c        if (checkpiDpjk(p)) goto 999
        call qqb_Waa(p,msq)
      elseif (kcase==kZgajet) then
        if (decayChannel() == decayQuarks) then
          call qqb_zaj_vdecay(p,msq)
        else
          call set_anomcoup(p)
          call qqb_zaj(p,msq)
        endif
      elseif (kcase==kZ2gajt) then
        call qqb_zaa_g(p,msq)
      elseif (kcase==kZga2jt) then
        call set_anomcoup(p)
        call qqb_zaj_g(p,msq)
c        call qqb_zaj_g_mad(p,msq1)
c        call qqb_zaj_g(p,msq)
c        do j=-4,4
c        do k=-4,4
c        if (msq(j,k)  /=  zip) write(6,*) j,k,msq(j,k),msq(j,k)/msq1(j,k)
c        enddo
c        enddo
c        stop
      elseif (kcase==kZbbmas) then
        call qqb_zbbm(p,msq)
      elseif (kcase==kZbbbar) then
        call qqb_zbb(p,msq)
      elseif (kcase==kZbbjet) then
        call qqb_zbb_g(p,msq)
      elseif (kcase==kWWqqbr) then
        call qqb_ww(p,msq)
        msq_sig_bsm(:,:,:) = 0._dp
        msq_int_bsm(:,:,:) = 0._dp
        msq_sbi_bsm(:,:,:) = 0._dp
        msq_sig_sm(:,:) = 0._dp
        msq_bkg_sm(:,:) = 0._dp
        msq_int_sm(:,:) = 0._dp
        msq_sbi_sm(:,:) = 0._dp
      elseif (kcase==kWWnpol) then
        call qqb_ww_unpol(p,msq)
      elseif (kcase==kWpWp2j) then
        call qqb_wpwp_qqb(p(1:12,:),msq(-5:5,-5:5))
      elseif (kcase==kWpWp3j) then
        call qqb_wpwp_qqb_g(p(1:12,:),msq(-5:5,-5:5))
      elseif (kcase==kWpmZjj) then
        call qqb_WZjj(p,msq)
      elseif (kcase==kWpmZbj) then
        call qqb_WZbj(p,msq)
      elseif (kcase==kWpmZbb) then
        call qqb_WZbb(p,msq)
      elseif (kcase==kWZbbar) then
        call qqb_wz(p,msq)
      elseif (kcase==kWpmZjj) then
        call qqb_WZjj(p,msq)
      elseif (kcase==kWpmZbj) then
        call qqb_WZbj(p,msq)
      elseif (kcase==kWpmZbb) then
        call qqb_WZbb(p,msq)
      elseif (kcase==kWZbbar) then
        call qqb_wz(p,msq)
      elseif (kcase==kZZlept) then
        call qqb_zz(p,msq)
        msq_sig_bsm(:,:,:) = 0._dp
        msq_int_bsm(:,:,:) = 0._dp
        msq_sbi_bsm(:,:,:) = 0._dp
        msq_sig_sm(:,:) = 0._dp
        msq_bkg_sm(:,:) = 0._dp
        msq_int_sm(:,:) = 0._dp
        msq_sbi_sm(:,:) = 0._dp
      elseif (kcase==kVVlept) then
        call qqb_vv(p,msq)
      elseif (kcase==kWHbbar) then
        call qqb_wh(p,msq)
      elseif (kcase==kWH1jet) then
        call qqb_WH1jet(p,msq)
      elseif (kcase==ktwojet) then
c         reweight=s(1,2)**2*reweight
        call qqb_twojet(p,msq)
        if (kewcorr /= knone) then
          msq_noew=msq
          if     (kewcorr == ksudakov) then
            call qqb_twojet_ew_sudakov(p,msq)
          elseif (kewcorr == kexact) then
            stop 'not implemented yet'
          endif
        endif
      elseif (kcase==ktwo_ew) then
c         reweight=s(1,2)**2*reweight
        call qqb_twojet_ew(p,msq)
        call qqb_twojet(p,msq_noew)
c      elseif (kcase==kthrjet) then
c        call qqb_3jet(p,msq)
      elseif (kcase==kdirgam) then
        call qqb_dirgam(p,msq)
      elseif (kcase==khflgam) then
        call qqb_hflgam(p,msq)
      elseif (kcase==kgamgam) then
        call qqb_gamgam(p,msq)
c        write(6,*) 'gamgam',msq(1,-1),msq(2,-2),msq(-1,1),msq(-2,2)
c        call qqb_gamgam_mad(p,msq)
c        write(6,*) 'gamgam',msq(1,-1),msq(2,-2),msq(-1,1),msq(-2,2)
c        pause
      elseif (kcase == kgg2gamjt) then
         call gg_2gam_g(p,msq)
      elseif (kcase==kgg2gam) then
         if(toploopgaga) then
            msq(:,:)=zip
c            call gg_2gam(p,msq)
            call gggaga_mt(p,msq(0,0))
         else
            call gg_2gam(p,msq)
         endif
      elseif (kcase==kgmgmjt) then
c      call checkgvec(+2, 0,2,p,qqb_gamgam_g,qqb_gmgmjt_gvec)
c      call checkgvec(-1, 0,2,p,qqb_gamgam_g,qqb_gmgmjt_gvec)
c      call checkgvec( 0,-1,1,p,qqb_gamgam_g,qqb_gmgmjt_gvec)
c      call checkgvec( 0, 2,1,p,qqb_gamgam_g,qqb_gmgmjt_gvec)
c      call checkgvec( 1,-1,5,p,qqb_gamgam_g,qqb_gmgmjt_gvec)
c      call checkgvec(-1, 1,5,p,qqb_gamgam_g,qqb_gmgmjt_gvec)
c        call qqb_gamgam_g(p,msq)
         call qqb_gmgmjt(p,msq)
      elseif (kcase==ktrigam) then
        call qqb_trigam(p,msq)
      elseif (kcase==kfourga) then
         call qqb_fourgam(p,msq)
      elseif (kcase==kgamjet) then
        call qqb_dirgam_g(p,msq)
      elseif (kcase==kWH__WW) then
        call qqb_wh_ww(p,msq)
      elseif (kcase==kWH__ZZ) then
        call qqb_wh_zz(p,msq)
      elseif (kcase==kWHgaga) then
        call qqb_wh_gaga(p,msq)
      elseif (kcase==kZHbbar) then
        call qqb_zh(p,msq)
      elseif (kcase==kZH1jet) then
        call qqb_ZH1jet(p,msq)
      elseif (kcase==kZH__WW) then
        call qqb_zh_ww(p,msq)
      elseif (kcase==kZH__ZZ) then
        call qqb_zh_zz(p,msq)
      elseif (kcase==kZHgaga) then
        call qqb_zh_gaga(p,msq)
      elseif (kcase==kggfus0) then
        call gg_h(p,msq)
      elseif (kcase==kggscl0) then
        call gg_h_scalar(p,msq)
      elseif (kcase==kHigaga) then
        call gg_hgamgam(p,msq)
      elseif (kcase==kHi_Zga) then
        call gg_hzgam(p,msq)
      elseif (kcase==kHi_Zaj) then
c        call checkgvec(+2, 0,2,p,gg_hzgamg,gg_hg_zgam_gvec)
c        call checkgvec(-1, 0,2,p,gg_hzgamg,gg_hg_zgam_gvec)
c        call checkgvec( 0,-1,1,p,gg_hzgamg,gg_hg_zgam_gvec)
c        call checkgvec( 0, 2,1,p,gg_hzgamg,gg_hg_zgam_gvec)
c        call checkgvec( 1,-1,6,p,gg_hzgamg,gg_hg_zgam_gvec)
c        call checkgvec(-1, 1,6,p,gg_hzgamg,gg_hg_zgam_gvec)
        call gg_hzgamg(p,msq)
      elseif (kcase==kHWW_4l) then
        call qqb_hww(p,msq)
      elseif (kcase==kHWW2lq) then
        call qqb_hww(p,msq)
      elseif (kcase==kHWW_tb) then
        call qqb_hww_tb(p,msq)
      elseif ((kcase==kHWWint) .or. (kcase==kHWWHpi)
     &   .or. (kcase==kggWW4l)) then
        call gg_ww_int(p,msq)
      elseif (kcase==kggWWbx) then
        msq(:,:)=zip
        call gg_WW(p,msq(0,0))
      elseif (kcase==kHZZ_4l) then
        call qqb_hzz(p,msq)
      elseif (kcase==kHZZ_tb) then
        msq_sig_bsm(:,:,:) = 0._dp
        msq_int_bsm(:,:,:) = 0._dp
        msq_sbi_bsm(:,:,:) = 0._dp
        if (bsm_higgs_scenario .eq. "eft") then 
        cx = cx_sm
        l = 1
        do i = 1, c6_nval
          do j = 1, ct_nval
            do k = 1, cg_nval
              c6 = c6_init + (i-1)*c6_step
              ct = ct_init + (j-1)*ct_step
              cg = cg_init + (k-1)*cg_step
              call gg_hzz_tb(p, msq_sig_bsm(l, :, :))
              call gg_zz_int(p, msq_int_bsm(l, :, :))
              call gg_zz_all(p, msq_sbi_bsm(l, :, :))
              l = l + 1
            enddo
          enddo
        enddo
        endif
        c6 = c6_sm
        ct = ct_sm
        cg = cg_sm
        msq_sig_sm(:,:) = 0._dp
        msq_bkg_sm(:,:) = 0._dp
        msq_int_sm(:,:) = 0._dp
        msq_sbi_sm(:,:) = 0._dp
        call gg_hzz_tb(p,msq_sig_sm)
        call gg_ZZ(p,msq_bkg_sm(0,0))
        call gg_zz_int(p,msq_int_sm)
        call gg_zz_all(p,msq_sbi_sm) 
        msq = msq_sig_sm
      elseif (kcase==kHVV_tb) then
        call gg_hvv_tb(p,msq)
      elseif (kcase==kggVV4l) then
        call gg_VV_all(p,msq)
      elseif (kcase==kggVVbx) then
        msq(:,:)=0._dp
        call gg_VV(p,msq(0,0))
      elseif (kcase==kHZZint) then
        msq_sig_bsm(:,:,:) = 0._dp
        msq_int_bsm(:,:,:) = 0._dp
        msq_sbi_bsm(:,:,:) = 0._dp
        if (bsm_higgs_scenario .eq. "eft") then 
        cx = cx_sm
        l = 1
        do i = 1, c6_nval
          do j = 1, ct_nval
            do k = 1, cg_nval
              c6 = c6_init + (i-1)*c6_step
              ct = ct_init + (j-1)*ct_step
              cg = cg_init + (k-1)*cg_step
              call gg_hzz_tb(p, msq_sig_bsm(l, :, :))
              call gg_zz_int(p, msq_int_bsm(l, :, :))
              call gg_zz_all(p, msq_sbi_bsm(l, :, :))
              l = l + 1
            enddo
          enddo
        enddo
        endif
        c6 = c6_sm
        ct = ct_sm
        cg = cg_sm
        msq_sig_sm(:,:) = 0._dp
        msq_bkg_sm(:,:) = 0._dp
        msq_int_sm(:,:) = 0._dp
        msq_sbi_sm(:,:) = 0._dp
        call gg_hzz_tb(p,msq_sig_sm)
        call gg_ZZ(p,msq_bkg_sm(0,0))
        call gg_zz_int(p,msq_int_sm)
        call gg_zz_all(p,msq_sbi_sm) 
        msq = msq_int_sm
      elseif (kcase==kHZZHpi) then
        cx = cx_sm
        c6 = c6_sm
        ct = ct_sm
        cg = cg_sm
        call gg_zz_Hpi(p,msq)
      elseif (kcase==kggZZ4l) then
        msq_sig_bsm(:,:,:) = 0._dp
        msq_int_bsm(:,:,:) = 0._dp
        msq_sbi_bsm(:,:,:) = 0._dp
        if (bsm_higgs_scenario .eq. "eft") then 
        cx = cx_sm
        l = 1
        do i = 1, c6_nval
          do j = 1, ct_nval
            do k = 1, cg_nval
              c6 = c6_init + (i-1)*c6_step
              ct = ct_init + (j-1)*ct_step
              cg = cg_init + (k-1)*cg_step
              call gg_hzz_tb(p, msq_sig_bsm(l, :, :))
              call gg_zz_int(p, msq_int_bsm(l, :, :))
              call gg_zz_all(p, msq_sbi_bsm(l, :, :))
              l = l + 1
            enddo
          enddo
        enddo
        endif
        c6 = c6_sm
        ct = ct_sm
        cg = cg_sm
        msq_sig_sm(:,:) = 0._dp
        msq_bkg_sm(:,:) = 0._dp
        msq_int_sm(:,:) = 0._dp
        msq_sbi_sm(:,:) = 0._dp
        call gg_hzz_tb(p,msq_sig_sm)
        call gg_ZZ(p,msq_bkg_sm(0,0))
        call gg_zz_int(p,msq_int_sm)
        call gg_zz_all(p,msq_sbi_sm) 
        msq = msq_sbi_sm
      elseif (kcase==kggZZbx) then
        msq_sig_bsm(:,:,:) = 0._dp
        msq_int_bsm(:,:,:) = 0._dp
        msq_sbi_bsm(:,:,:) = 0._dp
        if (bsm_higgs_scenario .eq. "eft") then 
        cx = cx_sm
        l = 1
        do i = 1, c6_nval
          do j = 1, ct_nval
            do k = 1, cg_nval
              c6 = c6_init + (i-1)*c6_step
              ct = ct_init + (j-1)*ct_step
              cg = cg_init + (k-1)*cg_step
              call gg_hzz_tb(p, msq_sig_bsm(l, :, :))
              call gg_zz_int(p, msq_int_bsm(l, :, :))
              call gg_zz_all(p, msq_sbi_bsm(l, :, :))
              l = l + 1
            enddo
          enddo
        enddo
        endif
        c6 = c6_sm
        ct = ct_sm
        cg = cg_sm
        msq_sig_sm(:,:) = 0._dp
        msq_bkg_sm(:,:) = 0._dp
        msq_int_sm(:,:) = 0._dp
        msq_sbi_sm(:,:) = 0._dp
        call gg_hzz_tb(p,msq_sig_sm)
        call gg_ZZ(p,msq_bkg_sm(0,0))
        call gg_zz_int(p,msq_int_sm)
        call gg_zz_all(p,msq_sbi_sm) 
        msq = msq_bkg_sm
      elseif (kcase==kppZZ4l) then
        msq_sig_bsm(:,:,:) = 0._dp
        msq_int_bsm(:,:,:) = 0._dp
        msq_sbi_bsm(:,:,:) = 0._dp
        msq_sig_sm(:,:) = 0._dp
        msq_bkg_sm(:,:) = 0._dp
        msq_int_sm(:,:) = 0._dp
        msq_sbi_sm(:,:) = 0._dp
        msq_qqb(:,:) = 0._dp
        cx = cx_sm
        c6 = c6_sm
        ct = ct_sm
        cg = cg_sm
        call gg_zz_all(p,msq_sbi_sm) 
        call qqb_zz(p,msq_qqb)
        msq = msq_sbi_sm
        msq = msq+msq_qqb 
      elseif (kcase==kHZZqgI) then
         call qg_Hint_ZZ(p,msq)
      elseif (kcase==kH_1jet) then
        call qqb_hg(p,msq)
      elseif (kcase==kttZbbl) then
        call qqbZtt(p,msq)
      elseif ((kcase==ktt_bbl) .or. (kcase==ktt_bbh)) then
        call qqb_QQbdk(p,msq)
      elseif (kcase==ktt_bbu) then
        call qqb_QQbdku(p,msq)
      elseif (kcase==kqq_ttg) then
       call qqb_QQbdk_g(p,msq)
      elseif (kcase==ktt_tot) then
        call qqb_QQb(p,msq)
        if (kewcorr /= knone) then
          msq_noew=msq
          if     (kewcorr == ksudakov) then
            call qqb_QQb_ew_sudakov(p,msq)
          elseif (kewcorr == kexact) then
            stop 'not implemented yet'
          endif
        endif
      elseif (kcase==kbb_tot) then
        call qqb_QQb(p,msq)
      elseif (kcase==kcc_tot) then
        call qqb_QQb(p,msq)
      elseif (kcase==ktt_glu) then
         call qqb_QQb_g(p,msq)
      elseif (kcase==ktopanom) then
        bbfrac = 0._dp
        call singletop2_tree(p,msq)
      elseif (kcase==kbq_tpq) then
        call bq_tpq(p,msq)
      elseif (kcase==kttdkay) then
        write(6,*) 'This process is not a leading order contribution'
        stop
      elseif (kcase==kt_bbar) then
      call qqb_tbbdk(p,msq)
      elseif (kcase==ktdecay) then
        write(6,*) 'This process is not a leading order contribution'
        stop
      elseif (kcase==kW_tndk) then
        call qqb_w_tndk(p,msq)
      elseif (kcase==kW_twdk) then
        call qqb_w_twdk(p,msq)
      elseif (kcase==kWtbwdk) then
        call qqb_wtbwdk(p,msq)
      elseif (kcase==kWtbndk) then
        call qqb_wtbndk(p,msq)
      elseif (kcase==ktottth) then
        call qqb_tottth(p,msq)
      elseif (kcase==kqq_tth) then
        call qqb_tth(p,msq)
      elseif (kcase==ktth_ww) then
        call qqb_tth(p,msq)
      elseif (kcase==kqq_ttz) then
        call qqb_ttz(p,msq)
      elseif (kcase==kqqtthz) then
        call qqb_ttz(p,msq)
      elseif (kcase==kqq_ttw) then
        call qqb_ttw(p,msq)
      elseif (kcase==khttjet) then
        call qqb_higgs(p,msq)
      elseif (kcase==khttscl) then
        call qqb_higgs_scalar(p,msq)
      elseif (kcase==kggfus1) then
        call gg_hg(p,msq)
      elseif (kcase==khjetma) then
        call qqb_higgs(p,msq)
      elseif (kcase==kHgagaj) then
        call gg_hgagag(p,msq)
      elseif (kcase==kHWWjet) then
        call gg_hWWg(p,msq)
      elseif (kcase==kHZZjet) then
        call gg_hZZg(p,msq)
      elseif (kcase==kHWW2jt) then
        call gg_hWWgg(p,msq)
      elseif (kcase==kHZZ2jt) then
        call gg_hZZgg(p,msq)
      elseif (kcase==kHWW3jt) then
        call gg_hWWggg(p,msq)
      elseif (kcase==kHZZ3jt) then
        call gg_hZZggg(p,msq)
      elseif (kcase==kattjet) then
        call qqb_higgs_odd(p,msq)
      elseif (kcase==kqq_Hqq) then
        call VV_hqq(p,msq)
      elseif (kcase==kqq_Hgg) then
        call VV_Hgaga(p,msq)
      elseif (kcase==kqqHqqg) then
        call VV_hqq_g(p,msq)
      elseif (kcase==kqq_HWW) then
        call VV_HWW(p,msq)
      elseif (kcase==kqq_HZZ) then
        call VV_HZZ(p,msq)
      elseif (kcase==ktautau) then
        call qqb_tautau(p,msq)
      elseif (kcase==kqg_tbq) then
        call qg_tbq(p,msq)
c--- Check of gvec routines
c      call checkgvec(+2,0,2,p,qg_tbq,qg_tbq_gvec)
c      call checkgvec(-1,0,2,p,qg_tbq,qg_tbq_gvec)
      elseif (kcase==kqqZZqq) then
        if (VVstrong) then
          call qq_ZZqqstrong(p,msq)
        else
          call qq_ZZqq(p,msq)
        endif
      elseif (kcase==kqqWWqq) then
        if (VVstrong) then
          call qq_WWqqstrong(p,msq)
        else
          call qq_WWqq(p,msq)
        endif
      elseif (kcase==kqqVVqq) then
        if (VVstrong) then
          write(6,*) 'Abort: kqqVVqq not yet implemented'
          stop
        else
          call qq_VVqq(p,msq)
        endif
      elseif (kcase==kqqWWss) then
        if (VVstrong) then
          call qq_WWssstrong(p,msq)
        else
          call qq_WWss(p,msq)
        endif
      elseif (kcase==kqqWZqq) then
        if (VVstrong) then
          call qq_WZqqstrong(p,msq)
        else
          call qq_WZqq(p,msq)
        endif
      elseif (kcase==kqgtbqq) then
        call qg_tbq_g(p,msq)
      elseif (kcase==k4ftwdk) then
        call qg_tbqdk(p,msq)
      elseif (kcase==k4ftjet) then
        call qg_tbqdk_g(p,msq)
      elseif (kcase==kqq_tbg) then
        call qq_tbg(p,msq)
c--- Check of gvec routines
c      call checkgvec(2,-1,5,p,qq_tbg,qq_tbg_gvec)
c      call checkgvec(-1,2,5,p,qq_tbg,qq_tbg_gvec)
      elseif (kcase==kqqtbgg) then
        call qq_tbg_g(p,msq)
      elseif (kcase==kepem3j) then
        call epem3j(p,msq)
      elseif (kcase==kgQ__ZQ) then
        call gQ_zQ(p,msq)
      elseif (kcase==kZccmas) then
        call qqb_zccm(p,msq)
      elseif (kcase==kggfus2) then
        call gg_hgg(p,msq)
      elseif (kcase==kgagajj) then
        call gg_hgg(p,msq)
      elseif (kcase==kh2jmas) then
        call gg_hgg_mass(p,msq)
        ! to also include the effect of bottom-quark loops
        ! call gg_hgg_mass_tb(p,msq)
      elseif (kcase==kh2jscl) then
        write(6,*) 'Process not fully implemented'
        stop
c        call gg_hgg_scalar(p,msq)
      elseif (kcase==kggfus3) then
        call gg_hggg(p,msq)
      elseif (kcase==kW_bjet) then
        call qqb_wbjet(p,msq)
      elseif (kcase==kWcjetg) then
        call qqb_w_cjet_massless_g(p,msq)
      elseif (kcase==kZ_bjet) then
        call qqb_zbjet(p,msq)
      elseif (kcase==kZbjetg) then
        call qqb_zbjet_g(p,msq)
      elseif (kcase==kH_tjet) then
        call qq_tchan_htq(p,msq)
      elseif (kcase==kH_tdkj) then
        call qq_tchan_htq_dk(p,msq)
      elseif (kcase==kZ_tjet) then
        call qq_tchan_ztq(p,msq)
      elseif (kcase==kZ_tdkj) then
         call qq_tchan_ztq_dk(p,msq)
      elseif (kcase==kZtdk2j) then
         call qq_tchan_ztqg_dk(p,msq)
      elseif (kcase==kZt2jet) then
        call qq_tchan_ztqg(p,msq)
      elseif (kcase==kHHpair) then
        call gg_HH(p,msq)
c        call gg_HH_phase(p,msq1)
c        write(6,*) 'msq(0,0),msq1(0,0),msq1(0,0)/msq(0,0)',msq(0,0),msq1(0,0),msq1(0,0)/msq(0,0)
      elseif (kcase==kdm_jet) then
         call qqb_dm_monojet(p,msq)
      elseif (kcase==kdm_gam) then
         call qqb_dm_monophot(p,msq)
      elseif ( kcase==kdm2jet) then
         call qqb_dm_monojet_g(p,msq)
      elseif ( kcase==kdm_gaj) then
         call qqb_dm_monophot_g(p,msq)
      elseif (kcase==kWW_jet) then
        call qqb_wwg(p,msq)
      elseif (kcase==kWW2jet) then
        call qqb_wwg_g(p,msq)
      elseif (kcase==kWZ_jet) then
        call qqb_wzg(p,msq)
      elseif (kcase==kZZ_jet) then
        call qqb_zzg(p,msq)
      elseif (kcase==kvlchk2) then
        call qqb_vol(p,msq)
        flux=one/vol(W,2)
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=1._dp
        fx2(-1)=1._dp
        nopdf=.true.
      elseif (kcase==kvlchk3) then
        call qqb_vol(p,msq)
        flux=one/vol(W,3)
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=2._dp
        fx2(-1)=2._dp
        nopdf=.true.
      elseif (kcase==kvlchk4) then
        taumin=0.0001_dp
        bbsqmax=W
        bbsqmin=0._dp
        call qqb_vol(p,msq)
        flux=one/vol(W,4)
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=4._dp*xx(1)
        fx2(-1)=2._dp/xx(2)
        nopdf=.true.
      elseif (kcase==kvlchk5) then
        taumin=0.0001_dp
        bbsqmax=W
        bbsqmin=0._dp
        call qqb_vol(p,msq)
        flux=one/vol(W,5)
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=4._dp
        fx2(-1)=4._dp
        nopdf=.true.
      elseif (kcase==kvlchk6) then
        call qqb_vol(p,msq)
        flux=one/vol(W,6)
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=6._dp*xx(1)
        fx2(-1)=4._dp/xx(2)
        nopdf=.true.
      elseif (kcase==kvlchk8) then
        call qqb_vol(p,msq)
        flux=one/vol(W,8)
        bbsqmax=W
        bbsqmin=0._dp
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=6._dp/xx(1)
        fx2(-1)=6._dp/xx(2)
        nopdf=.true.
      elseif (kcase==kvlchkm) then
        taumin=0.0001_dp
        bbsqmax=W
        bbsqmin=0._dp
        call qqb_vol(p,msq)
        flux=one/vol_mass(mb,W)
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=4._dp*xx(1)
        fx2(-1)=2._dp/xx(2)
        nopdf=.true.
      elseif (kcase==kvlchm3) then
        write(6,*) 'Process not fully implemented!'
        stop
        taumin=(2._dp*mt/sqrts)**2
        call qqb_vol(p,msq)
        flux=one/vol3_mass(mt,W)
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=2._dp
        fx2(-1)=2._dp
        nopdf=.true.
      elseif ((kcase==kvlchwt) .or. (kcase==kvlchwn)
     &   .or. (kcase==kvlchwg) .or. (kcase==kvlchwh)) then
        taumin=0.0001_dp
        call qqb_vol(p,msq)
        flux=one/vol_wt(W)
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=1._dp
        fx2(-1)=1._dp
        nopdf=.true.
      else
        write(6,*) 'Unimplemented process in lowint : kcase=',kcase
        stop
      endif
c code to find power of alpha-s for scale variation
      if ((doScalevar .and. currentPDF == 0 .and. bin) .and. (foundpow(currentPart) .eqv. .false.)) then
        if (itrial == 1) then
          msqtrial=maxval(msq)
          if (msqtrial == 0) goto 999
          as=as*two
          ason2pi=ason2pi*two
          ason4pi=ason4pi*two
          gsq=gsq*two
          itrial=itrial+1
          goto 66
        endif
        msqtrial=maxval(msq)/msqtrial
        alphaspow(currentPart) = -1
        if (abs(msqtrial-one) < 1.e-8) alphaspow(currentPart) = 0
        if (abs(msqtrial-two) < 1.e-8) alphaspow(currentPart) = 1
        if (abs(msqtrial-four) < 1.e-8) alphaspow(currentPart) = 2
        if (abs(msqtrial-eight) < 1.e-8) alphaspow(currentPart) = 3
        if (abs(msqtrial-16._dp) < 1.e-8) alphaspow(currentPart) = 4
        if (abs(msqtrial-32._dp) < 1.e-8) alphaspow(currentPart) = 5
        if (alphaspow(currentPart) == -1) then
          write(6,*) 'Unable to determine power of alpha-s for scale variation'
          stop
        endif
        as=as/two
        ason2pi=ason2pi/two
        ason4pi=ason4pi/two
        gsq=gsq/two
        foundpow(currentPart) = .true.
c        write(6,*) 'Found alpha-s power: ',alphaspow
        goto 66
      endif

c--- initialize a PDF set here, if calculating errors
  777 continue
      xmsq=0._dp
      xmsq_noew=0._dp

c--- calculate PDF's
      if (((kcase==kqg_tbq) .or. (kcase==k4ftwdk)) .and. (dynamicscale .eqv. .false.)) then
c--- for single top + b, make sure to use two different scales
         if (doPDFerrors .and. doScalevar) then
             error stop "simultaneous PDF and scale variation not supported for this process"
         endif

         call fdist(ih1,xx(1),facscale_H,fx1_H, 1)
         call fdist(ih2,xx(2),facscale_H,fx2_H, 2)
         call fdist(ih1,xx(1),facscale_L,fx1_L, 1)
         call fdist(ih2,xx(2),facscale_L,fx2_L, 2)

         do j=-nf,nf
           if (j == 0) then   ! heavy quark line has gluon init. state
             fx1(j)=fx1_H(j)
             fx2(j)=fx2_H(j)
           else
             fx1(j)=fx1_L(j)
             fx2(j)=fx2_L(j)
           endif
         enddo
      elseif ((kcase == kbq_tpq .or. kcase == kqg_tbq .or. kcase==ktopanom) .and. use_DDIS) then
          if (doPDFerrors .and. doScalevar) then
              error stop "simultaneous PDF and scale variation not supported for DDIS scales"
          endif

          call singletop2_scale_setup(p)
          ! onheavy or onlight does not matter for born kinematics
          call fdist(ih1,xx(1),facscale_beam1_isheavy_onheavy,fxb1,1)
          call fdist(ih2,xx(2),facscale_beam2_islight_onheavy,fx2,2)
          call fdist(ih1,xx(1),facscale_beam1_islight_onheavy,fx1,1)
          call fdist(ih2,xx(2),facscale_beam2_isheavy_onheavy,fxb2,2)
      else
c--- usual case
            if ((maxPDFsets > 0) .and. bin) then
              call fdist(ih1,xx(1),facscale,fx1,1)
              call fdist(ih2,xx(2),facscale,fx2,2)
              ! this covers the case of PDF error AND scale variation for central PDF
              if (doScalevar .and. currentPDF == 0 .and. bin) then
                  if (vetoscalevar) then
                    scaleup=max(facscale,jetptveto)*two
                    scaledn=min(facscale,jetptveto)/two
                  else
                    scaleup=facscale*two
                    scaledn=facscale/two
                  endif
                  call fdist(ih1,xx(1),scaleup,fx1up,1)
                  call fdist(ih2,xx(2),scaleup,fx2up,2)
                  call fdist(ih1,xx(1),scaledn,fx1dn,1)
                  call fdist(ih2,xx(2),scaledn,fx2dn,2)
                  xmsqvar(:)=zip
              endif
            else
              if (nopdf .eqv. .false.) then
                 call fdist(ih1,xx(1),facscale,fx1,1)
                 call fdist(ih2,xx(2),facscale,fx2,2)
                 if (doScalevar .and. currentPDF == 0 .and. bin) then
                     if (vetoscalevar) then
                       scaleup=max(facscale,jetptveto)*two
                       scaledn=min(facscale,jetptveto)/two
                     else
                       scaleup=facscale*two
                       scaledn=facscale/two
                     endif
                     call fdist(ih1,xx(1),scaleup,fx1up,1)
                     call fdist(ih2,xx(2),scaleup,fx2up,2)
                     call fdist(ih1,xx(1),scaledn,fx1dn,1)
                     call fdist(ih2,xx(2),scaledn,fx2dn,2)
                     xmsqvar(:)=zip
                 endif
              endif
            endif
      endif

      do j=-nflav,nflav
      do k=-nflav,nflav

      if ((selectpdfs(1,j) .eqv. .false.) .or. (selectpdfs(2,k) .eqv. .false.)) cycle

      if     ((kcase==kbq_tpq .or. kcase==ktopanom) .and. use_DDIS) then
c--- special case for dynamic scale in t-channel single top
        if     (abs(j) == 5) then
          xmsqjk=fxb1(j)*fx2(k)*msq(j,k)
      elseif (abs(k) == 5) then
          xmsqjk=fx1(j)*fxb2(k)*msq(j,k)
      else
        xmsqjk=0._dp
      endif
      elseif ((kcase==kqg_tbq) .and. use_DDIS) then
c--- special case for dynamic scale in t-channel single top
        if     (j == 0) then
          xmsqjk=fxb1(j)*fx2(k)*msq(j,k)
      elseif (k == 0) then
          xmsqjk=fx1(j)*fxb2(k)*msq(j,k)
      else
        xmsqjk=0._dp
      endif
      else
c--- DEFAULT
        xmsqjk=fx1(j)*fx2(k)*msq(j,k)
        if (doScalevar .and. currentPDF == 0 .and. bin) then
          xmsqvar(1)=xmsqvar(1)+fx1up(j)*fx2up(k)*msq(j,k)
          xmsqvar(2)=xmsqvar(2)+fx1dn(j)*fx2dn(k)*msq(j,k)
        endif
        if (kewcorr /= knone) xmsqjk_noew=fx1(j)*fx2(k)*msq_noew(j,k)
      endif

      xmsq=xmsq+xmsqjk
      if (kewcorr /= knone) xmsq_noew=xmsq_noew+xmsqjk_noew

      enddo
      enddo ! end loop over partons

      if (currentPDF == 0) then
        lowint=flux*pswt*xmsq/BrnRat
        if (kewcorr /= knone) lowint_noew=flux*pswt*xmsq_noew/BrnRat
      endif

c compute weights for scale variation
      if (doScalevar .and. currentPDF == 0 .and. bin) then
        if (abs(xmsq) > zip) then
          if (vetoscalevar) then
            scalereweight(1)=(alphas(max(scale,jetptveto)*two,amz,nlooprun)/as)**alphaspow(currentPart)
            scalereweight(2)=(alphas(min(scale,jetptveto)/two,amz,nlooprun)/as)**alphaspow(currentPart)
          else
            scalereweight(1)=(alphas(scale*two,amz,nlooprun)/as)**alphaspow(currentPart)
            scalereweight(2)=(alphas(scale/two,amz,nlooprun)/as)**alphaspow(currentPart)
          endif
          scalereweight(1)=scalereweight(1)*xmsqvar(1)/xmsq
          scalereweight(2)=scalereweight(2)*xmsqvar(2)/xmsq
          if (maxscalevar >= 6) then
            scalereweight(3)=scalereweight(1)*xmsq/xmsqvar(1)
            scalereweight(4)=scalereweight(2)*xmsq/xmsqvar(2)
            scalereweight(5)=xmsqvar(1)/xmsq
            scalereweight(6)=xmsqvar(2)/xmsq
          endif
          if (maxscalevar == 8) then
            scalereweight(7)=scalereweight(3)*xmsqvar(2)/xmsq
            scalereweight(8)=scalereweight(4)*xmsqvar(1)/xmsq
          endif
        else
          scalereweight(:)=1._dp
        endif
      endif

c--- loop over all PDF error sets, if necessary
      if ((maxPDFsets > 0) .and. bin) then
        if (currentPDF > 0) then
            pdfreweight(currentPDF) = (lowint - flux*pswt*xmsq/BrnRat)*wgt
        endif

        currentPDF=currentPDF+1
        if (currentPDF <= maxPDFsets) then
            if (doPDFAlphas) then
                ! only if as(mz) has changed recompute matrix element
                if ( abs(getalphas(zmass) - amz) > 1d-5 ) then
                    goto 771
                endif
            endif
            goto 777
        endif
      endif

      call getptildejet(0,pjet)

      call dotem(nvec,pjet,s)

      val=lowint*wgt
      val2=val**2
      if (ieee_is_nan(val)) then
         write(6,*) 'lowint val = ',val
         write(6,*) 'Discarding point with random variables',r
         lowint=zip
         val=zip
         goto 999
      endif

      if ((abs(val) > wtmax)) then
        wtmax=abs(val)
      endif

      if (bin) then
c--- for EW corrections, make additional weight available inside common block
        if (kewcorr /= knone) then
          wt_noew=lowint_noew*wgt
        endif
        if (newStyleHistograms) then
            call nplotter_new(pjet,val)
        else
            call nplotter(pjet,val,val2,0)
        endif
      endif

      if (includeTaucutgrid(0) .eqv. .false.) then
          lowint = 0._dp
      endif

      if (enable_reweight_user) then
          lowint = lowint * reweight_user(pjet)
      endif

      if (bin .and. \
        ( \
        kcase==kWWqqbr .or. \
        kcase==kZZlept .or. \
        kcase==kggZZ4l .or. \
        kcase==kHZZ_tb .or. \
        kcase==kggZZbx .or. \
        kcase==kHZZint .or. \
        kcase==kppZZ4l \
        )) then
        call output_csv(pjet,msq_sig_sm,msq_bkg_sm,msq_int_sm,msq_sbi_sm,msq_sig_bsm,msq_int_bsm,msq_sbi_bsm,val)
      endif

  888 continue ! just cuts did not pass

      return

 999  continue ! everything is zero / invalid
      lowint=0._dp

      return
      end


