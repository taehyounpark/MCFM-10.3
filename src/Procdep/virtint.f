!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function virtint(r,wgt)
        use PDFerrors
        use LHAPDF
        use Scalevar
        use ieee_arithmetic
        use singletop2_m, only : singletop2_z, singletop2_tree, singletop2_virt
        use singletop2_scale_m
        use singletop_virtint, only: singletop_virtintf => virtint
        use bbfrac_m, only : bbfrac
        use MCFMStorage
        use SCET
        use m_gencuts, only : enable_reweight_user, reweight_user
        use Multichannel
        use SafetyCuts, only : passed_smallnew
        use MCFMSetupPlots, only: nplotter_new
        use MCFMSettings
        use ggHwilson, only: expansionorder, Wilsonorder
        use ptveto, only: jetptveto
      implicit none
      real(dp):: virtint
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'noglue.f'
      include 'vegas_common.f'
      include 'sprods_com.f'
      include 'npart.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'agq.f'
      include 'PR_new.f'
      include 'PR_cs_new.f'
      include 'PR_h2j.f'
      include 'PR_twojet.f'
      include 'PR_stop.f'
      include 'PR_mix.f'
      include 'msq_cs.f'
      include 'msq_struc.f'
      include 'msq_mix.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'lc.f'
      include 'kprocess.f'
      include 'maxwt.f'
      include 'limits.f'
      include 'heavyflav.f'
      include 'nflav.f'
      include 'b0.f'
      include 'masses.f'
      include 'wts_bypart.f'
      include 'nores.f'
      include 'stopscales.f'
      include 'stopbmass.f'
      include 'ewcouple.f'
      include 'flags.f'
      include 'frag.f'
      include 'ipsgen.f'
      include 'phasemin.f'
      include 'xmin.f'
      include 'outputoptions.f'
      include 'dm_params.f'
      include 'runstring.f'
      include 'x1x2.f'
      include 'bypart.f'
      include 'energy.f'
      include 'first.f'
      include 'initialscales.f'
      include 'taucut.f'
      include 'ppmax.f'
      include 'kpart.f'
      include 'hbbparams.f'
      include 'ewcorr.f'
      include 'mpicommon.f'
      include 'scalevar.f'
      include 'couple.f'
      include 'nlooprun.f'
      include 'badpoint.f'
      include 'bqscale.f'
      include 'zcouple_cms.f'
      include 'debug.f'
      include 'cutoff.f'
      include 'nproc.f'
      include 'beamtype.f'
      include 'anomcoup.f'
      include 'lib/TensorReduction/Include/TRbadpoint.f' 
      include 'src/Inc/toploops.f'
#ifdef HAVE_RECOLA
      include 'nwz.f'
#endif
      real(dp), intent(in) :: wgt, r(mxdim)

      real(dp):: mqq(0:2,-nf:nf,-nf:nf),
     & msqx(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      real(dp):: ppmsqx(0:2,ppmax)
      real(dp):: msqx_cs(0:2,-nf:nf,-nf:nf)
      real(dp):: AP(-1:1,-1:1,3),APqg_mass,AP_mix(-1:1,-1:1,0:3,3),APqed(-1:1,-1:1,3)

      integer:: j,k,m,n,cs,ics,csmax,nvec,is,iq,ia
      integer itrial
      real(dp):: xmsqvar(8)
      real(dp):: p(mxpart,4),pjet(mxpart,4),W,xmsq,pttwo,
     & val,val2,fx1(-nf:nf),fx2(-nf:nf),fx1z(-nf:nf),fx2z(-nf:nf),xmsqt,
     & fx1_H(-nf:nf),fx2_H(-nf:nf),fx1_L(-nf:nf),fx2_L(-nf:nf),
     & fx1z_H(-nf:nf),fx2z_H(-nf:nf),fx1z_L(-nf:nf),fx2z_L(-nf:nf)
      real(dp):: xmsq_noew,msq_noew(-nf:nf,-nf:nf),virtint_noew
      real(dp):: pswt,xjac,msqqcd(-nf:nf,-nf:nf),
     & msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),msqvdk(-nf:nf,-nf:nf),
     & msqvdkW(-nf:nf,-nf:nf),
     & msq_qq,msq_aa,msq_aq,msq_qa,msq_qg,msq_gq,epcorr,msq_up,msq_dn,
     & msq1(-nf:nf,-nf:nf),msq0(-nf:nf,-nf:nf),msqm1(-nf:nf,-nf:nf),xlog,
     & bit2(-nf:nf,-nf:nf),bit1(-nf:nf,-nf:nf),bit0(-nf:nf,-nf:nf)
      real(dp) :: msqv_light(-nf:nf,-nf:nf), msqv_heavy(-nf:nf,-nf:nf) ! for ktopanom
      real(dp):: z,x1onz,x2onz,flux,omz,ptthree,
     & BrnRat,xmsq_old,tmp,wwidth_save,zwidth_save
      real(dp):: scaleup,scaledn,facscaleup,facscaledn
      integer:: rvcolourchoice
      logical:: bin,includedipole,WWjetcheckpiDpjk
      real(dp):: QandGint

      common/bin/bin
      common/BrnRat/BrnRat
      common/rvcolourchoice/rvcolourchoice
c      common/ggZZunstable/ggZZunstable
c      data p/56*0._dp/
      integer, save:: nshot=1
      external gg_ZZ,qqb_w1jet_vbis
!$omp threadprivate(/rvcolourchoice/)
!$omp threadprivate(nshot)

      real(dp) fxa
      common/photonpdf/fxa
!$omp threadprivate(/photonpdf/)

      if (nproc == 1610 .or. nproc == 1650) then
          virtint = singletop_virtintf(r,wgt)
          return
      endif

      QandGflag=.false.
      if (first) then
         first=.false.
c         write(*,*) case
         nshot=1
      endif

      virtint=0._dp
c--- ensure isolation code does not think this is fragmentation piece
      z_frag=0._dp

      p(:,:)=0._dp
      pjet(:,:)=0._dp

      if (doScalevar .and. bin) then
          scalereweight(:) = 1._dp
      endif

      W=sqrts**2

      currentPDF = 0
      if (maxPDFsets > 0) pdfreweight = 0._dp

      call gen_lops(r,p,pswt,*999)

      if (all(.not. ieee_is_nan(p(1:npart+2,:))) .eqv. .false.) then
          if (debug) then
              write(6,*) 'Discarding NaN or infinite phase space point'
          endif
          goto 999
      endif

      nvec=npart+2
      call dotem(nvec,p,s)

c--- (moved to includedipole) impose cuts on final state
c      call masscuts(p,*999)

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

c--- back-to-back check
      if ((kcase == kWWqqbr) .or. (kcase == kWZbbar) .or. (kcase == kZZlept)) then
        if (pttwo(3,4,p) < 1.e-3_dp) goto 999
      endif
c  47  continue ! 47 is the label to loop over single top beam contributions
      virtint=0._dp

c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead bail out
      includeTaucutgrid(0) = .true.
      if (includedipole(0,p) .eqv. .false.) then
        goto 999
      endif

c--- test to see whether we need Gflag and Qflag together
      if ( ((kcase==kW_2jet) .or. (kcase==kZ_2jet))
     &.and. (Qflag) .and. (Gflag) ) then
        QandGflag=.true.
        QandGint=0._dp
c--- first pass: Gflag
        Gflag=.true.
        Qflag=.false.
      endif

c--- restart from here when calculating with Qflag and Gflag
c--- (W+2 jet and Z+2 jet processes only)
   44 continue

      if (doScalevar .and. currentPDF == 0 .and. bin) then
        itrial=maxscalevar+1
      endif

   66 continue

  771 continue
      if (dynamicscale) then
          call scaleset(initscale,initfacscale,p)
      else
          call usescales(initscale,initfacscale)
      endif

      if (doPDFAlphas) then
          call updateAlphas(scale)
      endif

c adjust scales for scale variation
      if (doScalevar .and. currentPDF == 0 .and. bin) then
          if (vetoScalevar) then
            scaleup=max(scale,jetptveto)*2._dp
            scaledn=min(scale,jetptveto)/2._dp
            facscaleup=max(facscale,jetptveto)*2._dp
            facscaledn=min(facscale,jetptveto)/2._dp
          else
            scaleup=scale*2._dp
            scaledn=scale/2._dp
            facscaleup=facscale*2._dp
            facscaledn=facscale/2._dp
          endif
          if (itrial == 9) then
             call usescales(scaledn, facscaleup)
          elseif (itrial == 8) then
             call usescales(scaleup, facscaledn)
          elseif (itrial == 7) then
             call usescales(scale, facscaledn)
          elseif (itrial == 6) then
              call usescales(scale, facscaleup)
          elseif (itrial == 5) then
              call usescales(scaledn, facscale)
          elseif (itrial == 4) then
              call usescales(scaleup, facscale)
          elseif (itrial == 3) then
              call usescales(scaledn, facscaledn)
          elseif (itrial == 2) then
              call usescales(scaleup, facscaleup)
          endif
      endif

      xx(1)=-2._dp*p(1,4)/sqrts
      xx(2)=-2._dp*p(2,4)/sqrts
      if ((xx(1) >  1._dp) .or. (xx(2) >  1._dp)
     &.or.(xx(1) < xmin)   .or. (xx(2) < xmin)) then
         goto 999
      endif

      z=r(ndim)**2
! Exclude exceptional points
      if ((z == 0._dp) .or. (z == 1._dp)) goto 999
c      if (nshot == 1) z=0.95_dp
      xjac=two*sqrt(z)

      omz=1._dp-z

      flux=fbGeV2/(2._dp*xx(1)*xx(2)*W)
c--- for mlm study, divide by (Ecm)**2=W
c      if (runstring(1:3) == 'mlm') then
c      flux=flux/W
c      endif

c--- to test poles, we need colourchoice=0, but save real value
      if (nshot == 1) then
        rvcolourchoice=colourchoice
        colourchoice=0
      endif

   12 continue
c--- point to restart from when checking epsilon poles

c--- correction to epinv from AP subtraction when mu_FAC != mu_REN,
c--- corresponding to subtracting -1/epinv*Pab*log(musq_REN/musq_FAC)
      epcorr=epinv+2._dp*log(scale/facscale)

c--- for the case of virtual correction in the top quark decay,
c--- ('tdecay','ttdkay','Wtdkay') there are no extra initial-state
c--- contributions, so all these should be set to zero
      if  ( (kcase==ktdecay) .or. (kcase==kttdkay)
     & .or. (kcase==kWtdkay) .or. (kcase==ktt_ldk)
     & .or. (kcase==ktt_hdk) .or. (kcase==ktthWdk)
     & .or. (kcase==kdk_4ft) .or. (kcase==kttwldk)
     & .or. (kcase==kWHbbdk) .or. (kcase==kZHbbdk)
     & .or. (kcase==ktt_udk) .or. (kcase==kHWWdkW) ) then
        epcorr=0._dp
      endif

c--- for stop+b, splittings on light quark line produce a quark
      if ( (kcase==kqg_tbq) .or. (kcase==k4ftwdk)
     & .or.(kcase==kdk_4ft)) then
        epcorr=epinv+2._dp*log(renscale_L/facscale_L)
      endif

      AP(q,q,1)=+ason2pi*Cf*1.5_dp*epcorr
      AP(q,q,2)=+ason2pi*Cf*(-1._dp-z)*epcorr
      AP(q,q,3)=+ason2pi*Cf*2._dp/omz*epcorr
      AP(a,a,1)=+ason2pi*Cf*1.5_dp*epcorr
      AP(a,a,2)=+ason2pi*Cf*(-1._dp-z)*epcorr
      AP(a,a,3)=+ason2pi*Cf*2._dp/omz*epcorr

      AP(q,g,1)=0._dp
      AP(q,g,2)=ason2pi*Tr*(z**2+omz**2)*epcorr
      AP(q,g,3)=0._dp
      AP(a,g,1)=0._dp
      AP(a,g,2)=ason2pi*Tr*(z**2+omz**2)*epcorr
      AP(a,g,3)=0._dp

c--- modifications for running with mb>0
      if ( ((kcase==kW_twdk) .or. (kcase==kW_tndk))
     &  .and. (runstring(1:4) == 'mass')) then
      AP(q,g,2)=-ason2pi*Tr*(z**2+omz**2)*log(facscale**2/mbsq)
      AP(a,g,2)=AP(q,g,2)
      endif

      if ( ((kcase==kW_twdk) .or. (kcase==kW_tndk))
     &  .and. (nores) ) then
      AP(q,g,2)=0._dp
      AP(a,g,2)=AP(q,g,2)
      endif

c--- for stop+b, splittings on heavy quark line produce a gluon
      if ( (kcase==kqg_tbq) .or. (kcase==k4ftwdk)
     & .or.(kcase==kdk_4ft)) then
        epcorr=epinv+2._dp*log(renscale_H/facscale_H)
      endif

      AP(g,q,1)=0._dp
      AP(g,q,2)=ason2pi*Cf*(1._dp+omz**2)/z*epcorr
      AP(g,q,3)=0._dp
      AP(g,a,1)=0._dp
      AP(g,a,2)=ason2pi*Cf*(1._dp+omz**2)/z*epcorr
      AP(g,a,3)=0._dp

      AP(g,g,1)=+ason2pi*b0*epcorr
      AP(g,g,2)=+ason2pi*xn*2._dp*(1._dp/z+z*omz-2._dp)*epcorr
      AP(g,g,3)=+ason2pi*xn*2._dp/omz*epcorr

c--- for single top+b, make sure factors of alphas are correct
      if ( (kcase==kqg_tbq) .or. (kcase==k4ftwdk)
     & .or.(kcase==kdk_4ft)) then
        do is=1,3
c--- splittings on the light quark line make a (anti-)quark init. state
          AP(q,q,is)=AP(q,q,is)*(as_L/as)
          AP(a,a,is)=AP(a,a,is)*(as_L/as)
          AP(q,g,is)=AP(q,g,is)*(as_L/as)
          AP(a,g,is)=AP(a,g,is)*(as_L/as)
c--- splittings on the heavy quark line make a gluon init. state
          AP(g,g,is)=AP(g,g,is)*(as_H/as)
          AP(g,g,is)=AP(g,g,is)*(as_H/as)
          AP(g,q,is)=AP(g,q,is)*(as_H/as)
          AP(g,a,is)=AP(g,a,is)*(as_H/as)
      enddo
      endif

c--- AP functions for mixed EW contributions
      if     (kcase  ==  ktt_mix) then
         AP(:,:,:)=0._dp

      elseif (kcase  ==  ktwo_ew) then
         AP_mix(:,:,:,:) = 0._dp

         AP_mix(q,q,1:2,1)=+ason2pi*Cf*1.5_dp*epcorr
         AP_mix(q,q,1:2,2)=+ason2pi*Cf*(-1._dp-z)*epcorr
         AP_mix(q,q,1:2,3)=+ason2pi*Cf*2._dp/omz*epcorr
         AP_mix(a,a,1:2,1)=+ason2pi*Cf*1.5_dp*epcorr
         AP_mix(a,a,1:2,2)=+ason2pi*Cf*(-1._dp-z)*epcorr
         AP_mix(a,a,1:2,3)=+ason2pi*Cf*2._dp/omz*epcorr

         AP_mix(q,g,1:2,1)=0._dp
         AP_mix(q,g,1:2,2)=ason2pi*Tr*(z**2+omz**2)*epcorr
         AP_mix(q,g,1:2,3)=0._dp

         AP_mix(a,g,1:2,1)=0._dp
         AP_mix(a,g,1:2,2)=ason2pi*Tr*(z**2+omz**2)*epcorr
         AP_mix(a,g,1:2,3)=0._dp
      endif

c--- remove q -> g splittings for gg -> gam gam case
c--- (no longer necessary since inclusion of q+g->gam+gam+q contributions)
c      if (kcase == kgg2gam) then
c        AP(g,q,:)=0._dp
c        AP(g,a,:)=0._dp
c      endif

      if ( (kcase==kbq_tpq) .or. (kcase==kH_tjet)
     & .or.(kcase==kZ_tjet) .or. (kcase==kZ_tdkj)
     & .or. (kcase==kH_tdkj) .or. kcase==ktopanom ) then
        B1 = zip
        B2 = zip
      endif

      Q1 = zip
      Q2 = zip
      H1 = zip
      H2 = zip
      R1 = zip
      R2 = zip
      S1 = zip
      S2 = zip
      M1 = zip
      M2 = zip

      msqv = zip
      msq = zip
      msq_qq = zip
      msq_aq = zip
      msq_qa = zip

c--- Calculate the required matrix elements
      if     (kcase==kW_only) then
        call qqb_w(p,msq)
        call qqb_w_v(p,msqv)
        call qqb_w_z(p,z)
      elseif (kcase==kWln_ew) then
        call qqb_w(p,msq)
        call qqb_w_ew_v(p,msqv)
c        call qqb_w_ew_recola(p,msqv)
        call qqb_w_z_ew(p,z)
        AP(q,q,:)=AP(q,q,:)*abs(zesq)/gsq/Cf*Qu**2
        AP(a,a,:)=AP(a,a,:)*abs(zesq)/gsq/Cf*Qd**2
        AP(q,g,:)=AP(q,g,:)*abs(zesq)/gsq*(xn/tr)
        AP(a,g,:)=AP(a,g,:)*abs(zesq)/gsq*(xn/tr)
      elseif (kcase==kWln_aq) then
        call qqb_w(p,msq)
        msqv=0._dp
        call qphoton_wq_z(p,z)
        APqed(q,q,:)=AP(q,q,:)*abs(zesq)/gsq/Cf*Qu**2
        APqed(a,a,:)=AP(a,a,:)*abs(zesq)/gsq/Cf*Qd**2
        APqed(q,g,:)=AP(q,g,:)*abs(zesq)/gsq*(xn/tr)
        APqed(a,g,:)=AP(a,g,:)*abs(zesq)/gsq*(xn/tr)
      elseif (kcase==kW_1jet) then
        call qqb_w_g(p,msq)
        call qqb_w1jet_v(p,msqv)
        call qqb_w1jet_z(p,z)
      elseif (kcase==kWgaj_a) then
        call qqb_wgam(p,msq)
        msqv=0._dp
        call qphoton_wgamq_z(p,z)
        APqed(q,q,:)=AP(q,q,:)*abs(zesq)/gsq/Cf*Qu**2
        APqed(a,a,:)=AP(a,a,:)*abs(zesq)/gsq/Cf*Qd**2
        APqed(q,g,:)=AP(q,g,:)*abs(zesq)/gsq*(xn/tr)
        APqed(a,g,:)=AP(a,g,:)*abs(zesq)/gsq*(xn/tr)
      elseif (kcase==kWgajja) then
        call qqb_wgamg(p,msq)
        call qphoton_wgamq(p,msqqcd)
        msqv=0._dp
        call qphoton_wgamqg_z(p,z)
        APqed(q,q,:)=AP(q,q,:)*abs(zesq)/gsq/Cf*Qu**2
        APqed(a,a,:)=AP(a,a,:)*abs(zesq)/gsq/Cf*Qd**2
        APqed(q,g,:)=AP(q,g,:)*abs(zesq)/gsq*(xn/tr)
        APqed(a,g,:)=AP(a,g,:)*abs(zesq)/gsq*(xn/tr)
      elseif (kcase==kWga_ew) then
#ifdef HAVE_RECOLA
        call qqb_wgam(p,msq)
        call qqb_wgam_ew_recola(p,msqv)
        call qqb_wgam_z_ew(p,z)
        if (nwz == +1) then
        AP(q,q,:)=AP(q,q,:)*abs(zesq)/gsq/Cf*Qu**2
        AP(a,a,:)=AP(a,a,:)*abs(zesq)/gsq/Cf*Qd**2
        else
        AP(q,q,:)=AP(q,q,:)*abs(zesq)/gsq/Cf*Qd**2
        AP(a,a,:)=AP(a,a,:)*abs(zesq)/gsq/Cf*Qu**2
        endif
        AP(q,g,:)=AP(q,g,:)*abs(zesq)/gsq*(xn/tr)
        AP(a,g,:)=AP(a,g,:)*abs(zesq)/gsq*(xn/tr)
#else
        write(6,*) 'This process requires the code to be linked against Recola'
        write(6,*) 'and is not yet available in this version'
        stop
#endif
      elseif (kcase==kWgajew) then
#ifdef HAVE_RECOLA
        call qqb_wgamg(p,msq)
        call qqb_wgamg_ew_recola(p,msqv)
        call qqb_wgamg_z_ew(p,z)
        if (nwz == +1) then
        AP(q,q,:)=AP(q,q,:)*abs(zesq)/gsq/Cf*Qu**2
        AP(a,a,:)=AP(a,a,:)*abs(zesq)/gsq/Cf*Qd**2
        else
        AP(q,q,:)=AP(q,q,:)*abs(zesq)/gsq/Cf*Qd**2
        AP(a,a,:)=AP(a,a,:)*abs(zesq)/gsq/Cf*Qu**2
        endif
        AP(q,g,:)=AP(q,g,:)*abs(zesq)/gsq*(xn/tr)
        AP(a,g,:)=AP(a,g,:)*abs(zesq)/gsq*(xn/tr)
#else
        write(6,*) 'This process requires the code to be linked against Recola'
        write(6,*) 'and is not yet available in this version'
        stop
#endif
      elseif (kcase==kWgamma) then
        call qqb_wgam(p,msq)
        call qqb_wgam_v(p,msqv)
        call qqb_wgam_z(p,z)
      elseif (kcase==kWgajet) then
         call qqb_wgamg(p,msq)
         call qqb_wgamg_v(p,msqv)
         call qqb_wgamg_z(p,z)
      elseif (kcase==kW_cjet) then
        call qqb_w_cjet(p,msq)
        call qqb_w_cjet_v(p,msqv)
        call qqb_w_cjet_z(p,z)
      elseif (kcase==kWbfrmc) then
        call qqb_wbfromc(p,msq)
        call qqb_wbfromc_v(p,msqv)
        call qqb_wbfromc_z(p,z)
      elseif (kcase==kWbbmas) then
        call qqb_wbbm(p,msq)
        call qqb_wbbm_v(p,msqv)
        call qqb_wbbm_z(p,z)
      elseif (kcase==kWttmas) then
        call qqb_wbbm(p,msq)
        call qqb_wbbm_v(p,msqv)
        call qqb_wbbm_z(p,z)
      elseif (kcase==kWbbbar) then
        call qqb_wbb(p,msq)
        call qqb_wbb_v(p,msqv)
        call qqb_wbb_z(p,z)
      elseif (kcase==kW_2jet) then
        call qqb_wp2jetx_new(p,msq,mqq,ppmsqx,msqx_cs)
        call qqb_w2jet_v(p,msqv)
        call qqb_w2jet_z(p,z)
      elseif (kcase==kW_2gam) then
        stop
c         if (checkpiDpjk(p)) goto 999
c         call qqb_Waa(p,msq)
c         call qqb_Waa_v(p,msqv)
c         call qqb_Waa_z(p,z)
      elseif (kcase==kZ_only) then
        call qqb_z(p,msq)
        call qqb_z_v(p,msqv)
        call qqb_z_z(p,z)
        if (1 == 2) then
c Madgraph check: begin
        p(1,4)=-5.000000000000000E+002_dp
        p(1,1)=-0.000000000000000E+000_dp
        p(1,2)=-0.000000000000000E+000_dp
        p(1,3)=-5.000000000000000E+002_dp
        p(2,4)=-5.000000000000000E+002_dp
        p(2,1)=-0.000000000000000E+000_dp
        p(2,2)=-0.000000000000000E+000_dp
        p(2,3)= 5.000000000000000E+002_dp
        p(3,4)= 5.000000000000000E+002_dp
        p(3,1)= 4.510352131035046E+002_dp
        p(3,2)= 6.047968715179463E+001_dp
        p(3,3)=-2.071459485065961E+002_dp
        p(4,4)= 5.000000000000000E+002_dp
        p(4,1)=-4.510352131035046E+002_dp
        p(4,2)=-6.047968715179463E+001_dp
        p(4,3)= 2.071459485065961E+002_dp
        call qqb_z(p,msq)
        epinv=zip
        epinv2=zip
        call qqb_z_v(p,msq0)
        epinv=one
        epinv2=one
        call qqb_z_v(p,msq1)
        epinv=-one
        epinv2=-one
        call qqb_z_v(p,msqm1)
c remember extra log from expanding (Qsq/musq)^ep, with Qsq=1000^2 GeV
        xlog=two*log(1000._dp/scale)
        bit2(:,:)=(msq1(:,:)+msqm1(:,:)-two*msq0(:,:))/two
        bit1(:,:)=(msq1(:,:)-msqm1(:,:))/two+xlog*bit2(:,:)
        bit0(:,:)=msq0(:,:)+xlog*(bit1(:,:)-xlog*bit2(:,:))+xlog**2/two*bit2(:,:)
c translation from DR (MCFM) to tH-V (MG)
        bit0(:,:)=bit0(:,:)-CF*ason2pi*msq(:,:)
        write(6,*) 'msq(2,-2) ',msq(2,-2)
        write(6,*)
        write(6,*) 'msqv (2,-2)     ',bit0(2,-2)
        write(6,*) 'msqv (2,-2) ep  ',bit1(2,-2)
        write(6,*) 'msqv (2,-2) ep^2',bit2(2,-2)
        write(6,*)
        write(6,*) 'msqv (2,-2)      /(msq(2,-2)*ason2pi)',bit0(2,-2)/(msq(2,-2)*ason2pi)
        write(6,*) 'msqv (2,-2) ep   /(msq(2,-2)*ason2pi)',bit1(2,-2)/(msq(2,-2)*ason2pi)
        write(6,*) 'msqv (2,-2) ep^2 /(msq(2,-2)*ason2pi)',bit2(2,-2)/(msq(2,-2)*ason2pi)
        stop
c Madgraph check: end
        endif
      elseif (kcase==kZ_1jet) then
        call qqb_z1jet(p,msq)
        call qqb_z1jet_v(p,msqv)
        call qqb_z1jet_z(p,z)
        if (1 == 2) then
c Madgraph check: begin
        p(1,4)=-5.000000000000000E+002_dp
        p(1,1)=-0.000000000000000E+000_dp
        p(1,2)=-0.000000000000000E+000_dp
        p(1,3)=-5.000000000000000E+002_dp
        p(2,4)=-5.000000000000000E+002_dp
        p(2,1)=-0.000000000000000E+000_dp
        p(2,2)=-0.000000000000000E+000_dp
        p(2,3)= 5.000000000000000E+002_dp
        p(3,4)= 1.850618930902392E+002_dp
        p(3,1)= -5.101444150706673E+001_dp
        p(3,2)= -8.234967641309825E+001_dp
        p(3,3)= -1.576831057105458E+002_dp
        p(4,4)= 4.192381120256421E+002_dp
        p(4,1)= -1.165087601320307E+002_dp
        p(4,2)= 4.027108985766566E+002_dp
        p(4,3)= -3.199305378265266E+000_dp
        p(5,4)= 3.956999948841188E+002_dp
        p(5,1)= 1.675232016390975E+002_dp
        p(5,2)= -3.203612221635584E+002_dp
        p(5,3)= 1.608824110888112E+002_dp
        call qqb_z1jet(p,msq)
        epinv=zip
        epinv2=zip
        call qqb_z1jet_v(p,msq0)
        epinv=one
        epinv2=one
        call qqb_z1jet_v(p,msq1)
        epinv=-one
        epinv2=-one
        call qqb_z1jet_v(p,msqm1)
c remember extra log from expanding (Qsq/musq)^ep, with Qsq=1000^2 GeV
        xlog=two*log(1000._dp/scale)
        bit2(:,:)=(msq1(:,:)+msqm1(:,:)-two*msq0(:,:))/two
        bit1(:,:)=(msq1(:,:)-msqm1(:,:))/two+xlog*bit2(:,:)
        bit0(:,:)=msq0(:,:)+xlog*(bit1(:,:)-xlog*bit2(:,:))+xlog**2/two*bit2(:,:)
c now do UV subtraction that has been removed in virtual routine
        bit1(:,:)=bit1(:,:)-msq(:,:)*ason2pi*xn*(11._dp-2._dp*real(nflav-1,dp)/xn)/6._dp
        bit0(:,:)=bit0(:,:)-msq(:,:)*ason2pi*xn*(-1._dp)/6._dp
c translation from DR (MCFM) to tH-V (MG)
        bit0(:,:)=bit0(:,:)-(CF+XN/six)*ason2pi*msq(:,:)
        write(6,*) 'msq(2,-2) ',msq(2,-2)
        write(6,*)
        write(6,*) 'msqv (2,-2)     ',bit0(2,-2)
        write(6,*) 'msqv (2,-2) ep  ',bit1(2,-2)
        write(6,*) 'msqv (2,-2) ep^2',bit2(2,-2)
        write(6,*)
        write(6,*) 'msqv (2,-2)      /(msq(2,-2)*ason2pi)',bit0(2,-2)/(msq(2,-2)*ason2pi)
        write(6,*) 'msqv (2,-2) ep   /(msq(2,-2)*ason2pi)',bit1(2,-2)/(msq(2,-2)*ason2pi)
        write(6,*) 'msqv (2,-2) ep^2 /(msq(2,-2)*ason2pi)',bit2(2,-2)/(msq(2,-2)*ason2pi)
        stop
c Madgraph check: end
        endif
      elseif (kcase==kZ_2jet) then
c        call qqb_z2jetx(p,msq,mqq,msqx,msqx_cs)
        call qqb_z2jetx_new(p,msq,mqq,ppmsqx,msqx_cs)
        call qqb_z2jet_v(p,msqv)
        if (onlyaxial) then
          msq=0._dp
          mqq=0._dp
          ppmsqx=0._dp
          msqx_cs=0._dp
          msq_cs=0._dp
        else
          call qqb_z2jetx_new(p,msq,mqq,ppmsqx,msqx_cs)
        call qqb_z2jet_z(p,z)
        endif
      elseif (kcase==kZgamma) then
        ! Have not implemented the poles for Zgamma, it wouldn't test
        ! anything further..
        epinv=0._dp
        epinv2=0._dp
        nshot = 3
        call set_anomcoup(p)
        if (anomtgc .eqv. .false.) then
c Use old routines if in the SM since they're faster
          call qqb_zgam(p,msq)
          call qqb_zgam_v(p,msqv)
        else
          call qqb_zgam_v_new(p,msqv,msq)
        endif
        call qqb_zgam_z(p,z)
c        call gg_zgam(p,msqv(0,0)) ! additional gluon-gluon contribution
      elseif (kcase==kZ_2gam) then
        call qqb_zaa(p,msq)
        call qqb_zaa_v(p,msqv)
        call qqb_zaa_z(p,z)
      elseif (kcase==kZgajet) then
        call set_anomcoup(p)
c       call qqb_zaj(p,msq)
c       call qqb_zaj_v(p,msqv)
        call qqb_zaj_v_new(p,msqv,msq)
        call qqb_zaj_z(p,z)
      elseif (kcase==kZbbbar) then
        call qqb_zbb(p,msq)
        call qqb_zbb_v(p,msqv)
        call qqb_zbb_z(p,z)
      elseif (kcase==kWWqqbr) then
c gg contributions now included at NNLO
c        if (ggonly) then ! special catch for gg->WW piece only
c          msq(:,:)=0._dp
c          msqv(:,:)=0._dp
c        else
          call qqb_ww(p,msq)
          call qqb_ww_v(p,msqv)
          call qqb_ww_z(p,z)
          if (origkpart==ktodk) then
            call dkqqb_ww_v(p,msqvdk)
            do j=-nf,nf
            do k=-nf,nf
              msqv(j,k)=msqv(j,k)+msqvdk(j,k)
            enddo
            enddo
          endif
c        endif
c        if (omitgg .eqv. .false.) then ! do not compute if omitgg true
c          call gg_ww_int(p,msqvdk) ! additional gluon-gluon contribution
c          msqv(0,0)=msqvdk(0,0)
c        endif
      elseif (kcase==kWWqqdk) then
        msq(:,:)=0._dp
        call dkqqb_ww_v(p,msqv)
      elseif (kcase==kWZbbar) then
        call qqb_wz(p,msq)
        call qqb_wz_v(p,msqv)
        call qqb_wz_z(p,z)
      elseif (kcase==kZZlept) then
c gg contributions now included at NNLO
c        if (ggonly) then ! special catch for gg->ZZ piece only
c          msq(:,:)=0._dp
c          msqv(:,:)=0._dp
c        else
          call qqb_zz(p,msq)
          call qqb_zz_v(p,msqv)
          call qqb_zz_z(p,z)
c        endif
c        if (omitgg .eqv. .false.) then ! do not compute if omitgg true
c          call gg_ZZ_all(p,msqvdk) ! additional gluon-gluon contribution
c          msqv(0,0)=msqvdk(0,0)
c        endif
      elseif (kcase==kVVlept) then
          call qqb_vv(p,msq)
          call qqb_vv_v(p,msqv)
          call qqb_vv_z(p,z)
      elseif (kcase==kWHbbar) then
        call qqb_wh(p,msq)
        call qqb_wh_v(p,msqv)
        call qqb_wh_z(p,z)
        if (origkpart==ktodk) then
          if (mb < 1.e-6_dp) then
            call dkqqb_wh_v_massless(p,msqvdk)
          else
            call dkqqb_wh_v(p,msqvdk)
          endif
          msqv(:,:)=msqv(:,:)+msqvdk(:,:)
c--- correct for scaling with NLO partial width
          msqv(:,:)=msqv(:,:)+msq(:,:)*(GamHbb0/GamHbb1-one)
        endif

      elseif (kcase==kWHbbdk) then
        if (mb < 1.e-6_dp) then
          call dkqqb_wh_v_massless(p,msqv)
        else
          call dkqqb_wh_v(p,msqv)
        endif
c--- correct for scaling with NLO partial width
        call qqb_wh(p,msq)
        msq(:,:)=msq(:,:)*(GamHbb0/GamHbb1-one)
      elseif (kcase==kWH1jet) then
         if(toponly) then
            call qqb_WH1jet_v(p,msqv)
            msq(:,:)=zip
         else
            call qqb_WH1jet(p,msq)
            call qqb_WH1jet_v(p,msqv)
            call qqb_WH1jet_z(p,z)
         endif
      elseif (kcase==kWH__WW) then
        call qqb_wh_ww(p,msq)
        call qqb_wh_ww_v(p,msqv)
        call qqb_wh_z(p,z) ! nb: the same as above
      elseif (kcase==kWH__ZZ) then
        call qqb_wh_zz(p,msq)
        call qqb_wh_zz_v(p,msqv)
        call qqb_wh_z(p,z) ! nb: the same as above
      elseif (kcase==kWHgaga) then
        call qqb_wh_gaga(p,msq)
        call qqb_wh_gaga_v(p,msqv)
        call qqb_wh_z(p,z) ! nb: the same as above
      elseif (kcase==kZHbbar) then
        call qqb_zh(p,msq)
        call qqb_zh_v(p,msqv)
        call qqb_zh_z(p,z)
        if (origkpart==ktodk) then
          if (mb < 1.e-6_dp) then
            call dkqqb_zh_v_massless(p,msqvdk)
          else
            call dkqqb_zh_v(p,msqvdk)
          endif
          msqv(:,:)=msqv(:,:)+msqvdk(:,:)
c--- correct for scaling with NLO partial width
          msqv(:,:)=msqv(:,:)+msq(:,:)*(GamHbb0/GamHbb1-one)
        endif

      elseif (kcase==kZHbbdk) then
         if (mb < 1.e-6_dp) then
            call dkqqb_zh_v_massless(p,msqv)
         else
            call dkqqb_zh_v(p,msqv)
         endif
c--- correct for scaling with NLO partial width
         call qqb_zh(p,msq)
         msq(:,:)=msq(:,:)*(GamHbb0/GamHbb1-one)
      elseif (kcase==kZHgaga) then
        call qqb_zh_gaga(p,msq)
        call qqb_zh_gaga_v(p,msqv)
        call qqb_zh_z(p,z)      ! nb: the same as above
      elseif (kcase==kZH__WW) then
        call qqb_zh_ww(p,msq)
        call qqb_zh_ww_v(p,msqv)
        call qqb_zh_z(p,z) ! nb: the same as above
      elseif (kcase==kZH__ZZ) then
        call qqb_zh_zz(p,msq)
        call qqb_zh_zz_v(p,msqv)
        call qqb_zh_z(p,z) ! nb: the same as above
      elseif (kcase==kZH1jet) then
        if (toponly) then
          msq=zip
         call qqb_ZH1jet_v(p,msqv)
        else
          call qqb_ZH1jet(p,msq)
          call qqb_ZH1jet_v(p,msqv)
          call qqb_ZH1jet_z(p,z)
        endif
      elseif (kcase==ktwojet) then
        call qqb_twojet(p,msq)
        call qqb_twojet_v(p,msqv)
        write(6,*) 'This process not fully implemented'
        stop
      elseif (kcase==kggfus0) then
        call gg_h(p,msq)
        call gg_h_v(p,msqv)
        call gg_h_z(p,z)
      elseif (kcase==kHigaga) then
        call gg_hgamgam(p,msq)
        call gg_hgamgam_v(p,msqv)
        call gg_hgamgam_z(p,z)
      elseif (kcase==kHi_Zga) then
        call gg_hzgam(p,msq)
        call gg_hzgam_v(p,msqv)
        call gg_hzgam_z(p,z)
      elseif (kcase==kHWW_4l) then
        call qqb_hww(p,msq)
        call qqb_hww_v(p,msqv)
        call qqb_hww_z(p,z)
      elseif (kcase==kHWW2lq) then
        call qqb_hww(p,msq)
        call qqb_hww_v(p,msqv)
        call qqb_hww_z(p,z)
        if (origkpart==ktodk) then
          call dkqqb_hww_v(p,msqvdkW)
          msqv(:,:)=msqv(:,:)+msqvdkW(:,:)
        endif
      elseif (kcase==kHWWdkW) then
        msq(:,:)=0._dp
        call dkqqb_hww_v(p,msqv)
      elseif (kcase==kHZZ_4l) then
        call qqb_hzz(p,msq)
        call qqb_hzz_v(p,msqv)
        call qqb_hzz_z(p,z)
      elseif (kcase==kH_1jet) then
        call qqb_hg(p,msq)
        call qqb_hg_v(p,msqv)
        call qqb_hg_z(p,z)
      elseif (kcase==kHWWjet) then
        call gg_hWWg(p,msq)
        call gg_hWWg_v(p,msqv)
        call gg_hWWg_z(p,z)
      elseif (kcase==kHZZjet) then
        call gg_hZZg(p,msq)
        call gg_hZZg_v(p,msqv)
        call gg_hZZg_z(p,z)
      elseif (kcase==ktwo_ew) then
        wwidth_save=wwidth
        zwidth_save=zwidth
        wwidth=zip
        zwidth=zip
        call qqb_twojet_mix(p,msq)
        call qqb_twojet_ew_exact(p,msqv)
        call qqb_twojet_mix_z(p,z)
c        msq_mix(:,:,:)=zip ! uncomment to remove mixed terms if desired
        call qqb_twojet(p,msq_noew) ! no EW corrections
c---- include the following three lines to compute corrections wrt full LO
c        call qqb_twojet_ew(p,msq)
c        msq_noew=msq_noew+msq
c        msq=zip
        wwidth=wwidth_save
        zwidth=zwidth_save
      elseif (kcase==kdirgam) then
        call qqb_dirgam(p,msq)
        call qqb_dirgam_v(p,msqv)
        call qqb_dirgam_z(p,z)
      elseif (kcase==khflgam) then
        call qqb_hflgam(p,msq)
        call qqb_hflgam_v(p,msqv)
        call qqb_hflgam_z(p,z)
      elseif (kcase==kgamgam) then
        call qqb_gamgam(p,msq)
        call qqb_gamgam_v(p,msqv)
        call qqb_gamgam_z(p,z)
      elseif (kcase==kgg2gam) then
        call gg_2gam(p,msq)
        call gg_2gam_v(p,msqv)
        call gg_2gam_z(p,z)
      elseif (kcase==kgmgmjt) then
         call qqb_gmgmjt(p,msq)
         call qqb_gmgmjt_v(p,msqv)
         call qqb_gmgmjt_z(p,z)
      elseif (kcase==kgam_2j) then
c         call checksym_gam2jet_v(p)
         call qqb_gam2jx_new(p,msq,mqq,ppmsqx,msqx_cs)
         call qqb_gam2j_v(p,msqv)
         call qqb_gam2j_z(p,z)
      elseif (kcase==ktrigam) then
        call qqb_trigam(p,msq)
        call qqb_trigam_v(p,msqv)
        call qqb_trigam_z(p,z)
      elseif (kcase==kfourga) then
        call qqb_fourgam(p,msq)
        call qqb_fourgam_v(p,msqv)
        call qqb_fourgam_z(p,z)
      elseif ((kcase==ktt_bbl) .or. (kcase==ktt_bbh)) then
c Sanity check for tiny pt(top)/mt -- otherwise numerical Gram det. problems
        if (ptthree(3,4,5,p)/mt < cutoff) goto 999
        call qqb_QQbdk(p,msq)
        call qqb_QQbdk_v(p,msqv)
        call qqb_QQbdk_z(p,z)
        if (origkpart==ktodk) then
          call dkqqb_QQb_v(p,msqvdk)
          call dkWqqb_QQb_v(p,msqvdkW)
          do j=-nf,nf
          do k=-nf,nf
            msqv(j,k)=msqv(j,k)+msqvdk(j,k)+msqvdkW(j,k)
          enddo
          enddo
        endif
      elseif ((kcase==ktt_ldk) .or. (kcase==ktt_hdk)) then
        msq(:,:)=0._dp
        call dkqqb_QQb_v(p,msqv)
      elseif (kcase==ktthWdk) then
        msq(:,:)=0._dp
        call dkWqqb_QQb_v(p,msqv)
      elseif (kcase==ktt_bbu) then
        call qqb_QQbdku(p,msq)
        call qqb_QQbdku_v(p,msqv)
        call qqb_QQbdku_z(p,z)
        if (origkpart==ktodk) then
          call dkuqqb_QQb_v(p,msqvdk)
          do j=-nf,nf
          do k=-nf,nf
            msqv(j,k)=msqv(j,k)+msqvdk(j,k)
          enddo
          enddo
        endif
      elseif (kcase==ktt_udk) then
        msq(:,:)=0._dp
        call dkuqqb_QQb_v(p,msqv)
      elseif (kcase==ktt_mix) then
        call qqb_QQb_mix(p,msq)
        call qqb_QQb_mix_v(p,msqv)
        call qqb_QQb_mix_z(p,z)
        call qqb_QQb(p,msq_noew) ! matrix elements with no EW corrections
      elseif ((kcase==ktt_tot)
     &   .or. (kcase==kbb_tot)
     &   .or. (kcase==kcc_tot)) then
        call qqb_QQb(p,msq)
        call qqb_QQb_v(p,msqv)
        call qqb_QQb_z(p,z)
      elseif (kcase==ktopanom) then
        ! this has been debugged, but ocasionally there are unstable points
        ! due to the strong cancellation effects with epinv /= 0
        bbfrac = 0._dp
        msqv_heavy = 0._dp
        msqv_light = 0._dp

        call singletop2_scale_reset()
        call singletop2_scale_setup(p)

        badpoint = .false.
        call singletop2_tree(p,msq)
        call singletop2_virt(p,msqv_light,.true.,.false.)
        call singletop2_virt(p,msqv_heavy,.false.,.true.)
        msqv = msqv_light + msqv_heavy
        call singletop2_z(p,z)

        if (badpoint .eqv. .true.) then
            msq = 0._dp
            msqv = 0._dp
            msqv_light = 0._dp
            msqv_heavy = 0._dp
            goto 999
        endif
      elseif (kcase==kbq_tpq) then
        call bq_tpq(p,msq)
        call bq_tpq_v(p,msqv)
        call bq_tpq_z(p,z)
        if (origkpart==ktodk) then
          call bq_tpq_vdk(p,msqvdk)
          do j=-nf,nf
          do k=-nf,nf
            msqv(j,k)=msqv(j,k)+msqvdk(j,k)
          enddo
          enddo
        endif

      elseif (kcase==kttdkay) then
        msq(:,:)=0._dp
        call bq_tpq_vdk(p,msqv)
      elseif (kcase==kt_bbar) then
        call qqb_tbbdk(p,msq)
        call qqb_tbbdk_v(p,msqv)
        call qqb_tbbdk_z(p,z)
        if (origkpart==ktodk) then
          call dkqqb_tbbdk_v(p,msqvdk)
          do j=-nf,nf
          do k=-nf,nf
            msqv(j,k)=msqv(j,k)+msqvdk(j,k)
          enddo
          enddo
        endif
      elseif (kcase==ktdecay) then
        msq(:,:)=0._dp
        call dkqqb_tbbdk_v(p,msqv)
      elseif (kcase==kW_tndk) then
        call qqb_w_tndk(p,msq)
        call qqb_w_tndk_v(p,msqv)
        call qqb_w_tndk_z(p,z)
      elseif (kcase==kW_twdk) then
        call qqb_w_twdk(p,msq)
        call qqb_w_twdk_v(p,msqv)
        call qqb_w_twdk_z(p,z)
        if (origkpart==ktodk) then
          call dkqqb_w_twdk_v(p,msqvdk)
          msqv(:,:)=msqv(:,:)+msqvdk(:,:)
        endif
      elseif (kcase==kWtdkay) then
        msq(:,:)=0._dp
        call dkqqb_w_twdk_v(p,msqv)
      elseif (kcase==kqq_ttw) then
        call qqb_ttw(p,msq)
        call qqb_ttw_v(p,msqv)
        call qqb_ttw_z(p,z)
        if (origkpart==ktodk) then
          call dkqqb_ttw_v(p,msqvdk)
          msqv(:,:)=msqv(:,:)+msqvdk(:,:)
        endif
      elseif (kcase==kttwldk) then
        msq(:,:)=0._dp
        call dkqqb_ttw_v(p,msqv)
      elseif (kcase==kggfus1) then
        if (expansionorder < 2) then
            Wilsonorder = 0
        elseif (expansionorder == 2) then
            Wilsonorder = 1
        endif
        call gg_hg(p,msq)
        call gg_hg_v(p,msqv)
        call gg_hg_z(p,z)
        if (expansionorder > 0) then
            Wilsonorder = expansionorder
            call gg_hg(p,msqvdk)
            msqv = msqv + msqvdk - msq
            msqvdk = zip
        endif
      elseif (kcase==kHi_Zaj) then
        call gg_hzgamg(p,msq)
        call gg_hg_zgam_v(p,msqv)
        call gg_hg_zgam_z(p,z)
      elseif (kcase==khjetma) then
        call qqb_higgs(p,msq)
        call hjetmass_v(p,msqv)
        call hjetmass_z(p,z)
      elseif (kcase==kHgagaj) then
        call gg_hgagag(p,msq)
        call gg_hgagag_v(p,msqv)
        call gg_hg_z(p,z) ! nb: the same as above
      elseif (kcase==kggfus2) then
        Wilsonorder = 0
        call gg_hgg(p,msq)
        call gg_hgg_v(p,msqv)
        call gg_hgg_z(p,z)
      elseif (kcase==kgagajj) then
        Wilsonorder = 0
        call gg_hgg(p,msq)
        call gg_hgg_v(p,msqv)
        call gg_hgg_z(p,z)
      elseif (kcase==kHWW2jt) then
        call gg_hWWgg(p,msq)
        call gg_hWWgg_v(p,msqv)
        call gg_hWWgg_z(p,z)
      elseif (kcase==kHZZ2jt) then
        call gg_hZZgg(p,msq)
        call gg_hZZgg_v(p,msqv)
        call gg_hZZgg_z(p,z)
      elseif (kcase==kqq_Hqq) then
        call VV_hqq(p,msq)
        call VV_hqq_v(p,msqv)
        call VV_hqq_z(p,z)
      elseif (kcase==kqq_Hgg) then
        call VV_Hgaga(p,msq)
        call VV_Hgaga_v(p,msqv)
        call VV_Hgaga_z(p,z)
      elseif (kcase==kqq_HWW) then
        call VV_HWW(p,msq)
        call VV_HWW_v(p,msqv)
        call VV_HWW_z(p,z)
      elseif (kcase==kqq_HZZ) then
        call VV_HZZ(p,msq)
        call VV_HZZ_v(p,msqv)
        call VV_HZZ_z(p,z)
      elseif (kcase==kqg_tbq) then
        call qg_tbq(p,msq)
        call qg_tbq_v(p,msqv)
        call qg_tbq_z(p,z)
      elseif (kcase==k4ftwdk) then
        call qg_tbqdk(p,msq)
        call qg_tbqdk_v(p,msqv)
        call qg_tbqdk_z(p,z)
        if (origkpart==ktodk) then
          call dkqg_tbqdk_v(p,msqvdk)
          msqv(:,:)=msqv(:,:)+msqvdk(:,:)
        endif
      elseif (kcase==kdk_4ft) then
        msq(:,:)=0._dp
        call dkqg_tbqdk_v(p,msqv)
      elseif (kcase==kqq_tbg) then
c--- do not include initial-state subtractions since we are only
c--- calculating corrections on the heavy quark line (for now)
        do j=-1,1
        do k=-1,1
        AP(j,k,1)=0._dp
        AP(j,k,2)=0._dp
        AP(j,k,3)=0._dp
        enddo
        enddo
        call qq_tbg(p,msq)
        call qq_tbg_v(p,msqv)
        call qq_tbg_z(p,z)
      elseif (kcase==kH_tjet) then
        call qq_tchan_htq(p,msq)
        call qq_tchan_htq_v(p,msqv)
        if (pvbadpoint) then
          goto 999
        endif
        call qq_tchan_htq_z(p,z)
      elseif (kcase==kH_tdkj) then
        call qq_tchan_htq_dk(p,msq)
        call qq_tchan_htq_dk_v(p,msqv)
        if (pvbadpoint) then
          goto 999
        endif
        call qq_tchan_htq_dk_z(p,z)
      elseif (kcase==kZ_tjet) then
        call qq_tchan_ztq(p,msq)
        call qq_tchan_ztq_v(p,msqv)
        if (pvbadpoint) then
          goto 999
        endif
        call qq_tchan_ztq_z(p,z)
      elseif (kcase==kZ_tdkj) then
        call qq_tchan_ztq_dk(p,msq)
        call qq_tchan_ztq_dk_v(p,msqv)
        if (pvbadpoint) then
          goto 999
        endif
        call qq_tchan_ztq_dk_z(p,z)

c -- NLO corrections to decay ME -- can include flag for this at some stage...
c         do j=-nf,nf
c         do k=-nf,nf
c         msqv=msqv+msqvdk
c         enddo
c         enddo
c         call qq_tchan_ztq_z(p,z)

      elseif (kcase==kepem3j) then
c--- do not include initial-state subtractions since we are only
c--- calculating corrections on the quark line
        do j=-1,1
        do k=-1,1
        AP(j,k,1)=0._dp
        AP(j,k,2)=0._dp
        AP(j,k,3)=0._dp
        enddo
        enddo
        call epem3j(p,msq)
        call epem3j_v(p,msqv)
        call epem3j_z(p,z)
      elseif (kcase==kgQ__ZQ) then
        call gQ_zQ(p,msq)
        call gQ_zQ_v(p,msqv)
        call gQ_zQ_z(p,z)
      elseif (kcase==kZ_bjet) then
        call qqb_zbjet(p,msq)
        call qqb_zbjet_v(p,msqv)
        call qqb_zbjet_z(p,z)
      elseif (kcase==kW_bjet) then
        call qqb_Wbjet(p,msq)
        call qqb_Wbjet_v(p,msqv)
        call qqb_Wbjet_z(p,z)
        APqg_mass=-ason2pi*Tr*(z**2+omz**2)*log(facscale**2/mbsq)
      elseif (kcase==kWcsbar) then
        call qqb_w(p,msq)
        call qqb_w_v(p,msqv)
        call qqb_w_z(p,z)
      elseif (kcase==kWcs_ms) then
c--- massive subtraction term only
        call qqb_w(p,msq)
      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0._dp
      enddo
      enddo
        do j=-1,1
        do k=-1,1
        AP(j,k,1)=0._dp
        AP(j,k,2)=0._dp
        AP(j,k,3)=0._dp
        enddo
        enddo
        AP(q,g,2)=-ason2pi*Tr*(z**2+omz**2)*log(facscale**2/mcsq)

      elseif (kcase==kdm_jet) then
         call qqb_dm_monojet_v(p,msqv)
         call qqb_dm_monojet(p,msq)
         if(dm_mediator=='gluonO') then
            call gg_dm_monojet_z(p,z)
         else
            call qqb_dm_monojet_z(p,z)
         endif

      elseif(kcase==kdm_gam) then
         call qqb_dm_monophot_v(p,msqv)
         call qqb_dm_monophot(p,msq)
         call qqb_dm_monophot_z(p,z)

      elseif (kcase==kWW_jet) then
        if (WWjetcheckpiDpjk(p)) goto 999 ! for numerical safety
        call qqb_wwg(p,msq)
        call qqb_wwg_v(p,msqv)
        call qqb_wwg_z(p,z)

      elseif (kcase==kWZ_jet) then
        if (WWjetcheckpiDpjk(p)) goto 999 ! for numerical safety
        call qqb_wzg(p,msq)
        call qqb_wzg_v(p,msqv)
        call qqb_wzg_z(p,z)

      elseif (kcase==kZZ_jet) then
        if (WWjetcheckpiDpjk(p)) goto 999 ! for numerical safety
        call qqb_zzg(p,msq)
        call qqb_zzg_v(p,msqv)
        call qqb_zzg_z(p,z)

      else
        write (*,*) "Undefined kcase in virtint.f"
        write(6,*) 'Abort in virtint'
      stop
      endif

c--- explicitly remove factor of LO if we are only interested in coefficient
      if (coeffonly) then
        msqv(:,:)=msqv(:,:)-msq(:,:)
      endif

      !currentPDF = 0
      currentNd = 0

c--- initialize a PDF set here, if calculating errors
  777 continue
      xmsq=0._dp
      fx1z(:)=0._dp
      fx2z(:)=0._dp

      if (z > xx(1)) x1onz=xx(1)/z
      if (z > xx(2)) x2onz=xx(2)/z
      xmsq_noew=0._dp

c--- calculate PDF's
c--- Note: if computing PDF errors then the whole block of code should
c--- be OMP critical;  to avoid unnecessary duplication of code the
c--- block should be shuffled off to a routine.  Current solution is lazy!!!
      if (kcase==ktopanom) then
        ! assume singletop2_scale_setup has already been called above
        call fdist(ih1,xx(1),facscale_beam1_islight_onlight, singletop2_pdfs(i_beam1_light,:),1)
        call fdist(ih2,xx(2),facscale_beam2_isheavy_onheavy, singletop2_pdfs(i_beam2_heavy,:),2)

        call fdist(ih1,xx(1),facscale_beam1_isheavy_onheavy, singletop2_pdfs(i_beam1_heavy,:),1)
        call fdist(ih2,xx(2),facscale_beam2_islight_onlight, singletop2_pdfs(i_beam2_light,:),2)

        if (z > xx(1)) then
            call fdist(ih1,x1onz,facscale_beam1_islight_onlight, singletop2_pdfs(i_beam1_light_z,:),1)
        endif

        if (z > xx(2)) then
            call fdist(ih2,x2onz,facscale_beam2_isheavy_onheavy, singletop2_pdfs(i_beam2_heavy_z,:),2)
        endif

        if (z > xx(1)) then
            call fdist(ih1,x1onz,facscale_beam1_isheavy_onheavy, singletop2_pdfs(i_beam1_heavy_z,:),1)
        endif

        if (z > xx(2)) then
            call fdist(ih2,x2onz,facscale_beam2_islight_onlight, singletop2_pdfs(i_beam2_light_z,:),2)
        endif

        fx1 = singletop2_pdfs(i_beam1_light,:)
        fx2 = singletop2_pdfs(i_beam2_heavy,:)

        fx1z = singletop2_pdfs(i_beam1_light_z,:)
        fx2z = singletop2_pdfs(i_beam2_heavy_z,:)

      else
          if ((maxPDFsets > 0).and. bin) then
        if ( (kcase==kqg_tbq) .or. (kcase==k4ftwdk)
     &   .or.(kcase==kdk_4ft)) then
c--- for single top + b, make sure to use two different scales
          call fdist(ih1,xx(1),facscale_H,fx1_H,1)
          call fdist(ih2,xx(2),facscale_H,fx2_H,2)
          call fdist(ih1,xx(1),facscale_L,fx1_L,1)
          call fdist(ih2,xx(2),facscale_L,fx2_L,2)
          fx1z_H(:)=0._dp
          fx2z_H(:)=0._dp
          fx1z_L(:)=0._dp
          fx2z_L(:)=0._dp
          if (z > xx(1)) then
            call fdist(ih1,x1onz,facscale_H,fx1z_H,1)
            call fdist(ih1,x1onz,facscale_L,fx1z_L,2)
          endif
          if (z > xx(2)) then
            call fdist(ih2,x2onz,facscale_H,fx2z_H,1)
            call fdist(ih2,x2onz,facscale_L,fx2z_L,2)
          endif
        else
c--- usual case
          call fdist(ih1,xx(1),facscale,fx1,1)
          if (kcase == kWgajja) call photonpdffix(fx1,fxa)
          call fdist(ih2,xx(2),facscale,fx2,2)
          if (kcase == kWgajja) call photonpdffix(fx2,fxa)
          if (z > xx(1)) then
            call fdist(ih1,x1onz,facscale,fx1z,1)
            if (kcase == kWgajja) call photonpdffix(fx1z,fxa)
          endif
          if (z > xx(2)) then
            call fdist(ih2,x2onz,facscale,fx2z,2)
            if (kcase == kWgajja) call photonpdffix(fx2z,fxa)
          endif
        endif
      else
        if ( (kcase==kqg_tbq) .or. (kcase==k4ftwdk)
     &   .or.(kcase==kdk_4ft)) then
c--- for single top + b, make sure to use two different scales
          call fdist(ih1,xx(1),facscale_H,fx1_H,1)
          call fdist(ih2,xx(2),facscale_H,fx2_H,2)
          call fdist(ih1,xx(1),facscale_L,fx1_L,1)
          call fdist(ih2,xx(2),facscale_L,fx2_L,2)
          fx1z_H(:)=0._dp
          fx2z_H(:)=0._dp
          fx1z_L(:)=0._dp
          fx2z_L(:)=0._dp
          if (z > xx(1)) then
            call fdist(ih1,x1onz,facscale_H,fx1z_H,1)
            call fdist(ih1,x1onz,facscale_L,fx1z_L,2)
          endif
          if (z > xx(2)) then
            call fdist(ih2,x2onz,facscale_H,fx2z_H,1)
            call fdist(ih2,x2onz,facscale_L,fx2z_L,2)
          endif
        else
c--- usual case
          call fdist(ih1,xx(1),facscale,fx1,1)
          if (kcase == kWgajja) call photonpdffix(fx1,fxa)
          call fdist(ih2,xx(2),facscale,fx2,2)
          if (kcase == kWgajja) call photonpdffix(fx2,fxa)
          if (z > xx(1)) then
            call fdist(ih1,x1onz,facscale,fx1z,1)
            if (kcase == kWgajja) call photonpdffix(fx1z,fxa)
          endif
          if (z > xx(2)) then
            call fdist(ih2,x2onz,facscale,fx2z,2)
            if (kcase == kWgajja) call photonpdffix(fx2z,fxa)
          endif
        endif

      endif
      endif

      do j=-nflav,nflav
      do k=-nflav,nflav

      if ((.not. selectpdfs(1,j)) .or. (.not. selectpdfs(2,k))) cycle

      if ((kcase==kWcsbar).and.(j  /=  4).and.(k  /=  4)) cycle

      tmp=xmsq

c--- The variables R1 and R2 provide the Regular and Plus pieces associated
c--- with radiation from leg 1 (R1(a,b,c,cs,is)) and leg 2 (R2(a,b,c,cs,is))
c--- In each case the parton labelling is using the normal QM notation of
c--- putting everything backward
c---       emitted line after emission =    a
c---       emitter before emission     =    b
c---       spectator                   =    c
c--- There is no label for he or she who is emitted.
c--- Note that in general each piece will be composed of many different
c--- dipole contributions

c--- special block for EW orrections to Wgamma production
      if  ((kcase==kWga_ew) .or. (kcase==kWgajew) .or. (kcase==kWln_ew)) then

      if ((j > 0) .and. (k<0)) then
      xmsq=xmsq+(msqv(j,k)
     & + msq(j,k)*(    AP(q,q,1)-AP(q,q,3)+Q1(q,q,a,1)-Q1(q,q,a,3)
     &                +AP(a,a,1)-AP(a,a,3)+Q2(a,a,q,1)-Q2(a,a,q,3)))
     &                *fx1(j)*fx2(k)
     & +(msq(j,k)*(AP(q,q,2)+AP(q,q,3)+Q1(q,q,a,2)+Q1(q,q,a,3)))*fx1z(j)/z*fx2(k)
     & +(msq(j,k)*(AP(a,a,2)+AP(a,a,3)+Q2(a,a,q,2)+Q2(a,a,q,3)))*fx1(j)*fx2z(k)/z

c      if (j==2 .and. k==-1) then
c        write(6,*) 'AP1  ',msq(j,k)*AP(q,q,1),msq(j,k)*Q1(q,q,a,1)
c        write(6,*) 'AP2  ',msq(j,k)*AP(a,a,1),msq(j,k)*Q2(a,a,q,1)
c        write(6,*) 'CT  ',msq(j,k)*(AP(q,q,1)+Q1(q,q,a,1)
c     &                             +AP(a,a,1)+Q2(a,a,q,1))
c        write(6,*) 'virt',msqv(j,k)
c        write(6,*) 'qqb sum ',msq(j,k)*(AP(q,q,1)+Q1(q,q,a,1)
c     &                                 +AP(a,a,1)+Q2(a,a,q,1))
c     &                   +msqv(j,k)
c      endif

      elseif ((j < 0) .and. (k>0)) then
c--QbarQ
      xmsq=xmsq+(msqv(j,k)
     & +msq(j,k)*(    AP(a,a,1)-AP(a,a,3)+Q1(a,a,q,1)-Q1(a,a,q,3)
     &               +AP(q,q,1)-AP(q,q,3)+Q2(q,q,a,1)-Q2(q,q,a,3)))
     &               *fx1(j)*fx2(k)
     & +(msq(j,k)*(AP(a,a,3)+AP(a,a,2)+Q1(a,a,q,3)+Q1(a,a,q,2)))*fx1z(j)/z*fx2(k)
     & +(msq(j,k)*(AP(q,q,3)+AP(q,q,2)+Q2(q,q,a,3)+Q2(q,q,a,2)))*fx1(j)*fx2z(k)/z

c      if (j == -1 .and. k == 2) then
c        write(6,*) 'AP1  ',msq(j,k)*AP(a,a,1),msq(j,k)*Q1(a,a,q,1)
c        write(6,*) 'AP2  ',msq(j,k)*AP(q,q,1),msq(j,k)*Q2(q,q,a,1)
c        write(6,*) 'CT  ',msq(j,k)*(AP(a,a,1)+Q1(a,a,q,1)
c     &                             +AP(q,q,1)+Q2(q,q,a,1))
c        write(6,*) 'virt',msqv(j,k)
c        write(6,*) 'qbq sum',msqv(j,k)
c     &   +msq(j,k)*(AP(a,a,1)+Q1(a,a,q,1)
c     &             +AP(q,q,1)+Q2(q,q,a,1))
c      endif

       endif

      if (kcase==kWgajew) then

        if     ((j > 0) .and. (k == g)) then
       xmsq=xmsq+(msqv(j,g)
     & +msq(j,g)*(+AP(q,q,1)-AP(q,q,3)+Q1(q,q,g,1)-Q1(q,q,g,3)))*fx1(j)*fx2(g)
     & +(msq(j,g)*(AP(q,q,2)+AP(q,q,3)+Q1(q,q,g,2)+Q1(q,q,g,3)))*fx1z(j)/z*fx2(g)
        elseif ((j < 0) .and. (k == g)) then
       xmsq=xmsq+(msqv(j,g)
     & +msq(j,g)*(+AP(a,a,1)-AP(a,a,3)+Q1(a,a,g,1)-Q1(a,a,g,3)))*fx1(j)*fx2(g)
     & +(msq(j,g)*(AP(a,a,2)+AP(a,a,3)+Q1(a,a,g,2)+Q1(a,a,g,3)))*fx1z(j)/z*fx2(g)
        elseif ((j == g) .and. (k > 0 )) then
       xmsq=xmsq+(msqv(g,k)
     & +msq(g,k)*(+AP(q,q,1)-AP(q,q,3)+Q2(q,q,g,1)-Q2(q,q,g,3)))*fx1(g)*fx2(k)
     & +(msq(g,k)*(AP(q,q,2)+AP(q,q,3)+Q2(q,q,g,2)+Q2(q,q,g,3)))*fx1(g)*fx2z(k)/z
        elseif ((j == g) .and. (k < 0 )) then
       xmsq=xmsq+(msqv(g,k)
     & +msq(g,k)*(+AP(a,a,1)-AP(a,a,3)+Q2(a,a,g,1)-Q2(a,a,g,3)))*fx1(g)*fx2(k)
     & +(msq(g,k)*(AP(a,a,2)+AP(a,a,3)+Q2(a,a,g,2)+Q2(a,a,g,3)))*fx1(g)*fx2z(k)/z
        endif

      endif

c--- special block for photon-induced corrections to W or Wgamma production
      elseif  ((kcase==kWln_aq) .or. (kcase==kWgaj_a)) then

       if     ((j == 0) .and. (k /= 0)) then
       msq_up=msq(+2,k)+msq(+4,k)+msq(-2,k)+msq(-4,k)
       msq_dn=msq(+1,k)+msq(+3,k)+msq(+5,k)+msq(-1,k)+msq(-3,k)+msq(-5,k)
       xmsq=xmsq+((two/three)**2*msq_up*(APqed(q,g,2)+Q1(q,g,a,2))
     &           +(one/three)**2*msq_dn*(APqed(q,g,2)+Q1(q,g,a,2)))*fx1z(g)/z*fx2(k)
       elseif ((j /= 0) .and. (k == 0)) then
       msq_up=msq(j,+2)+msq(j,+4)+msq(j,-2)+msq(j,-4)
       msq_dn=msq(j,+1)+msq(j,+3)+msq(j,+5)+msq(j,-1)+msq(j,-3)+msq(j,-5)
       xmsq=xmsq+((two/three)**2*msq_up*(APqed(q,g,2)+Q2(q,g,a,2))
     &           +(one/three)**2*msq_dn*(APqed(q,g,2)+Q2(q,g,a,2)))*fx1(j)*fx2z(g)/z
       endif

c--- special block for photon-induced corrections to Wgamma production
      elseif  (kcase==kWgajja) then

c j == 1 is a photon, k == 0 is a gluon
       if     ((j == 1) .and. (k == 0)) then
       msq_up=msq(+2,k)+msq(+4,k)+msq(-2,k)+msq(-4,k)
       msq_dn=msq(+1,k)+msq(+3,k)+msq(+5,k)+msq(-1,k)+msq(-3,k)+msq(-5,k)
       xmsq=xmsq+((two/three)**2*msq_up*(APqed(q,g,2)+Q1(q,g,a,2))
     &           +(one/three)**2*msq_dn*(APqed(q,g,2)+Q1(q,g,a,2)))*fx1z(1)/z*fx2(0)

       msq_gq=msqqcd(0,-1)+msqqcd(0,-2)+msqqcd(0,-3)+msqqcd(0,-4)+msqqcd(0,-5)
     &       +msqqcd(0,+1)+msqqcd(0,+2)+msqqcd(0,+3)+msqqcd(0,+4)+msqqcd(0,+5)

       xmsq=xmsq+msq_gq*(AP(q,g,2)+Q2(q,g,g,2))*fx1(1)*fx2z(0)/z

       elseif ((j == 0) .and. (k == 1)) then
       msq_up=msq(j,+2)+msq(j,+4)+msq(j,-2)+msq(j,-4)
       msq_dn=msq(j,+1)+msq(j,+3)+msq(j,+5)+msq(j,-1)+msq(j,-3)+msq(j,-5)
       xmsq=xmsq+((two/three)**2*msq_up*(APqed(q,g,2)+Q2(q,g,a,2))
     &           +(one/three)**2*msq_dn*(APqed(q,g,2)+Q2(q,g,a,2)))*fx1(0)*fx2z(1)/z

       msq_qg=msqqcd(-1,0)+msqqcd(-2,0)+msqqcd(-3,0)+msqqcd(-4,0)+msqqcd(-5,0)
     &       +msqqcd(+1,0)+msqqcd(+2,0)+msqqcd(+3,0)+msqqcd(+4,0)+msqqcd(+5,0)

       xmsq=xmsq+msq_qg*(AP(q,g,2)+Q1(q,g,g,2))*fx1z(0)/z*fx2(1)

       endif

c--- SUM BY COLOUR STRUCTURES: mixed QCD/EW effects for dijets only
      elseif  (kcase==ktwo_ew) then
       xmsq=xmsq+fx1(j)*fx2(k)*(msqv(j,k)+msq(j,k))
c      write(6,*) j,k,'-> msqv = ',fx1(j)*fx2(k)*(
c     & msqv(j,k)+msq(j,k))
c      tmp=xmsq

c--- quark-quark or antiquark-antiquark
      if (  ((j > 0).and.(k > 0))
     & .or. ((j < 0).and.(k < 0))) then
c      write(6,*) 'virtint: j,k,msqv=',j,k,fx1(j)*fx2(k)*msqv(j,k)
      xmsqt=zip
      do cs=0,3
      xmsqt=xmsqt
     & +msq_mix(cs,j,k)*(AP_mix(q,q,cs,1)-AP_mix(q,q,cs,3)
     &                 +M1(q,q,q,cs,1)-M1(q,q,q,cs,3)
     &                 +AP_mix(q,q,cs,1)-AP_mix(q,q,cs,3)
     &                 +M2(q,q,q,cs,1)-M2(q,q,q,cs,3))*fx1(j)*fx2(k)
     & +(msq_mix(cs,j,k)*(AP_mix(q,q,cs,2)+AP_mix(q,q,cs,3)
     &                 +M1(q,q,q,cs,2)+M1(q,q,q,cs,3))
     & +msq_mix(cs,g,k)*(AP_mix(g,q,cs,2)
     & +M1(g,q,q,cs,2)))*fx1z(j)/z*fx2(k)
     & +(msq_mix(cs,j,k)*(AP_mix(q,q,cs,2)+AP_mix(q,q,cs,3)
     &                  +M2(q,q,q,cs,2)+M2(q,q,q,cs,3))
     & +msq_mix(cs,j,g)*(AP_mix(g,q,cs,2)
     & +M2(g,q,q,cs,2)))*fx1(j)*fx2z(k)/z
      enddo
      xmsq=xmsq+xmsqt

c      write(6,*) 'virtint: j,k,  ct=',j,k,xmsqt
c      if ((j  >=  3) .and. (k  >=  3)) then
c      write(6,*) 'epinv,epinv2',epinv,epinv2
c      write(6,*) 'j,k,fx1(j)*fx2(k)*msqv(j,k),xmsqt,sum',
c     & j,k,fx1(j)*fx2(k)*msqv(j,k),xmsqt,fx1(j)*fx2(k)*msqv(j,k)+xmsqt
c      pause
c      endif


c--- quark-antiquark or antiquark-quark
      elseif (  ((j > 0).and.(k < 0))
     &     .or. ((j < 0).and.(k > 0))) then
c      write(6,*) 'virtint: j,k,msqv=',j,k,fx1(j)*fx2(k)*msqv(j,k)
      xmsqt=zip
      do cs=0,3
      if (j  /=  -k) then
c------ non-identical quarks
      xmsqt=xmsqt
     & +msq_mix(cs,j,k)*(AP_mix(q,q,cs,1)-AP_mix(q,q,cs,3)
     &                 +M1(q,q,a,cs,1)-M1(q,q,a,cs,3)
     &                 +AP_mix(a,a,cs,1)-AP_mix(a,a,cs,3)
     &                 +M2(a,a,q,cs,1)-M2(a,a,q,cs,3))*fx1(j)*fx2(k)
     & +(msq_mix(cs,j,k)*(AP_mix(q,q,cs,2)+AP_mix(q,q,cs,3)
     &                  +M1(q,q,a,cs,2)+M1(q,q,a,cs,3))
     & + msq_mix(cs,g,k)*(AP_mix(g,q,cs,2)
     & +M1(g,q,a,cs,2)))*fx1z(j)/z*fx2(k)
     & +(msq_mix(cs,j,k)*(AP_mix(a,a,cs,2)+AP_mix(a,a,cs,3)
     &                  +M2(a,a,q,cs,2)+M2(a,a,q,cs,3))
     & + msq_mix(cs,j,g)*(AP_mix(g,a,cs,2)
     & +M2(g,a,q,cs,2)))*fx1(j)*fx2z(k)/z
      else
c------ identical quarks
      xmsqt=xmsqt
     & +msq_mix(cs,j,k)*(AP_mix(q,q,cs,1)-AP_mix(q,q,cs,3)
     &                 +M1(q,q,a,cs+isame,1)-M1(q,q,a,cs+isame,3)
     &                 +AP_mix(a,a,cs,1)-AP_mix(a,a,cs,3)
     &                 +M2(a,a,q,cs+isame,1)-M2(a,a,q,cs+isame,3))
     &                  *fx1(j)*fx2(k)
     & +(msq_mix(cs,j,k)*(AP_mix(q,q,cs,2)+AP_mix(q,q,cs,3)
     &                  +M1(q,q,a,cs+isame,2)+M1(q,q,a,cs+isame,3))
     & + msq_mix(cs,g,k)*(AP_mix(g,q,cs,2)+M1(g,q,a,cs+isame,2)))
     &                  *fx1z(j)/z*fx2(k)
     & +(msq_mix(cs,j,k)*(AP_mix(a,a,cs,2)+AP_mix(a,a,cs,3)
     &                  +M2(a,a,q,cs+isame,2)+M2(a,a,q,cs+isame,3))
     & + msq_mix(cs,j,g)*(AP_mix(g,a,cs,2)+M2(g,a,q,cs+isame,2)))
     &                  *fx1(j)*fx2z(k)/z
      endif
      enddo
      xmsq=xmsq+xmsqt
c      write(6,*) 'virtint: j,k,  ct=',j,k,xmsqt
c      if (j > 0) then
c      write(6,*) 'j,k,msqv(j,k)',j,k,msqv(j,k)
c      write(6,*) 'xmsqt',xmsqt
c      pause
c      endif

      elseif (j  ==  g) then
c--gQ
       if    (k > 0) then
       xmsqt=zip
       do cs=0,3
       msq_aq=msq_mix(cs,-1,k)+msq_mix(cs,-2,k)
     &       +msq_mix(cs,-3,k)+msq_mix(cs,-4,k)
       msq_qq=msq_mix(cs,+1,k)+msq_mix(cs,+2,k)
     &       +msq_mix(cs,+3,k)+msq_mix(cs,+4,k)
       xmsqt=xmsqt
     & +(msq_aq*(AP_mix(a,g,cs,2)+M1(a,g,q,cs,2))
     &  +msq_qq*(AP_mix(q,g,cs,2)+M1(q,g,q,cs,2)))*fx1z(g)/z*fx2(k)
       enddo
       xmsq=xmsq+xmsqt

c--gQbar
       elseif (k<0) then
       xmsqt=zip
       do cs=0,3
       msq_qa=msq_mix(cs,+1,k)+msq_mix(cs,+2,k)
     &       +msq_mix(cs,+3,k)+msq_mix(cs,+4,k)
       msq_aa=msq_mix(cs,-1,k)+msq_mix(cs,-2,k)
     &       +msq_mix(cs,-3,k)+msq_mix(cs,-4,k)
       xmsqt=xmsqt
     & +(msq_qa*(AP_mix(q,g,cs,2)+M1(q,g,a,cs,2))
     &  +msq_aa*(AP_mix(a,g,cs,2)+M1(a,g,a,cs,2)))*fx1z(g)/z*fx2(k)
       enddo
       xmsq=xmsq+xmsqt

       endif
c--Qg
      elseif (k  ==  g) then
       if     (j>0) then
       xmsqt=zip
       do cs=0,3
       msq_qa=msq_mix(cs,j,-1)+msq_mix(cs,j,-2)
     &       +msq_mix(cs,j,-3)+msq_mix(cs,j,-4)
       msq_qq=msq_mix(cs,j,+1)+msq_mix(cs,j,+2)
     &       +msq_mix(cs,j,+3)+msq_mix(cs,j,+4)
       xmsqt=xmsqt
     & +(msq_qa*(AP_mix(a,g,cs,2)+M2(a,g,q,cs,2))
     &  +msq_qq*(AP_mix(q,g,cs,2)+M2(q,g,q,cs,2)))*fx1(j)*fx2z(g)/z
       enddo
       xmsq=xmsq+xmsqt

c--Qbarg
       elseif (j<0) then
       xmsqt=zip
       do cs=0,3
       msq_aq=msq_mix(cs,j,+1)+msq_mix(cs,j,+2)
     &       +msq_mix(cs,j,+3)+msq_mix(cs,j,+4)
       msq_aa=msq_mix(cs,j,-1)+msq_mix(cs,j,-2)
     &       +msq_mix(cs,j,-3)+msq_mix(cs,j,-4)
       xmsqt=xmsqt
     & +(msq_aq*(AP_mix(q,g,cs,2)+M2(q,g,a,cs,2))
     &  +msq_aa*(AP_mix(a,g,cs,2)+M2(a,g,a,cs,2)))*fx1(j)*fx2z(g)/z
       enddo
       xmsq=xmsq+xmsqt

       endif

      endif

c--- END ktwo_ew
c--- SUM BY COLOUR STRUCTURES: H+2jets only
      elseif  ( (kcase==kggfus2) .or. (kcase==kgagajj)
     &     .or. (kcase==kHWW2jt) .or. (kcase==kHZZ2jt)) then
       xmsq=xmsq+fx1(j)*fx2(k)*(
     & msqv(j,k)+msq(j,k))
c      write(6,*) j,k,'-> msqv = ',fx1(j)*fx2(k)*(
c     & msqv(j,k)+msq_cs(0,j,k)+msq_cs(1,j,k)+msq_cs(2,j,k))
c      tmp=xmsq

c--- quark-quark or antiquark-antiquark
      if (  ((j > 0).and.(k > 0))
     & .or. ((j < 0).and.(k < 0))) then
c      write(6,*) 'virtint: j,k,msqv=',j,k,fx1(j)*fx2(k)*msqv(j,k)
      if (j == k) then
        m=+1
        n=+1
        csmax=6
      else
        m=abs(j)
        n=abs(k)
        csmax=6
      endif
      xmsqt=0._dp
      do cs=1,csmax
      xmsqt=xmsqt
     & +msq_struc(cs,m,n)*(AP(q,q,1)-AP(q,q,3)
     &                 +H1(q,q,q,cs,1)-H1(q,q,q,cs,3)
     &                 +AP(q,q,1)-AP(q,q,3)
     &                 +H2(q,q,q,cs,1)-H2(q,q,q,cs,3))*fx1(j)*fx2(k)
     & +(msq_struc(cs,m,n)*(AP(q,q,2)+AP(q,q,3)
     &                  +H1(q,q,q,cs,2)+H1(q,q,q,cs,3))
     & +msq_struc(cs,g,k)*(AP(g,q,2)+H1(g,q,q,cs,2)))*fx1z(j)/z*fx2(k)
     & +(msq_struc(cs,m,n)*(AP(q,q,2)+AP(q,q,3)
     &                  +H2(q,q,q,cs,2)+H2(q,q,q,cs,3))
     & +msq_struc(cs,j,g)*(AP(g,q,2)+H2(g,q,q,cs,2)))*fx1(j)*fx2z(k)/z
      enddo
      xmsq=xmsq+xmsqt
c      write(6,*) 'virtint: j,k,  ct=',j,k,xmsqt

c--- quark-antiquark or antiquark-quark
      elseif (  ((j > 0).and.(k < 0))
     &     .or. ((j < 0).and.(k > 0))) then
c      write(6,*) 'virtint: j,k,msqv=',j,k,fx1(j)*fx2(k)*msqv(j,k)
      if (j == -k) then
        m=+1
        n=-1
        csmax=7
      else
        m=abs(j)
        n=-abs(k)
        csmax=6
      endif
      xmsqt=0._dp
      do cs=1,csmax
c      if ((cs > 3) .and. (cs < 7)) goto 67
c      do cs=4,6
      xmsqt=xmsqt
     & +msq_struc(cs,m,n)*(AP(q,q,1)-AP(q,q,3)
     &                 +H1(q,q,a,cs,1)-H1(q,q,a,cs,3)
     &                 +AP(a,a,1)-AP(a,a,3)
     &                 +H2(a,a,q,cs,1)-H2(a,a,q,cs,3))*fx1(j)*fx2(k)
     & +(msq_struc(cs,m,n)*(AP(q,q,2)+AP(q,q,3)
     &                  +H1(q,q,a,cs,2)+H1(q,q,a,cs,3))
     & + msq_struc(cs,g,k)*(AP(g,q,2)+H1(g,q,a,cs,2)))*fx1z(j)/z*fx2(k)
     & +(msq_struc(cs,m,n)*(AP(a,a,2)+AP(a,a,3)
     &                  +H2(a,a,q,cs,2)+H2(a,a,q,cs,3))
     & + msq_struc(cs,j,g)*(AP(g,a,2)+H2(g,a,q,cs,2)))*fx1(j)*fx2z(k)/z
c   67 continue
      enddo
      xmsq=xmsq+xmsqt
c      write(6,*) 'virtint: j,k,  ct=',j,k,xmsqt

c--- gluon-gluon
      elseif ((j == g) .and. (k == g)) then
c      write(6,*) 'virtint: j,k,msqv=',j,k,fx1(j)*fx2(k)*msqv(j,k)
      xmsqt=0._dp
c--- loop up to 3 to cancel poles from gggg
c      do cs=1,3
      do cs=1,6
      msq_qg=msq_struc(cs,+5,g)+msq_struc(cs,+4,g)+msq_struc(cs,+3,g)
     &      +msq_struc(cs,+2,g)+msq_struc(cs,+1,g)
     &      +msq_struc(cs,-5,g)+msq_struc(cs,-4,g)+msq_struc(cs,-3,g)
     &      +msq_struc(cs,-2,g)+msq_struc(cs,-1,g)
      msq_gq=msq_struc(cs,g,+5)+msq_struc(cs,g,+4)+msq_struc(cs,g,+3)
     &      +msq_struc(cs,g,+2)+msq_struc(cs,g,+1)
     &      +msq_struc(cs,g,-5)+msq_struc(cs,g,-4)+msq_struc(cs,g,-3)
     &      +msq_struc(cs,g,-2)+msq_struc(cs,g,-1)
      xmsqt=xmsqt
     & +msq_struc(cs,g,g)*(AP(g,g,1)-AP(g,g,3)
     &                 +H1(g,g,g,cs,1)-H1(g,g,g,cs,3)
     &                 +AP(g,g,1)-AP(g,g,3)
     &                 +H2(g,g,g,cs,1)-H2(g,g,g,cs,3))*fx1(g)*fx2(g)
     & +(msq_struc(cs,g,g)*(AP(g,g,2)+AP(g,g,3)
     &                  +H1(g,g,g,cs,2)+H1(g,g,g,cs,3))
     &  +msq_qg*(AP(q,g,2)+H1(q,g,g,cs,2)))*fx1z(g)/z*fx2(g)
     & +(msq_struc(cs,g,g)*(AP(g,g,2)+AP(g,g,3)
     &                  +H2(g,g,g,cs,2)+H2(g,g,g,cs,3))
     &  +msq_gq*(AP(q,g,2)+H2(q,g,g,cs,2)))*fx1(g)*fx2z(g)/z
      enddo
      xmsq=xmsq+xmsqt
c      write(6,*) 'virtint: j,k,  ct=',j,k,xmsqt

c--- quark-gluon and anti-quark gluon
      elseif ((j  /=  0) .and. (k == g)) then
c      write(6,*) 'virtint: j,k,msqv=',j,k,fx1(j)*fx2(k)*msqv(j,k)
      m=+1
      n=0
      xmsqt=0._dp
      do cs=4,6
      xmsqt=xmsqt
     &+ msq_struc(cs,m,g)*(AP(q,q,1)-AP(q,q,3)
     &                 +H1(q,q,g,cs,1)-H1(q,q,g,cs,3)
     &                 +H2(g,g,q,cs,1)-H2(g,g,q,cs,3)
     &                 +AP(g,g,1)-AP(g,g,3))*fx1(j)*fx2(g)
     &+(msq_struc(cs,m,g)*(AP(q,q,2)+AP(q,q,3)
     &                    +H1(q,q,g,cs,2)+H1(q,q,g,cs,3))
     & +msq_struc(cs,g,g)*(AP(g,q,2)+H1(g,q,g,cs,2)))*fx1z(j)/z*fx2(g)
     &+(msq_struc(cs,m,g)*(AP(g,g,2)+AP(g,g,3)
     &                    +H2(g,g,q,cs,2)+H2(g,g,q,cs,3))
     & +msq_struc(cs,m,-m)*(AP(a,g,2)+H2(a,g,q,cs,2)))*fx1(j)*fx2z(g)/z
      enddo
      do cs=1,3
      msq_qa=msq_struc(cs,m,-1)+msq_struc(cs,m,-2)+msq_struc(cs,m,-3)
     &      +msq_struc(cs,m,-4)+msq_struc(cs,m,-5)
      msq_qq=msq_struc(cs,m,+1)+msq_struc(cs,m,+2)+msq_struc(cs,m,+3)
     &      +msq_struc(cs,m,+4)+msq_struc(cs,m,+5)
      xmsqt=xmsqt
     &+(msq_struc(cs,g,g)*(AP(g,q,2)+H1(g,q,g,cs,2)))*fx1z(j)/z*fx2(g)
     &+(+msq_qa*(AP(a,g,2)+H2(a,g,q,cs,2))
     &  +msq_qq*(AP(q,g,2)+H2(q,g,q,cs,2)))*fx1(j)*fx2z(g)/z
      enddo
      xmsqt=xmsqt
     & +msq_struc(iqr,m,-m)*(AP(a,g,2)+H2(a,g,q,iqr,2))*fx1(j)*fx2z(g)/z
      xmsq=xmsq+xmsqt
c      write(6,*) 'virtint: j,k,  ct=',j,k,xmsqt
c      write(6,*) 'virtint: j,k, SUM=',j,k,xmsqt+fx1(j)*fx2(k)*msqv(j,k)

c--- gluon-quark and gluon anti-quark
      elseif ((j == 0) .and. (k  /=  0)) then
c      write(6,*) 'virtint: j,k,msqv=',j,k,fx1(j)*fx2(k)*msqv(j,k)
      m=0
      n=+1
      xmsqt=0._dp
      do cs=4,6
      xmsqt=xmsqt
     & +msq_struc(cs,g,n)*(AP(g,g,1)-AP(g,g,3)
     &                 +H1(g,g,q,cs,1)-H1(g,g,q,cs,3)
     &                 +AP(q,q,1)-AP(q,q,3)
     &                 +H2(q,q,g,cs,1)-H2(q,q,g,cs,3))*fx1(g)*fx2(k)
     & +(msq_struc(cs,g,n)*(AP(g,g,2)+AP(g,g,3)
     &                 +H1(g,g,q,cs,2)+H1(g,g,q,cs,3))
     &  +msq_struc(cs,-n,n)*(AP(a,g,2)+H1(a,g,q,cs,2)))*fx1z(g)/z*fx2(k)
     & +(msq_struc(cs,g,n)*(AP(q,q,2)+AP(q,q,3)
     &                 +H2(q,q,g,cs,2)+H2(q,q,g,cs,3))
     &  +msq_struc(cs,g,g)*(AP(g,q,2)+H2(g,q,g,cs,2)))*fx1(g)*fx2z(k)/z
      enddo
      do cs=1,3
      msq_aq=msq_struc(cs,-1,n)+msq_struc(cs,-2,n)+msq_struc(cs,-3,n)
     &      +msq_struc(cs,-4,n)+msq_struc(cs,-5,n)
      msq_qq=msq_struc(cs,+1,n)+msq_struc(cs,+2,n)+msq_struc(cs,+3,n)
     &      +msq_struc(cs,+4,n)+msq_struc(cs,+5,n)
      xmsqt=xmsqt
     & +(+msq_aq*(AP(a,g,2)+H1(a,g,q,cs,2))
     &   +msq_qq*(AP(q,g,2)+H1(q,g,q,cs,2)))*fx1z(g)/z*fx2(k)
     & +(+msq_struc(cs,g,g)*(AP(g,q,2)+H2(g,q,g,cs,2)))*fx1(g)*fx2z(k)/z
      enddo
      xmsqt=xmsqt
     & +msq_struc(iqr,n,-n)*(AP(a,g,2)+H1(a,g,q,iqr,2))*fx1z(g)/z*fx2(k)
      xmsq=xmsq+xmsqt
c      write(6,*) 'virtint: j,k,  ct=',j,k,xmsqt

      endif

      elseif (kcase==ktwojet)then

      xmsq=xmsq+fx1(j)*fx2(k)*(
     & msqv(j,k)+msq_cs(0,j,k)+msq_cs(1,j,k)+msq_cs(2,j,k))

c--- END COLORLESS + 2JET
c--- SUM BY COLOUR STRUCTURES AND FINAL STATES: 2 jets only
      if      ((j == 0) .and. (k == 0)) then
      tmp=xmsq
      do cs=0,2
      msq_qg=real(nf,dp)*msqx(cs,+1,g,+1,g)
      msq_gq=real(nf,dp)*msqx(cs,g,+1,+1,g)
      xmsq=xmsq
     & +msqx(cs,g,g,q,a)*real(nf,dp)*(AP(g,g,1)-AP(g,g,3)
     &     +S1(g,g,g,qf_af,cs,1)-S1(g,g,g,qf_af,cs,3)
     &     +AP(g,g,1)-AP(g,g,3)
     &     +S2(g,g,g,qf_af,cs,1)-S2(g,g,g,qf_af,cs,3))*fx1(g)*fx2(g)
     & +msqx(cs,g,g,q,a)*real(nf,dp)*(AP(g,g,2)+AP(g,g,3)
     &     +S1(g,g,g,qf_af,cs,3)+S1(g,g,g,qf_af,cs,2))*fx1z(g)/z*fx2(g)
     & +msqx(cs,g,g,q,a)*real(nf,dp)*(AP(g,g,2)+AP(g,g,3)
     &     +S2(g,g,g,qf_af,cs,3)+S2(g,g,g,qf_af,cs,2))*fx1(g)*fx2z(g)/z
     & +msq_qg*(AP(q,g,2)+S1(q,g,g,qf_gf,cs,2)
     &         +AP(a,g,2)+S1(a,g,g,af_gf,cs,2))*fx1z(g)/z*fx2(g)
     & +msq_gq*(AP(q,g,2)+S2(q,g,g,qf_gf,cs,2)
     &         +AP(a,g,2)+S2(a,g,g,af_gf,cs,2))*fx1(g)*fx2z(g)/z
c     & +msqx(cs,g,g,g,g)*(AP(g,g,1)-AP(g,g,3)
c     &     +S1(g,g,g,gf_gf,cs,1)-S1(g,g,g,gf_gf,cs,3)
c     &     +AP(g,g,1)-AP(g,g,3)
c     &     +S2(g,g,g,gf_gf,cs,1)-S2(g,g,g,gf_gf,cs,3))*fx1(g)*fx2(g)
c     & +msqx(cs,g,g,g,g)*(AP(g,g,2)+AP(g,g,3)
c     &     +S1(g,g,g,gf_gf,cs,3)+S1(g,g,g,gf_gf,cs,2))*fx1z(g)/z*fx2(g)
c     & +msqx(cs,g,g,g,g)*(AP(g,g,2)+AP(g,g,3)
c     &     +S2(g,g,g,gf_gf,cs,3)+S2(g,g,g,gf_gf,cs,2))*fx1(g)*fx2z(g)/z
      enddo
      write(6,*) '_v: ',tmp
      write(6,*) '_z: ',xmsq-tmp
      endif

      elseif ((kcase==kW_2jet) .or. (kcase==kZ_2jet)
     &   .or. (kcase==kW_bjet) .or. (kcase==kZ_bjet)
     &   .or. (kcase==ktt_bbl) .or. (kcase==ktt_bbh)
     &   .or. (kcase==ktt_bbu)
     &   .or. (kcase==ktt_tot) .or. (kcase==kbb_tot)
     &   .or. (kcase==kcc_tot)) then
c--- SUM BY COLOUR STRUCTURES: W/Z + 2 jet only

      xmsq=xmsq+fx1(j)*fx2(k)*(
     & msqv(j,k)+msq_cs(0,j,k)+msq_cs(1,j,k)+msq_cs(2,j,k))
c      write(6,*) j,k,'-> msqv = ',fx1(j)*fx2(k)*(
c     & msqv(j,k)+msq_cs(0,j,k)+msq_cs(1,j,k)+msq_cs(2,j,k))
c      tmp=xmsq


      if ((j > 0) .and. (k>0)) then
      do cs=0,2
c--- handle the Z+b jet case where identity of 5 and 6 are switched
c---   (interchange colour structures 1 and 2)
      ics=cs
      if ( (kcase==kZ_bjet) .and. (k == +flav)
     &     .and. (cs > 0) ) ics=3-cs
      xmsq=xmsq
     & +msq_cs(ics,j,k)*(AP(q,q,1)-AP(q,q,3)
     &                  +R1(q,q,q,cs,1)-R1(q,q,q,cs,3)
     &                  +AP(q,q,1)-AP(q,q,3)
     &                  +R2(q,q,q,cs,1)-R2(q,q,q,cs,3))*fx1(j)*fx2(k)
     & +(msq_cs(ics,j,k)*(AP(q,q,2)+AP(q,q,3)
     &                   +R1(q,q,q,cs,2)+R1(q,q,q,cs,3))
     & + msq_cs(ics,g,k)*(AP(g,q,2)+R1(g,q,q,cs,2)))*fx1z(j)/z*fx2(k)
     & +(msq_cs(ics,j,k)*(AP(q,q,2)+AP(q,q,3)
     &                   +R2(q,q,q,cs,2)+R2(q,q,q,cs,3))
     & + msq_cs(ics,j,g)*(AP(g,q,2)+R2(g,q,q,cs,2)))*fx1(j)*fx2z(k)/z
      enddo
      elseif ((j < 0) .and. (k<0)) then
      do cs=0,2
c--- handle the Z+b jet case where identity of 5 and 6 are switched
c---   (interchange colour structures 1 and 2)
      ics=cs
      if ( (kcase==kZ_bjet) .and. (k == -flav)
     &     .and. (cs > 0) ) ics=3-cs
      xmsq=xmsq
     & +msq_cs(ics,j,k)*(AP(a,a,1)-AP(a,a,3)
     &                 +R1(a,a,a,cs,1)-R1(a,a,a,cs,3)
     &                 +AP(a,a,1)-AP(a,a,3)
     &                 +R2(a,a,a,cs,1)-R2(a,a,a,cs,3))*fx1(j)*fx2(k)
     & +(msq_cs(ics,j,k)*(AP(a,a,2)+AP(a,a,3)
     &                  +R1(a,a,a,cs,2)+R1(a,a,a,cs,3))
     & + msq_cs(ics,g,k)*(AP(g,a,2)+R1(g,a,a,cs,2)))*fx1z(j)/z*fx2(k)
     & +(msq_cs(ics,j,k)*(AP(a,a,2)+AP(a,a,3)
     &                  +R2(a,a,a,cs,2)+R2(a,a,a,cs,3))
     & + msq_cs(ics,j,g)*(AP(g,a,2)+R2(g,a,a,cs,2)))*fx1(j)*fx2z(k)/z

      enddo
      elseif ((j > 0) .and. (k<0)) then
c--- handle the Z+b jet case where identity of 5 and 6 are switched
c---   (interchange integrated CT's for q-qbar and qbar-q)
      iq=q
      ia=a
      if ((kcase==kZ_bjet) .and. (k == -flav)) then
        iq=a
        ia=q
      endif
      do cs=0,2
      xmsq=xmsq
     & +msq_cs(cs,j,k)*(AP(q,q,1)-AP(q,q,3)
     &                 +R1(iq,iq,ia,cs,1)-R1(iq,iq,ia,cs,3)
     &                 +AP(a,a,1)-AP(a,a,3)
     &                 +R2(ia,ia,iq,cs,1)-R2(ia,ia,iq,cs,3)
     &                   )*fx1(j)*fx2(k)
     & +(msq_cs(cs,j,k)*(AP(q,q,2)+AP(q,q,3)
     &                  +R1(iq,iq,ia,cs,2)+R1(iq,iq,ia,cs,3))
     & + msq_cs(cs,g,k)*(AP(g,q,2)+R1(g,q,a,cs,2)))*fx1z(j)/z*fx2(k)
     & +(msq_cs(cs,j,k)*(AP(a,a,2)+AP(a,a,3)
     &                  +R2(ia,ia,iq,cs,2)+R2(ia,ia,iq,cs,3))
     & + msq_cs(cs,j,g)*(AP(g,a,2)+R2(g,a,q,cs,2)))*fx1(j)*fx2z(k)/z

      enddo
      elseif ((j < 0) .and. (k>0)) then
c--- handle the Z+b jet case where identity of 5 and 6 are switched
c---   (interchange integrated CT's for q-qbar and qbar-q)
      iq=q
      ia=a
      if ((kcase==kZ_bjet) .and. (j == -flav)) then
        iq=a
        ia=q
      endif
      do cs=0,2
      xmsq=xmsq
     & +msq_cs(cs,j,k)*(AP(a,a,1)-AP(a,a,3)
     &                 +R1(ia,ia,iq,cs,1)-R1(ia,ia,iq,cs,3)
     &                 +AP(q,q,1)-AP(q,q,3)
     &                 +R2(iq,iq,ia,cs,1)-R2(iq,iq,ia,cs,3)
     &                   )*fx1(j)*fx2(k)
     & +(msq_cs(cs,j,k)*(AP(a,a,2)+AP(a,a,3)
     &                  +R1(ia,ia,iq,cs,2)+R1(ia,ia,iq,cs,3))
     & + msq_cs(cs,g,k)*(AP(g,a,2)+R1(g,a,q,cs,2)))*fx1z(j)/z*fx2(k)
     & +(msq_cs(cs,j,k)*(AP(q,q,2)+AP(q,q,3)
     &                  +R2(iq,iq,ia,cs,2)+R2(iq,iq,ia,cs,3))
     & + msq_cs(cs,j,g)*(AP(g,q,2)+R2(g,q,a,cs,2)))*fx1(j)*fx2z(k)/z

      enddo
      elseif ((j == g) .and. (k == g)) then
      do cs=0,2
      msq_qg=msq_cs(cs,+5,g)+msq_cs(cs,+4,g)+msq_cs(cs,+3,g)
     &      +msq_cs(cs,+2,g)+msq_cs(cs,+1,g)
     &      +msq_cs(cs,-5,g)+msq_cs(cs,-4,g)+msq_cs(cs,-3,g)
     &      +msq_cs(cs,-2,g)+msq_cs(cs,-1,g)
      msq_gq=msq_cs(cs,g,+5)+msq_cs(cs,g,+4)+msq_cs(cs,g,+3)
     &      +msq_cs(cs,g,+2)+msq_cs(cs,g,+1)
     &      +msq_cs(cs,g,-5)+msq_cs(cs,g,-4)+msq_cs(cs,g,-3)
     &      +msq_cs(cs,g,-2)+msq_cs(cs,g,-1)
      xmsq=xmsq
     & +msq_cs(cs,g,g)*(AP(g,g,1)-AP(g,g,3)
     &                 +R1(g,g,g,cs,1)-R1(g,g,g,cs,3)
     &                 +AP(g,g,1)-AP(g,g,3)
     &                 +R2(g,g,g,cs,1)-R2(g,g,g,cs,3))*fx1(g)*fx2(g)
     & +(msq_cs(cs,g,g)*(AP(g,g,2)+AP(g,g,3)
     &                  +R1(g,g,g,cs,3)+R1(g,g,g,cs,2))
     & + msq_qg*(AP(q,g,2)+R1(q,g,g,cs,2)))*fx1z(g)/z*fx2(g)
     & +(msq_cs(cs,g,g)*(AP(g,g,2)+AP(g,g,3)
     &                  +R2(g,g,g,cs,3)+R2(g,g,g,cs,2))
     & + msq_gq*(AP(q,g,2)+R2(q,g,g,cs,2)))*fx1(g)*fx2z(g)/z
      enddo
      elseif ((j == g) .and. (k > 0)) then
c--- special case for W+bj - remove b-PDF contribution
      if ((kcase==kW_bjet) .and. (k  /=  5)) then
      xmsq=xmsq+(msq(-5,k)+msq(+5,k))*APqg_mass*fx1z(g)/z*fx2(k)
      else
      do cs=0,2
      msq_aq=msq_cs(cs,-1,k)+msq_cs(cs,-2,k)+msq_cs(cs,-3,k)
     &      +msq_cs(cs,-4,k)+msq_cs(cs,-5,k)
      msq_qq=msq_cs(cs,+1,k)+msq_cs(cs,+2,k)+msq_cs(cs,+3,k)
     &      +msq_cs(cs,+4,k)+msq_cs(cs,+5,k)
      xmsq=xmsq
     & +msq_cs(cs,g,k)*(AP(g,g,1)-AP(g,g,3)
     &                 +R1(g,g,q,cs,1)-R1(g,g,q,cs,3)
     &                 +AP(q,q,1)-AP(q,q,3)
     &                 +R2(q,q,g,cs,1)-R2(q,q,g,cs,3))*fx1(g)*fx2(k)
     & +(msq_cs(cs,g,k)*(AP(g,g,2)+AP(g,g,3)
     &                  +R1(g,g,q,cs,2)+R1(g,g,q,cs,3))
     & + msq_aq*(AP(a,g,2)+R1(a,g,q,cs,2))
     & + msq_qq*(AP(q,g,2)+R1(q,g,q,cs,2)))*fx1z(g)/z*fx2(k)
     & +(msq_cs(cs,g,k)*(AP(q,q,2)+AP(q,q,3)
     &                +R2(q,q,g,cs,2)+R2(q,q,g,cs,3))
     & + msq_cs(cs,g,g)*(AP(g,q,2)+R2(g,q,g,cs,2)))*fx1(g)*fx2z(k)/z

      enddo
      endif

      elseif ((j == g) .and. (k < 0)) then
c--- special case for W+bj - remove b-PDF contribution
      if ((kcase==kW_bjet) .and. (k  /=  -5)) then
      xmsq=xmsq+(msq(-5,k)+msq(+5,k))*APqg_mass*fx1z(g)/z*fx2(k)
      else
      do cs=0,2
      msq_qa=msq_cs(cs,+1,k)+msq_cs(cs,+2,k)+msq_cs(cs,+3,k)
     &      +msq_cs(cs,+4,k)+msq_cs(cs,+5,k)
      msq_aa=msq_cs(cs,-1,k)+msq_cs(cs,-2,k)+msq_cs(cs,-3,k)
     &      +msq_cs(cs,-4,k)+msq_cs(cs,-5,k)
      xmsq=xmsq
     & +msq_cs(cs,g,k)*(AP(g,g,1)-AP(g,g,3)
     &                 +R1(g,g,a,cs,1)-R1(g,g,a,cs,3)
     &                 +AP(a,a,1)-AP(a,a,3)
     &                 +R2(a,a,g,cs,1)-R2(a,a,g,cs,3))*fx1(g)*fx2(k)
     & +(msq_cs(cs,g,k)*(AP(g,g,2)+AP(g,g,3)
     &                  +R1(g,g,a,cs,2)+R1(g,g,a,cs,3))
     & + msq_qa*(AP(q,g,2)+R1(q,g,a,cs,2))
     & + msq_aa*(AP(a,g,2)+R1(a,g,a,cs,2)))*fx1z(g)/z*fx2(k)
     & +(msq_cs(cs,g,k)*(AP(a,a,2)+AP(a,a,3)
     &                  +R2(a,a,g,cs,2)+R2(a,a,g,cs,3))
     & + msq_cs(cs,g,g)*(AP(g,a,2)+R2(g,a,g,cs,2)))*fx1(g)*fx2z(k)/z

      enddo
      endif

      elseif ((j > 0) .and. (k == g)) then
c--- special case for W+bj - remove b-PDF contribution
      if ((kcase==kW_bjet) .and. (j  /=  5)) then
      xmsq=xmsq+(msq(j,-5)+msq(j,+5))*APqg_mass*fx1(j)*fx2z(g)/z
      else
      do cs=0,2
      msq_qa=msq_cs(cs,j,-1)+msq_cs(cs,j,-2)+msq_cs(cs,j,-3)
     &      +msq_cs(cs,j,-4)+msq_cs(cs,j,-5)
      msq_qq=msq_cs(cs,j,+1)+msq_cs(cs,j,+2)+msq_cs(cs,j,+3)
     &      +msq_cs(cs,j,+4)+msq_cs(cs,j,+5)
       xmsq=xmsq
     &+ msq_cs(cs,j,g)*(AP(q,q,1)-AP(q,q,3)
     &                 +R1(q,q,g,cs,1)-R1(q,q,g,cs,3)
     &                 +R2(g,g,q,cs,1)-R2(g,g,q,cs,3)
     &                 +AP(g,g,1)-AP(g,g,3))*fx1(j)*fx2(g)
     &+(msq_cs(cs,j,g)*(AP(q,q,2)+AP(q,q,3)
     &                 +R1(q,q,g,cs,2)+R1(q,q,g,cs,3))
     &+ msq_cs(cs,g,g)*(AP(g,q,2)+R1(g,q,g,cs,2)))*fx1z(j)/z*fx2(g)
     &+(msq_cs(cs,j,g)*(AP(g,g,2)+AP(g,g,3)
     &                 +R2(g,g,q,cs,2)+R2(g,g,q,cs,3))
     &+ msq_qa*(AP(a,g,2)+R2(a,g,q,cs,2))
     &+ msq_qq*(AP(q,g,2)+R2(q,g,q,cs,2)))*fx1(j)*fx2z(g)/z

      enddo
      endif

      elseif ((j < 0) .and. (k == g)) then
c--- special case for W+bj - remove b-PDF contribution
      if ((kcase==kW_bjet) .and. (j  /=  -5)) then
      xmsq=xmsq+(msq(j,-5)+msq(j,+5))*APqg_mass*fx1(j)*fx2z(g)/z
      else
      do cs=0,2
      msq_aq=msq_cs(cs,j,+1)+msq_cs(cs,j,+2)+msq_cs(cs,j,+3)
     &      +msq_cs(cs,j,+4)+msq_cs(cs,j,+5)
      msq_aa=msq_cs(cs,j,-1)+msq_cs(cs,j,-2)+msq_cs(cs,j,-3)
     &      +msq_cs(cs,j,-4)+msq_cs(cs,j,-5)
       xmsq=xmsq
     & + msq_cs(cs,j,g)*(AP(a,a,1)-AP(a,a,3)
     &                  +R1(a,a,g,cs,1)-R1(a,a,g,cs,3)
     &                  +AP(g,g,1)-AP(g,g,3)
     &                  +R2(g,g,a,cs,1)-R2(g,g,a,cs,3))*fx1(j)*fx2(g)
     & +(msq_cs(cs,j,g)*(AP(a,a,2)+AP(a,a,3)
     &                  +R1(a,a,g,cs,2)+R1(a,a,g,cs,3))
     & + msq_cs(cs,g,g)*(AP(g,a,2)+R1(g,a,g,cs,2)))*fx1z(j)/z*fx2(g)
     & +(msq_cs(cs,j,g)*(AP(g,g,2)+AP(g,g,3)
     &                  +R2(g,g,a,cs,2)+R2(g,g,a,cs,3))
     & + msq_aq*(AP(q,g,2)+R2(q,g,a,cs,2))
     & + msq_aa*(AP(a,g,2)+R2(a,g,a,cs,2)))*fx1(j)*fx2z(g)/z

      enddo
      endif
      endif

c--- END 2 JET

c--- FOUR FLAVOR SINGLE TOP
      elseif ((kcase==kqg_tbq) .or. (kcase==k4ftwdk)
     &    .or.(kcase==kdk_4ft)) then

c--- SPECIAL SUM FOR SINGLE TOP + B CASE
c--QQ
      if     ((j > 0) .and. (k>0)) then
      xmsq=xmsq
     & +(msq(g,k)*(AP(g,q,2)+Q1(g,q,q,2)))*fx1z_H(j)/z*fx2_L(k)
     & +(msq(j,g)*(AP(g,q,2)+Q2(g,q,q,2)))*fx1_L(j)*fx2z_H(k)/z
c--QbarQbar
      elseif ((j < 0) .and. (k<0)) then
      xmsq=xmsq
     & +(msq(g,k)*(AP(g,a,2)+Q1(g,a,a,2)))*fx1z_H(j)/z*fx2_L(k)
     & +(msq(j,g)*(AP(g,a,2)+Q2(g,a,a,2)))*fx1_L(j)*fx2z_H(k)/z
c--QQbar
      elseif ((j > 0) .and. (k<0)) then
      xmsq=xmsq
     & +(msq(g,k)*(AP(g,q,2)+Q1(g,q,a,2)))*fx1z_H(j)/z*fx2_L(k)
     & +(msq(j,g)*(AP(g,a,2)+Q2(g,a,q,2)))*fx1_L(j)*fx2z_H(k)/z

      elseif ((j < 0) .and. (k>0)) then
c--QbarQ
      xmsq=xmsq
     & +(msq(g,k)*(AP(g,a,2)+Q1(g,a,q,2)))*fx1z_H(j)/z*fx2_L(k)
     & +(msq(j,g)*(AP(g,q,2)+Q2(g,q,a,2)))*fx1_L(j)*fx2z_H(k)/z

      elseif ((j == g) .and. (k==g)) then
c--gg
       msq_qg=msq(+5,g)+msq(+4,g)+msq(+3,g)+msq(+2,g)+msq(+1,g)
     &       +msq(-5,g)+msq(-4,g)+msq(-3,g)+msq(-2,g)+msq(-1,g)
       msq_gq=msq(g,+5)+msq(g,+4)+msq(g,+3)+msq(g,+2)+msq(g,+1)
     &       +msq(g,-5)+msq(g,-4)+msq(g,-3)+msq(g,-2)+msq(g,-1)
       xmsq=xmsq
     &  +(msq_qg*(AP(q,g,2)+Q1(q,g,g,2)))*fx1z_L(g)/z*fx2_H(g)
     &  +(msq_gq*(AP(q,g,2)+Q2(q,g,g,2)))*fx1_H(g)*fx2z_L(g)/z

      elseif (j == g) then
c--gQ
       if    (k > 0) then
       msq_aq=msq(-1,k)+msq(-2,k)+msq(-3,k)+msq(-4,k)+msq(-5,k)
       msq_qq=msq(+1,k)+msq(+2,k)+msq(+3,k)+msq(+4,k)+msq(+5,k)
       xmsq=xmsq+(msqv(g,k)
     & +msq(g,k)*(one+AP(g,g,1)-AP(g,g,3)+Q1(g,g,q,1)-Q1(g,g,q,3)
     &               +AP(q,q,1)-AP(q,q,3)+Q2(q,q,g,1)-Q2(q,q,g,3)))
     &               *fx1_H(g)*fx2_L(k)
     & +(msq(g,k)*(AP(g,g,2)+AP(g,g,3)+Q1(g,g,q,2)+Q1(g,g,q,3))
     & +   msq_aq*(AP(a,g,2)+Q1(a,g,q,2))
     & +   msq_qq*(AP(q,g,2)+Q1(q,g,q,2)))*fx1z_H(g)/z*fx2_L(k)
     & +(msq(g,k)*(AP(q,q,2)+AP(q,q,3)+Q2(q,q,g,2)+Q2(q,q,g,3))
     & + msq(g,g)*(AP(g,q,2)+Q2(g,q,g,2)))*fx1_H(g)*fx2z_L(k)/z
c--gQbar
       elseif (k<0) then
       msq_qa=msq(+1,k)+msq(+2,k)+msq(+3,k)+msq(+4,k)+msq(+5,k)
       msq_aa=msq(-1,k)+msq(-2,k)+msq(-3,k)+msq(-4,k)+msq(-5,k)
       xmsq=xmsq+(msqv(g,k)
     & +msq(g,k)*(one+AP(g,g,1)-AP(g,g,3)+Q1(g,g,a,1)-Q1(g,g,a,3)
     &               +AP(a,a,1)-AP(a,a,3)+Q2(a,a,g,1)-Q2(a,a,g,3)))
     &               *fx1_H(g)*fx2_L(k)
     & +(msq(g,k)*(AP(g,g,2)+AP(g,g,3)+Q1(g,g,a,2)+Q1(g,g,a,3))
     & +   msq_qa*(AP(q,g,2)+Q1(q,g,a,2))
     & +   msq_aa*(AP(a,g,2)+Q1(a,g,a,2)))*fx1z_H(g)/z*fx2_L(k)
     & +(msq(g,k)*(AP(a,a,2)+AP(a,a,3)+Q2(a,a,g,2)+Q2(a,a,g,3))
     & + msq(g,g)*(AP(g,a,2)+Q2(g,a,g,2)))*fx1_H(g)*fx2z_L(k)/z
       endif
c--Qg
      elseif (k == g) then
       if     (j>0) then
       msq_qa=msq(j,-1)+msq(j,-2)+msq(j,-3)+msq(j,-4)+msq(j,-5)
       msq_qq=msq(j,+1)+msq(j,+2)+msq(j,+3)+msq(j,+4)+msq(j,+5)
       xmsq=xmsq+(msqv(j,g)
     & +msq(j,g)*(one
     &               +AP(q,q,1)-AP(q,q,3)+Q1(q,q,g,1)-Q1(q,q,g,3)
     &               +AP(g,g,1)-AP(g,g,3)+Q2(g,g,q,1)-Q2(g,g,q,3)))
     &               *fx1_L(j)*fx2_H(g)
     & +(msq(j,g)*(AP(q,q,2)+AP(q,q,3)+Q1(q,q,g,2)+Q1(q,q,g,3))
     & + msq(g,g)*(AP(g,q,2)+Q1(g,q,g,2)))*fx1z_L(j)/z*fx2_H(g)
     & +(msq(j,g)*(AP(g,g,2)+AP(g,g,3)+Q2(g,g,q,2)+Q2(g,g,q,3))
     & +   msq_qa*(AP(a,g,2)+Q2(a,g,q,2))
     & +   msq_qq*(AP(q,g,2)+Q2(q,g,q,2)))*fx1_L(j)*fx2z_H(g)/z
c--Qbarg
       elseif (j<0) then
       msq_aq=msq(j,+1)+msq(j,+2)+msq(j,+3)+msq(j,+4)+msq(j,+5)
       msq_aa=msq(j,-1)+msq(j,-2)+msq(j,-3)+msq(j,-4)+msq(j,-5)
       xmsq=xmsq+(msqv(j,g)
     & +msq(j,g)*(one+AP(a,a,1)-AP(a,a,3)+Q1(a,a,g,1)-Q1(a,a,g,3)
     &               +AP(g,g,1)-AP(g,g,3)+Q2(g,g,a,1)-Q2(g,g,a,3)))
     &                *fx1_L(j)*fx2_H(g)
     & +(msq(j,g)*(AP(a,a,2)+AP(a,a,3)+Q1(a,a,g,2)+Q1(a,a,g,3))
     & + msq(g,g)*(AP(g,a,2)+Q1(g,a,g,2)))*fx1z_L(j)/z*fx2_H(g)
     & +(msq(j,g)*(AP(g,g,2)+AP(g,g,3)+Q2(g,g,a,3)+Q2(g,g,a,2))
     & + msq_aq*(AP(q,g,2)+Q2(q,g,a,2))
     & + msq_aa*(AP(a,g,2)+Q2(a,g,a,2)))*fx1_L(j)*fx2z_H(g)/z
       endif
      endif

c--- END FOUR FLAVOR SINGLE TOP
      elseif (kcase == ktopanom) then ! let's just have a special case
        ! be very careful with modifications, since we don't follow fully the conventions used below

        if (j /= 0 .and. abs(k) == 5 ) then
            xmsq = xmsq + (msqv(j,k) + msq(j,k)) * singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy,k)

            ! WARNING: msqv_light and msqv_heavy do not get the msq
            ! subtraction on coeffonly! either add it, or use them for
            ! debugging only

c           xmsq = xmsq + msqv_light(j,k) * singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy,k)
            call singletop2_fillAP(z, i_beam1_light, AP)
            xmsq = xmsq + msq(j,k)*(AP(q,q,1)-AP(q,q,3)+B1(q,q,b,1)-B1(q,q,b,3)) *
     &              singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy,k)
            xmsq = xmsq + msq(j,k)*(AP(q,q,2)+AP(q,q,3)+B1(q,q,b,2)+B1(q,q,b,3)) *
     &              singletop2_pdfs(i_beam1_light_z,j)/z*singletop2_pdfs(i_beam2_heavy,k)

c           xmsq = xmsq + msqv_heavy(j,k) * singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy,k)
            call singletop2_fillAP(z, i_beam2_heavy, AP)
            xmsq = xmsq + msq(j,k)*(AP(q,q,1)-AP(q,q,3)+B2(b,b,q,1)-B2(b,b,q,3)) *
     &              singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy,k)
            xmsq = xmsq + msq(j,k)*(AP(q,q,2)+AP(q,q,3)+B2(b,b,q,2)+B2(b,b,q,3)) *
     &              singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy_z,k)/z
        elseif (abs(j) == 5 .and. k /= 0) then
            xmsq = xmsq + (msq(j,k) + msqv(j,k)) * singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light,k)

c           xmsq = xmsq + msqv_heavy(j,k) * singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light,k)
            call singletop2_fillAP(z, i_beam1_heavy, AP)
            xmsq = xmsq + msq(j,k)*(AP(q,q,1)-AP(q,q,3)+B1(b,b,q,1)-B1(b,b,q,3)) *
     &              singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light,k)
            xmsq = xmsq + msq(j,k)*(AP(q,q,2)+AP(q,q,3)+B1(b,b,q,2)+B1(b,b,q,3)) *
     &              singletop2_pdfs(i_beam1_heavy_z,j)/z * singletop2_pdfs(i_beam2_light,k)

c           xmsq = xmsq + msqv_light(j,k) * singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light,k)
            call singletop2_fillAP(z, i_beam2_light, AP)
            xmsq = xmsq + msq(j,k)*(AP(q,q,1)-AP(q,q,3)+B2(q,q,b,1)-B2(q,q,b,3)) *
     &              singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light,k)
            xmsq = xmsq + msq(j,k)*(AP(q,q,2)+AP(q,q,3)+B2(q,q,b,2)+B2(q,q,b,3)) *
     &              singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light_z,k)/z
        elseif (abs(j) /= 5 .and. k==0) then ! qg channel contribution
            call singletop2_fillAP(z, i_beam2_heavy, AP)
            xmsq = xmsq + (msq(j,+5) + msq(j,-5))*(AP(q,g,2) + B2(q,g,q,2)) *
     &              singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy_z,k)/z
        elseif (j==0 .and. abs(k) /= 5) then
            call singletop2_fillAP(z, i_beam1_heavy, AP)
            xmsq = xmsq + (msq(+5,k) + msq(-5,k))*(AP(q,g,2) + B1(q,g,q,2)) *
     &              singletop2_pdfs(i_beam1_heavy_z,j)/z * singletop2_pdfs(i_beam2_light,k)
        elseif (abs(j)==5 .and. k==0) then
            call singletop2_fillAP(z, i_beam2_light, AP)
            xmsq = xmsq + (sum(msq(j,1:5)) + sum(msq(j,-5:-1))) * (AP(q,g,2) + Q2(q,g,q,2)) *
     &              singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light_z,k)/z
        elseif (j==0 .and. abs(k)==5) then
            call singletop2_fillAP(z, i_beam1_light, AP)
            xmsq = xmsq + (sum(msq(1:5,k)) + sum(msq(-5:-1,k))) * (AP(q,g,2) + Q1(q,g,q,2)) *
     &          singletop2_pdfs(i_beam1_light_z,j)/z * singletop2_pdfs(i_beam2_heavy,k)
        else
c           write (*,*) "ERROR: bad channel for ktopanom in virtint.f: ", j,k
c           write(6,*) 'Abort in virtint'
c            stop
        endif

      else

c--- SUM BY TOTAL MATRIX ELEMENTS: everything else
c--- special code to remove the Q+Qb -> Z+Q+Qb contribution
      if ((j*k == -flav*flav) .and. (kcase==kgQ__ZQ)) then
        xmsq=xmsq-(
     & + msq(g,k)*(AP(g,q,2)+Q1(g,q,q,2))*fx1z(j)/z*fx2(k)
     & + msq(j,g)*(AP(g,q,2)+Q2(g,q,q,2))*fx1(j)*fx2z(k)/z)
      endif
c--- special code to remove the b+b(bar) -> W+t+b(bar) contribution
      if ( (abs(j*k) == nf*nf) .and.
     &     ((kcase==kW_tndk) .or. (kcase==kW_twdk)) ) then
        xmsq=xmsq-(
     & + msq(g,k)*(AP(g,q,2)+Q1(g,q,q,2))*fx1z(j)/z*fx2(k)
     & + msq(j,g)*(AP(g,q,2)+Q2(g,q,q,2))*fx1(j)*fx2z(k)/z)
      endif
c--QQ
      if     ((j > 0) .and. (k>0)) then
      if    ((((kcase==kbq_tpq) .or. (kcase==kH_tjet)
     &     .or.(kcase==kZ_tjet) .or. (kcase==kZ_tdkj)
     &     .or.(kcase==kH_tdkj) .or. (kcase==ktopanom)) .and. (j == 5))) then
      xmsq=xmsq+(msqv(j,k)
     & + msq(j,k)*(one+AP(q,q,1)-AP(q,q,3)+B1(b,b,q,1)-B1(b,b,q,3)
     &                +AP(q,q,1)-AP(q,q,3)+B2(q,q,b,1)-B2(q,q,b,3)))
     &                *fx1(j)*fx2(k)
     & +(msq(j,k)*(AP(q,q,2)+AP(q,q,3)+B1(b,b,q,2)+B1(b,b,q,3))
     & + msq(g,k)*(AP(g,q,2)+Q1(g,q,q,2)))*fx1z(j)/z*fx2(k)
     & +(msq(j,k)*(AP(q,q,2)+AP(q,q,3)+B2(q,q,b,2)+B2(q,q,b,3))
     & + msq(j,g)*(AP(g,q,2)+Q2(g,q,q,2)))*fx1(j)*fx2z(k)/z
      elseif((((kcase==kbq_tpq) .or. (kcase==kH_tjet)
     &    .or. (kcase==kZ_tjet) .or. (kcase==kZ_tdkj)
     &    .or. (kcase==kH_tdkj) .or. (kcase==ktopanom)) .and. (k == 5))) then
      xmsq=xmsq+(msqv(j,k)
     & + msq(j,k)*(one+AP(q,q,1)-AP(q,q,3)+B1(q,q,b,1)-B1(q,q,b,3)
     &                +AP(q,q,1)-AP(q,q,3)+B2(b,b,q,1)-B2(b,b,q,3)))
     &                *fx1(j)*fx2(k)
     & +(msq(j,k)*(AP(q,q,2)+AP(q,q,3)+B1(q,q,b,2)+B1(q,q,b,3))
     & + msq(g,k)*(AP(g,q,2)+Q1(g,q,q,2)))*fx1z(j)/z*fx2(k)
     & +(msq(j,k)*(AP(q,q,2)+AP(q,q,3)+B2(b,b,q,2)+B2(b,b,q,3))
     & + msq(j,g)*(AP(g,q,2)+Q2(g,q,q,2)))*fx1(j)*fx2z(k)/z
      else
      xmsq=xmsq+(msqv(j,k)
     & + msq(j,k)*(one+AP(q,q,1)-AP(q,q,3)+Q1(q,q,q,1)-Q1(q,q,q,3)
     &                +AP(q,q,1)-AP(q,q,3)+Q2(q,q,q,1)-Q2(q,q,q,3)))
     &                *fx1(j)*fx2(k)
     & +(msq(j,k)*(AP(q,q,2)+AP(q,q,3)+Q1(q,q,q,2)+Q1(q,q,q,3))
     & + msq(g,k)*(AP(g,q,2)+Q1(g,q,q,2)))*fx1z(j)/z*fx2(k)
     & +(msq(j,k)*(AP(q,q,2)+AP(q,q,3)+Q2(q,q,q,2)+Q2(q,q,q,3))
     & + msq(j,g)*(AP(g,q,2)+Q2(g,q,q,2)))*fx1(j)*fx2z(k)/z

      endif
c--QbarQbar
      elseif ((j < 0) .and. (k<0)) then
      if    ((((kcase==kbq_tpq) .or. (kcase==kH_tjet)
     &     .or.(kcase==kZ_tjet) .or. (kcase==kZ_tdkj)
     &    .or. (kcase==kH_tdkj) .or. (kcase==ktopanom)) .and. (j == -5))) then
      xmsq=xmsq+(msqv(j,k)
     & + msq(j,k)*(one+AP(a,a,1)-AP(a,a,3)+B1(b,b,a,1)-B1(b,b,a,3)
     &                +AP(a,a,1)-AP(a,a,3)+B2(a,a,b,1)-B2(a,a,b,3)))
     &                *fx1(j)*fx2(k)
     & +(msq(j,k)*(AP(a,a,2)+AP(a,a,3)+B1(b,b,a,2)+B1(b,b,a,3))
     & + msq(g,k)*(AP(g,a,2)+Q1(g,a,a,2)))*fx1z(j)/z*fx2(k)
     & +(msq(j,k)*(AP(a,a,2)+AP(a,a,3)+B2(a,a,b,2)+B2(a,a,b,3))
     & + msq(j,g)*(AP(g,a,2)+Q2(g,a,a,2)))*fx1(j)*fx2z(k)/z
      elseif((((kcase==kbq_tpq) .or. (kcase==kH_tjet)
     &     .or.(kcase==kZ_tjet) .or. (kcase==kZ_tdkj)
     &     .or.(kcase==kH_tdkj) .or. (kcase==ktopanom)) .and. (k == -5))) then
      xmsq=xmsq+(msqv(j,k)
     & + msq(j,k)*(one+AP(a,a,1)-AP(a,a,3)+B1(a,a,b,1)-B1(a,a,b,3)
     &                +AP(a,a,1)-AP(a,a,3)+B2(b,b,a,1)-B2(b,b,a,3)))
     &                *fx1(j)*fx2(k)
     & +(msq(j,k)*(AP(a,a,2)+AP(a,a,3)+B1(a,a,b,2)+B1(a,a,b,3))
     & + msq(g,k)*(AP(g,a,2)+Q1(g,a,a,2)))*fx1z(j)/z*fx2(k)
     & +(msq(j,k)*(AP(a,a,2)+AP(a,a,3)+B2(b,b,a,2)+B2(b,b,a,3))
     & + msq(j,g)*(AP(g,a,2)+Q2(g,a,a,2)))*fx1(j)*fx2z(k)/z
      else
      xmsq=xmsq+(msqv(j,k)
     & + msq(j,k)*(one+AP(a,a,1)-AP(a,a,3)+Q1(a,a,a,1)-Q1(a,a,a,3)
     &                +AP(a,a,1)-AP(a,a,3)+Q2(a,a,a,1)-Q2(a,a,a,3)))
     &                *fx1(j)*fx2(k)
     & +(msq(j,k)*(AP(a,a,2)+AP(a,a,3)+Q1(a,a,a,2)+Q1(a,a,a,3))
     & + msq(g,k)*(AP(g,a,2)+Q1(g,a,a,2)))*fx1z(j)/z*fx2(k)
     & +(msq(j,k)*(AP(a,a,2)+AP(a,a,3)+Q2(a,a,a,2)+Q2(a,a,a,3))
     & + msq(j,g)*(AP(g,a,2)+Q2(g,a,a,2)))*fx1(j)*fx2z(k)/z

      endif
c--QQbar
      elseif ((j > 0) .and. (k<0)) then
      xmsq=xmsq+(msqv(j,k)
     & + msq(j,k)*(one+AP(q,q,1)-AP(q,q,3)+Q1(q,q,a,1)-Q1(q,q,a,3)
     &                +AP(a,a,1)-AP(a,a,3)+Q2(a,a,q,1)-Q2(a,a,q,3)))
     &                *fx1(j)*fx2(k)
     & +(msq(j,k)*(AP(q,q,2)+AP(q,q,3)+Q1(q,q,a,3)+Q1(q,q,a,2))
     & + msq(g,k)*(AP(g,q,2)+Q1(g,q,a,2)))*fx1z(j)/z*fx2(k)
     & +(msq(j,k)*(AP(a,a,2)+AP(a,a,3)+Q2(a,a,q,3)+Q2(a,a,q,2))
     & + msq(j,g)*(AP(g,a,2)+Q2(g,a,q,2)))*fx1(j)*fx2z(k)/z

      elseif ((j < 0) .and. (k>0)) then
c--QbarQ
      xmsq=xmsq+(msqv(j,k)
     & +msq(j,k)*(one+AP(a,a,1)-AP(a,a,3)+Q1(a,a,q,1)-Q1(a,a,q,3)
     &               +AP(q,q,1)-AP(q,q,3)+Q2(q,q,a,1)-Q2(q,q,a,3)))
     &               *fx1(j)*fx2(k)
     & +(msq(j,k)*(AP(a,a,3)+AP(a,a,2)+Q1(a,a,q,3)+Q1(a,a,q,2))
     & + msq(g,k)*(AP(g,a,2)+Q1(g,a,q,2)))*fx1z(j)/z*fx2(k)
     & +(msq(j,k)*(AP(q,q,3)+AP(q,q,2)+Q2(q,q,a,3)+Q2(q,q,a,2))
     & + msq(j,g)*(AP(g,q,2)+Q2(g,q,a,2)))*fx1(j)*fx2z(k)/z

      elseif ((j == g) .and. (k==g)) then
c--gg
      if (kcase==kZbbbar) then
       xmsq=xmsq+msqv(g,g)*fx1(g)*fx2(g)
       do cs=0,2
       xmsq=xmsq+(
     & +msq_cs(cs,g,g)*(one
     &  +AP(g,g,1)-AP(g,g,3)+R1(g,g,g,cs,1)-R1(g,g,g,cs,3)
     &  +AP(g,g,1)-AP(g,g,3)+R2(g,g,g,cs,1)-R2(g,g,g,cs,3))
     &                 )*fx1(g)*fx2(g)
     & +msq_cs(cs,g,g)*(AP(g,g,2)+AP(g,g,3)
     &                 +R1(g,g,g,cs,2)+R1(g,g,g,cs,3))*fx1z(g)/z*fx2(g)
     & +msq_cs(cs,g,g)*(AP(g,g,2)+AP(g,g,3)
     &                 +R2(g,g,g,cs,2)+R2(g,g,g,cs,3))*fx1(g)*fx2z(g)/z
        enddo
      else
       msq_qg=msq(+5,g)+msq(+4,g)+msq(+3,g)+msq(+2,g)+msq(+1,g)
     &       +msq(-5,g)+msq(-4,g)+msq(-3,g)+msq(-2,g)+msq(-1,g)
       msq_gq=msq(g,+5)+msq(g,+4)+msq(g,+3)+msq(g,+2)+msq(g,+1)
     &       +msq(g,-5)+msq(g,-4)+msq(g,-3)+msq(g,-2)+msq(g,-1)
       xmsq=xmsq+(msqv(g,g)
     &  +msq(g,g)*(one+AP(g,g,1)-AP(g,g,3)+Q1(g,g,g,1)-Q1(g,g,g,3)
     &                +AP(g,g,1)-AP(g,g,3)+Q2(g,g,g,1)-Q2(g,g,g,3)))
     &                *fx1(g)*fx2(g)
     &  +(msq(g,g)*(AP(g,g,2)+AP(g,g,3)+Q1(g,g,g,2)+Q1(g,g,g,3))
     &  +   msq_qg*(AP(q,g,2)+Q1(q,g,g,2)))*fx1z(g)/z*fx2(g)
     &  +(msq(g,g)*(AP(g,g,2)+AP(g,g,3)+Q2(g,g,g,2)+Q2(g,g,g,3))
     &  +   msq_gq*(AP(q,g,2)+Q2(q,g,g,2)))*fx1(g)*fx2z(g)/z

      endif
      elseif (j == g) then
c--gQ
       if    (k > 0) then
       msq_aq=msq(-1,k)+msq(-2,k)+msq(-3,k)+msq(-4,k)+msq(-5,k)
       msq_qq=msq(+1,k)+msq(+2,k)+msq(+3,k)+msq(+4,k)+msq(+5,k)
       xmsq=xmsq+(msqv(g,k)
     & +msq(g,k)*(one+AP(g,g,1)-AP(g,g,3)+Q1(g,g,q,1)-Q1(g,g,q,3)
     &               +AP(q,q,1)-AP(q,q,3)+Q2(q,q,g,1)-Q2(q,q,g,3)))
     &               *fx1(g)*fx2(k)
     & +(msq(g,k)*(AP(g,g,2)+AP(g,g,3)+Q1(g,g,q,2)+Q1(g,g,q,3))
     & +   msq_aq*(AP(a,g,2)+Q1(a,g,q,2))
     & +   msq_qq*(AP(q,g,2)+Q1(q,g,q,2)))*fx1z(g)/z*fx2(k)
     & +(msq(g,k)*(AP(q,q,2)+AP(q,q,3)+Q2(q,q,g,2)+Q2(q,q,g,3))
     & + msq(g,g)*(AP(g,q,2)+Q2(g,q,g,2)))*fx1(g)*fx2z(k)/z

       if     ((kcase==kbq_tpq) .and. (k  /=  5)
     &   .and. (masslessb .eqv. .false.)) then
c--- replace subtraction with a b in initial state by massive sub
        xmsq=xmsq-(
     & +   msq_aq*(AP(a,g,2)+Q1(a,g,q,2)
     &            +ason2pi*Tr*(z**2+omz**2)*log(facscale**2/mbsq))
     & +   msq_qq*(AP(q,g,2)+Q1(q,g,q,2)
     &            +ason2pi*Tr*(z**2+omz**2)*log(facscale**2/mbsq))
     &            )*fx1z(g)/z*fx2(k)
       endif
c--gQbar
       elseif (k<0) then
       msq_qa=msq(+1,k)+msq(+2,k)+msq(+3,k)+msq(+4,k)+msq(+5,k)
       msq_aa=msq(-1,k)+msq(-2,k)+msq(-3,k)+msq(-4,k)+msq(-5,k)
       xmsq=xmsq+(msqv(g,k)
     & +msq(g,k)*(one+AP(g,g,1)-AP(g,g,3)+Q1(g,g,a,1)-Q1(g,g,a,3)
     &               +AP(a,a,1)-AP(a,a,3)+Q2(a,a,g,1)-Q2(a,a,g,3)))
     &               *fx1(g)*fx2(k)
     & +(msq(g,k)*(AP(g,g,2)+AP(g,g,3)+Q1(g,g,a,2)+Q1(g,g,a,3))
     & +   msq_qa*(AP(q,g,2)+Q1(q,g,a,2))
     & +   msq_aa*(AP(a,g,2)+Q1(a,g,a,2)))*fx1z(g)/z*fx2(k)
     & +(msq(g,k)*(AP(a,a,2)+AP(a,a,3)+Q2(a,a,g,2)+Q2(a,a,g,3))
     & + msq(g,g)*(AP(g,a,2)+Q2(g,a,g,2)))*fx1(g)*fx2z(k)/z

       if     ((kcase==kbq_tpq) .and. (k  /=  -5)
     &   .and. (masslessb .eqv. .false.)) then
c--- replace subtraction with a b in initial state by massive sub
        xmsq=xmsq-(
     & +   msq_qa*(AP(q,g,2)+Q1(q,g,a,2)
     &            +ason2pi*Tr*(z**2+omz**2)*log(facscale**2/mbsq))
     & +   msq_aa*(AP(a,g,2)+Q1(a,g,a,2)
     &            +ason2pi*Tr*(z**2+omz**2)*log(facscale**2/mbsq))
     &            )*fx1z(g)/z*fx2(k)
       endif
       endif
c--Qg
      elseif (k == g) then
       if     (j>0) then
       msq_qa=msq(j,-1)+msq(j,-2)+msq(j,-3)+msq(j,-4)+msq(j,-5)
       msq_qq=msq(j,+1)+msq(j,+2)+msq(j,+3)+msq(j,+4)+msq(j,+5)
       xmsq=xmsq+(msqv(j,g)
     & +msq(j,g)*(one
     &               +AP(q,q,1)-AP(q,q,3)+Q1(q,q,g,1)-Q1(q,q,g,3)
     &               +AP(g,g,1)-AP(g,g,3)+Q2(g,g,q,1)-Q2(g,g,q,3)))
     &               *fx1(j)*fx2(g)
     & +(msq(j,g)*(AP(q,q,2)+AP(q,q,3)+Q1(q,q,g,2)+Q1(q,q,g,3))
     & + msq(g,g)*(AP(g,q,2)+Q1(g,q,g,2)))*fx1z(j)/z*fx2(g)
     & +(msq(j,g)*(AP(g,g,2)+AP(g,g,3)+Q2(g,g,q,2)+Q2(g,g,q,3))
     & +   msq_qa*(AP(a,g,2)+Q2(a,g,q,2))
     & +   msq_qq*(AP(q,g,2)+Q2(q,g,q,2)))*fx1(j)*fx2z(g)/z

       if     ((kcase==kbq_tpq) .and. (j  /=  5)
     &   .and. (masslessb .eqv. .false.)) then
c--- replace subtraction with a b in initial state by massive sub
        xmsq=xmsq-(
     & +   msq_qa*(AP(a,g,2)+Q2(a,g,q,2)
     &            +ason2pi*Tr*(z**2+omz**2)*log(facscale**2/mbsq))
     & +   msq_qq*(AP(q,g,2)+Q2(q,g,q,2)
     &            +ason2pi*Tr*(z**2+omz**2)*log(facscale**2/mbsq))
     &            )*fx1(j)*fx2z(g)/z
       endif
c--Qbarg
       elseif (j<0) then
       msq_aq=msq(j,+1)+msq(j,+2)+msq(j,+3)+msq(j,+4)+msq(j,+5)
       msq_aa=msq(j,-1)+msq(j,-2)+msq(j,-3)+msq(j,-4)+msq(j,-5)
       xmsq=xmsq+(msqv(j,g)
     & +msq(j,g)*(one+AP(a,a,1)-AP(a,a,3)+Q1(a,a,g,1)-Q1(a,a,g,3)
     &               +AP(g,g,1)-AP(g,g,3)+Q2(g,g,a,1)-Q2(g,g,a,3)))
     &                *fx1(j)*fx2(g)
     & +(msq(j,g)*(AP(a,a,2)+AP(a,a,3)+Q1(a,a,g,2)+Q1(a,a,g,3))
     & + msq(g,g)*(AP(g,a,2)+Q1(g,a,g,2)))*fx1z(j)/z*fx2(g)
     & +(msq(j,g)*(AP(g,g,2)+AP(g,g,3)+Q2(g,g,a,3)+Q2(g,g,a,2))
     & + msq_aq*(AP(q,g,2)+Q2(q,g,a,2))
     & + msq_aa*(AP(a,g,2)+Q2(a,g,a,2)))*fx1(j)*fx2z(g)/z

       if     ((kcase==kbq_tpq) .and. (j  /=  -5)
     &   .and. (masslessb .eqv. .false.)) then
c--- replace subtraction with a b in initial state by massive sub
        xmsq=xmsq-(
     & +   msq_aq*(AP(q,g,2)+Q2(q,g,a,2)
     &            +ason2pi*Tr*(z**2+omz**2)*log(facscale**2/mbsq))
     & +   msq_aa*(AP(a,g,2)+Q2(a,g,a,2)
     &            +ason2pi*Tr*(z**2+omz**2)*log(facscale**2/mbsq))
     &            )*fx1(j)*fx2z(g)/z
       endif
       endif
      endif

c      if (xmsq-tmp  /=  0._dp) write(6,*) j,k,'-> msqc = ',xmsq-tmp

      endif

c--- subtract off LO (we don't want it) for Wcs_ms case and for
c--- the comparison with C. Oleari's e+e- --> QQbg calculation
c      if ((kcase==kWcs_ms) .or. (runstring(1:5) == 'carlo')) then
c--- (MCFM_original)  if (kcase==kWcs_ms) then

c     if ((kcase==kWcs_ms) .or. (purevirt)) then
c       xmsq=xmsq-msq(j,k)*fx1(j)*fx2(k)
c     endif

      if (kewcorr /= knone) xmsq_noew=xmsq_noew+fx1(j)*fx2(k)*msq_noew(j,k)

      enddo
      enddo

c compute weights for scale variation, looping if necessary
      if (doscalevar .and. bin) then
        itrial=itrial-1
        if (itrial > 0) then
          xmsqvar(itrial)=xmsq
          goto 66
        endif
        if (abs(xmsq) > zip) then
          scalereweight(1:maxscalevar)=xmsqvar(1:maxscalevar)/xmsq
        else
          scalereweight(:)=1._dp
        endif
      endif
      ! at this point itrial = 0 and scalereweight is set when doscalevar = .true.

      if (currentPDF == 0) then
        virtint=flux*xjac*pswt*xmsq/BrnRat
        if (kewcorr /= knone) virtint_noew=flux*xjac*pswt*xmsq_noew/BrnRat
      endif

c--- loop over all PDF error sets, if necessary
      if ((maxPDFsets > 0) .and. bin) then
          if (currentPDF > 0) then
              pdfreweight(currentPDF) = (virtint - flux*xjac*pswt*xmsq/BrnRat)*wgt
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
      ! at this point currentPDF = maxPDFsets when doPDFerrors

c--- code to check that epsilon poles cancel
      if (nshot == 1) then
        if ((xmsq == 0._dp) .or. ieee_is_nan(xmsq)) goto 999
        xmsq_old=xmsq
        nshot=nshot+1
        epinv=0._dp
        epinv2=0._dp
        goto 12
      elseif (nshot == 2) then
        nshot=nshot+1

        if ((abs(xmsq_old/xmsq-1._dp) > 1.e-7_dp) .or. ieee_is_nan(xmsq_old/xmsq)) then
!$omp master
          if (rank == 0) then
            write(6,*) 'epsilon fails to cancel'
            write(6,*) 'xmsq (epinv=large) = ',xmsq_old
            write(6,*) 'xmsq (epinv=zero ) = ',xmsq
            write(6,*) 'fractional difference',xmsq/xmsq_old-1d0
            call flush(6)
            !pause
          endif
!$omp end master
          colourchoice=rvcolourchoice
        else
!$omp master
          if (rank == 0) then
            write(6,*) 'Poles cancelled!'
c            write(6,*) 'xmsq (epinv=large) = ',xmsq_old
c            write(6,*) 'xmsq (epinv=zero ) = ',xmsq
c            write(6,*) 'fractional difference',xmsq/xmsq_old-1d0
            call flush(6)
          endif
!$omp end master
          colourchoice=rvcolourchoice
        endif
      endif

      call getptildejet(0,pjet)

      call dotem(nvec,pjet,s)

      if (ieee_is_nan(virtint) .or. (.not. ieee_is_finite(virtint))) then
         write(6,*) 'Found virtint=',virtint
         !write(6,*) 'random #s',r
         virtint=zip
         goto 999
      endif

      val=virtint*wgt
      val2=val**2

      if (abs(val) > wtmax) then
        wtmax=abs(val)
      endif

      if (bin) then
c--- for EW corrections, make additional weight available inside common block
        if (kewcorr /= knone) then
          wt_noew=virtint_noew*wgt
        endif
        if (newStyleHistograms) then
            call nplotter_new(pjet,val)
        else
            call nplotter(pjet,val,val2,0)
        endif
      endif


      ! QandGflag part is also affected by this, so have it before
      if (includeTaucutgrid(0) .eqv. .false.) then
          virtint = 0._dp
      endif


c--- handle special case of Qflag and Gflag
      if (QandGflag) then
        QandGint=QandGint+virtint
        if ((Gflag) .and. (.not.(Qflag))) then
c--- go back for second pass (Qflag)
          Qflag=.true.
          Gflag=.false.
          goto 44
        else
c--- return both to .true. and assign value to virtint (to return to VEGAS)
          Qflag=.true.
          Gflag=.true.
          virtint=QandGint
        endif
      endif

      if (enable_reweight_user) then
          virtint = virtint * reweight_user(pjet)
      endif

      return

 999  continue
      virtint=0._dp

c--- safety catch
      if (QandGflag) then
        Qflag=.true.
        Gflag=.true.
      endif

      return
      end

