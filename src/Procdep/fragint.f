!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function fragint(r,wgt)
          use ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
          use PDFerrors
        use LHAPDF
          use Scalevar
          use MCFMStorage
          use SCET
          use m_gencuts, only : enable_reweight_user, reweight_user
          use MCFMSetupPlots, only: nplotter_new
          use MCFMSettings
      implicit none
      real(dp):: fragint

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'debug.f'
      include 'vegas_common.f'
      include 'frag.f'
      include 'npart.f'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'sprods_com.f'
      include 'kprocess.f'
      include 'noglue.f'
      include 'masses.f'
      include 'maxwt.f'
      include 'wts_bypart.f'
      include 'nflav.f'
      include 'ipsgen.f'
      include 'dm_params.f'
      include 'lastphot.f'
      include 'x1x2.f'
      include 'bypart.f'
      include 'energy.f'
      include 'initialscales.f'
      include 'phasemin.f'
      include 'xmin.f'
      include 'qcdcouple.f'
      include 'couple.f'
      include 'nlooprun.f'
      include 'beamtype.f'
      integer:: j,k,sgnj,sgnk,nvec,itrial
      real(dp):: r(mxdim),wgt,pswt,alphas,msqtrial,xmsqvar(2),
     & p(mxpart,4),flux,
     & BrnRat,pjet(mxpart,4),val,val2,
     & xmsq,xmsqjk,W,msq(-nf:nf,-nf:nf),fx1(-nf:nf),fx2(-nf:nf),
     & fx1up(-nf:nf),fx2up(-nf:nf),fx1dn(-nf:nf),fx2dn(-nf:nf),
     & msqdips(-nf:nf,-nf:nf),p_phys(mxpart,4),
     & wt34,wt345,wtprop,s34,s345,dot,wtips(4)
      logical:: bin,includedipole
      external qqb_w_g,qqb_z1jet,qqb_dirgam,qqb_2j_t,qqb_2j_s,
     & qqb_z2jetx,qqb_zaj,qqb_dm_monojet,qqb_gmgmjt,qqb_dirgam_g,
     & qqb_trigam_g
      common/bin/bin
      common/BrnRat/BrnRat

c--- statement function
      wtprop(s34,wmass,wwidth)=(s34-wmass**2)**2+(wmass*wwidth)**2

      fragint=0._dp

      W=sqrts**2

      p(:,:)=0._dp
      pjet(:,:)=0._dp

      currentNd = 0
      currentPDF=0

      if (maxPDFsets > 0) pdfreweight(:) = 0._dp

c----------------------------- GENERATE PHASE SPACE ---------------------------


c--- need to do something special for W_2gam and Z_2gam due to ipsgen
      if ((kcase==kW_2gam) .or. (kcase==kZ_2gam)) then
        npart=4
        if  (ipsgen == 1) then
            call gen_Vphotons_jets(r,2,0,p,pswt,*999)
        elseif (((ipsgen == 2) .and. (kcase==kZ_2gam)) .or.
     &          ((ipsgen == 3) .and. (kcase==kW_2gam))) then
            call gen_Vphotons_jets_dkrad(r,2,0,p,pswt,*999)
        else
           write(6,*) 'Parameter ipsgen not allowed'
           write(6,*) 'ipsgen = ',ipsgen
           stop
        endif
        if (kcase==kW_2gam) then
c          if (vetow_2gam(p)) goto 999 ! partition PS according to ipsgen
          s34=2._dp*dot(p,3,4)
          s345=s34+2._dp*dot(p,3,5)+2._dp*dot(p,4,5)
          wt34=wtprop(s34,wmass,wwidth)
          wt345=wtprop(s345,wmass,wwidth)
          wtips(1)=wt345
          wtips(3)=wt34
          pswt=pswt*wtips(ipsgen)/(wtips(1)+wtips(3))
        endif
      else
c--- otherwise, use same PS generation as at LO
        call gen_lops(r,p,pswt,*999)
      endif

      if (all(.not. ieee_is_nan(p(1:npart+2,:))) .eqv. .false.) then
          if (debug) then
              write(6,*) 'Discarding NaN or infinite phase space point'
          endif
          goto 999
      endif

      z_frag=r(ndim)
      frag=.true.

c--------------------------------- PHASE SPACE CUTS ---------------------------

      nvec=npart+2
      call dotem(nvec,p,s)

c----reject event if any s(i,j) is too small
c      call smalls(s,npart,*999)
c----reject event if any tau is too small
      call smalltau(p,npart,*999)

c--- the generated phase space point is the one that should be used
c--- in the calculation of the matrix elements;
c--- the physical momenta (p_phys) correspond to rescaling one of the photons
c--- "lastphot" by z_frag; this array should be used for cuts, plotting, etc.
      p_phys(:,:)=p(:,:)
      p_phys(lastphot,:)=z_frag*p(lastphot,:)

c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead bail out
      includeTaucutgrid(0) = .true.
      if (includedipole(0,p_phys) .eqv. .false.) then
        goto 999
      endif

c--- cut on z_frag
      if((z_frag < 0.0001_dp) .or. (z_frag > 1._dp)) goto 999

  771 continue
      if (dynamicscale) then
         call scaleset(initscale,initfacscale,p_phys)
      else
          call usescales(initscale,initfacscale)
      endif

      if (doPDFAlphas) then
          call updateAlphas(scale)
      endif

      xx(1)=-2._dp*p(1,4)/sqrts
      xx(2)=-2._dp*p(2,4)/sqrts

c      if (debug) write(*,*) 'Reconstructed x1,x2 ',xx(1),xx(2)
      if ( (xx(1) > one)  .or. (xx(2) > one)
     & .or.(xx(1) < xmin) .or. (xx(2) < xmin)) goto 999

      if ((doScalevar .and. currentPDF == 0 .and. bin) .and. (foundpow(currentPart) .eqv. .false.)) then
        itrial=1
      endif
   66 continue

c-------------------------- CALCULATE MATRIX ELEMENTS ---------------------------

c--------------------------------------------------------
c------------ NEED TO PASS P, P_PHYS TO ALL FRAGDIPS ROUTINES NOW
c---------------------------------------------------------
c-----------------------------------------------------------

      if     (kcase==kWgamma) then
         call qqb_wgam_frag(p,msq)
         call qqb_wgam_fragdips(p,p_phys,qqb_w_g,msqdips)
      elseif (kcase==kZgamma) then
         call qqb_zgam_frag(p,msq)
         call qqb_zgam_fragdips(p,p_phys,qqb_z1jet,msqdips)
      elseif (kcase==kdirgam) then
c         call qqb_dirgam_frag(p,msq)
c         call qqb_dirgam_fragdips(p,qqb_2j_t,qqb_2j_s,msqdips)
         msqdips(:,:)=0._dp
c======== new format
         call qqb_dirgam_frag_combo(p,p_phys,msq)
      elseif (kcase==kgamgam) then
        call qqb_gamgam_frag(p,msq)
        call qqb_gamgam_fragdips(p,p_phys,qqb_dirgam,msqdips)
      elseif (kcase==ktrigam) then
        call qqb_trigam_frag(p,msq)
        call qqb_trigam_fragdips(p,p_phys,qqb_gmgmjt,msqdips)
      elseif (kcase==kfourga) then
         call qqb_fourgam_frag(p,msq)
         call qqb_fourgam_fragdips(p,p_phys,qqb_trigam_g,msqdips)
      elseif (kcase==kgmgmjt) then
         msqdips(:,:)=0._dp
c====== new format
        call qqb_gmgmjt_frag_combo(p,p_phys,msq)
c        call qqb_gmgmjt_fragdips(p,p_phys,msq,qqb_dirgam_g)
      elseif(kcase==kZ_2gam) then
         call qqb_zaa_frag(p,msq)
         call qqb_zaa_fragdips(p,p_phys,qqb_zaj,msqdips)
      elseif(kcase==kZgajet) then
         call qqb_zaj_frag(p,msq)
         call qqb_zaj_fragdips(p,p_phys,qqb_z2jetx,msqdips)
      elseif(kcase==kdm_gam) then
         call qqb_dm_monophot_frag(p,msq)
         call qqb_dm_monophot_fragdips(p,p_phys,qqb_dm_monojet,msqdips)
      elseif(kcase==kW_2gam) then
         stop
c         call qqb_waa_frag_combo(p,p_phys,msq)
         msqdips(:,:)=0._dp
      else
        write(6,*) 'Fragmentation MEs not available for this process.'
        write(6,*) 'kcase = ',kcase
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

c--------------------------------------- INCLUDE PDF ---------------------------

      currentPDF = 0
      currentNd = 0

c--- do not calculate the flux if we're only checking the volume
      flux=fbGeV2/(2._dp*xx(1)*xx(2)*W)

c--- initialize a PDF set here, if calculating errors
  777 continue
      xmsq=0._dp

c--- calculate PDF's
      if ((maxPDFsets > 0) .and. bin) then
        call fdist(ih1,xx(1),facscale,fx1,1)
        call fdist(ih2,xx(2),facscale,fx2,2)
        ! this covers the case of PDF error AND scale variation for central PDF
        if (doScalevar .and. currentPDF == 0 .and. bin) then
            call fdist(ih1,xx(1),facscale*two,fx1up,1)
            call fdist(ih2,xx(2),facscale*two,fx2up,2)
            call fdist(ih1,xx(1),facscale/two,fx1dn,1)
            call fdist(ih2,xx(2),facscale/two,fx2dn,2)
            xmsqvar(:)=zip
        endif
      else
        call fdist(ih1,xx(1),facscale,fx1,1)
        call fdist(ih2,xx(2),facscale,fx2,2)
        if (doScalevar .and. currentPDF == 0 .and. bin) then
          call fdist(ih1,xx(1),facscale*two,fx1up,1)
          call fdist(ih2,xx(2),facscale*two,fx2up,2)
          call fdist(ih1,xx(1),facscale/two,fx1dn,1)
          call fdist(ih2,xx(2),facscale/two,fx2dn,2)
          xmsqvar(:)=zip
        endif
      endif

      do j=-nflav,nflav
      do k=-nflav,nflav

      if (ggonly) then
      if ((j /= 0) .or. (k /= 0)) cycle
      endif

      if (gqonly) then
      if (((j==0).and.(k==0)) .or. ((j /= 0).and.(k /= 0))) cycle
      endif

      if (noglue) then
      if ((j==0) .or. (k==0)) cycle
      endif

      if (omitgg) then
      if ((j==0) .and. (k==0)) cycle
      endif

c--- sum of fragmentation contribution and integrated fragmentation dipoles
      xmsqjk=fx1(j)*fx2(k)*(msq(j,k)+msqdips(j,k))
      if (doScalevar .and. currentPDF == 0 .and. bin) then
        xmsqvar(1)=xmsqvar(1)+fx1up(j)*fx2up(k)*(msq(j,k)+msqdips(j,k))
        xmsqvar(2)=xmsqvar(2)+fx1dn(j)*fx2dn(k)*(msq(j,k)+msqdips(j,k))
      endif

      xmsq=xmsq+xmsqjk

      if     (j > 0) then
        sgnj=+1
      elseif (j < 0) then
        sgnj=-1
      else
        sgnj=0
      endif
      if     (k > 0) then
        sgnk=+1
      elseif (k < 0) then
        sgnk=-1
      else
        sgnk=0
      endif

      enddo
      enddo

      if (currentPDF == 0) then
        fragint=flux*pswt*xmsq/BrnRat
      endif

c compute weights for scale variation
      if (doScalevar .and. currentPDF == 0 .and. bin) then
        if (abs(xmsq) > zip) then
          scalereweight(1)=(alphas(scale*two,amz,nlooprun)/as)**alphaspow(currentPart)
          scalereweight(2)=(alphas(scale/two,amz,nlooprun)/as)**alphaspow(currentPart)
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
          scalereweight(:)=zip
        endif
      endif

c--- loop over all PDF error sets, if necessary
      if ((maxPDFsets > 0) .and. bin) then
          if (currentPDF > 0) then
              pdfreweight(currentPDF) = (fragint - flux*pswt*xmsq/BrnRat)*wgt
          endif

          currentPDF = currentPDF + 1
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

      val=fragint*wgt
      val2=val**2
      if (ieee_is_nan(val)) then
         write(6,*) 'fragint val = ',val
         write(6,*) 'Discarding point with random variables',r
         fragint=zip
         val=zip
         goto 999
      endif

      if ((abs(val) > wtmax)) then
        wtmax=abs(val)
      endif


      if (bin) then
        if (newStyleHistograms) then
            call nplotter_new(pjet,val)
        else
          call nplotter(pjet,val,val2,0)
        endif
      endif

c --- Check weights :
c      if (unweight) then
c       write(6,*) 'Test weight ',val,' against max ',wtmax
c        wtabs = abs(val)
c        if (ran2() < (wtabs/wtmax)) then
c         write(6,*) 'Keep event with weight',val
c          if (wtabs<wtmax) then
c            newwt = 1._dp
c          else
c            newwt = wtabs/wtmax
c          endif
c          if (newwt > 1.0_dp) then
c            write(6,*) 'WARNING : fragint : event with |weight| > 1.',
c     +            ' |weight| = ',newwt
c          endif
c ---     just in case the weight was negative :
c          newwt = newwt*sign(1._dp,val)
c          call nplotter(pjet,newwt,newwt,0)
c ---     DSW. If I'm storing the event, I need to make a decision
c ---     about the flavours :
c          call decide_flavour(pflav,pbarflav)
c          call storeevent(pjet,newwt,pflav,pbarflav)
c        endif
c      endif

      if (enable_reweight_user) then
          fragint = fragint * reweight_user(pjet)
      endif

      return

 999  continue
      fragint=0._dp

      return
      end

