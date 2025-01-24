!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_ZH1jet_v(p,msq)
      implicit none
      include 'types.f'

c***********************************************************************
c     Authors: R.K. Ellis and John Campbell                            *
c     May, 2001.                                                       *
c     Matrix element for Z + jet production                            *
c     in order alpha_s^2                                               *
c     averaged over initial colours and spins                          *
c     q(-p1)+qbar(-p2)-->Z^+(l(p3)+a(p4))+H(p5,p6)+g(p7)               *
c***********************************************************************

c Note that some Yukawa-type contributions are included, and axial-anomaly diagrams,
c but we neglect box contributions with external legs (g^* -> q q~, g, Z, H) that decouple

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'hdecaymode.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'epinv.f'
      include 'scheme.f'
      include 'nflav.f'
c      include 'cutoff.f'
      include 'hbbparams.f'
      include 'noglue.f'
      include 'zcouple_cms.f'
      include 'blha.f'
      include 'ewinput.f'
      integer:: j,k,ig
      real(dp):: msq(-nf:nf,-nf:nf),msq0(-nf:nf,-nf:nf),
     & p(mxpart,4),fac,s34,s127,virt5_VH,virt5_ZHtop,subuv,
     & msqhtautau,msqhbb,msqhgamgam,hdecay
      real(dp):: qqbZgLL,qqbZgRR,qqbZgLR,qqbZgRL
      real(dp):: gqZqLL,gqZqRR,gqZqLR,gqZqRL
      real(dp):: qgZqLL,qgZqRR,qgZqLR,qgZqRL
      real(dp):: qbqZgLL,qbqZgRR,qbqZgLR,qbqZgRL
      real(dp):: gqbZqbLL,gqbZqbRR,gqbZqbLR,gqbZqbRL
      real(dp):: qbgZqbLL,qbgZqbRR,qbgZqbLR,qbgZqbRL
      complex(dp):: qqbZgLLtop,qqbZgRRtop,qqbZgLRtop,qqbZgRLtop
      complex(dp):: gqZqLLtop,gqZqRRtop,gqZqLRtop,gqZqRLtop
      complex(dp):: qgZqLLtop,qgZqRRtop,qgZqLRtop,qgZqRLtop
      complex(dp):: qbqZgLLtop,qbqZgRRtop,qbqZgLRtop,qbqZgRLtop
      complex(dp):: gqbZqbLLtop,gqbZqbRRtop,gqbZqbLRtop,gqbZqbRLtop
      complex(dp):: qbgZqbLLtop,qbgZqbRRtop,qbgZqbLRtop,qbgZqbRLtop
      complex(dp):: prop,proptop,virt5ax,zsin2winv
      complex(dp):: qqbZgLLax,qqbZgRRax,qqbZgLRax,qqbZgRLax
      complex(dp):: qbqZgLLax,qbqZgRRax,qbqZgLRax,qbqZgRLax
      complex(dp):: qgZqLLax,qgZqRRax,qgZqLRax,qgZqRLax
      complex(dp):: gqbZqbLLax,gqbZqbRRax,gqbZqbLRax,gqbZqbRLax
      complex(dp):: gqZqLLax,gqZqRRax,gqZqLRax,gqZqRLax
      complex(dp):: qbgZqbLLax,qbgZqbRRax,qbgZqbLRax,qbgZqbRLax
      real(dp):: fac_top,sH!,cutoff_orig
      integer,parameter::
     & iqqbgLL(5)=(/1,2,3,4,7/),iqqbgLR(5)=(/1,2,4,3,7/),
     & iqgqLL(5)=(/1,7,3,4,2/),iqgqLR(5)=(/1,7,4,3,2/),
     & igqqLL(5)=(/2,7,3,4,1/),igqqLR(5)=(/2,7,4,3,1/)
      integer,parameter::
     & iqqbgLL9(5)=(/1,2,3,4,9/),iqqbgLR9(5)=(/1,2,4,3,9/),
     & iqgqLL9(5)=(/1,9,3,4,2/),iqgqLR9(5)=(/1,9,4,3,2/),
     & igqqLL9(5)=(/2,9,3,4,1/),igqqLR9(5)=(/2,9,4,3,1/)
      include 'cplx.h'
      common/virt5ax/virt5ax
!$omp threadprivate(/virt5ax/)

      scheme='dred'
c--set msq=0 to initialize
      msq(:,:)=zip

c      cutoff_orig=cutoff

      if(hdecaymode /= 'wpwm') then
         ig=7
      else
         ig=9
      endif

c--calculate spinor and dot-products (using BDK type notation)
      call spinoru(ig,p,za,zb)

c--- calculate lowest order
      msq0=zip
      if (toponly .eqv. .false.) then
        call qqb_ZH1jet(p,msq0)
      endif

c----UV counterterm contains the finite renormalization to arrive
c----at MS bar scheme.
      subuv=ason2pi*xn*(epinv*(11._dp-two*real(nflav,dp)/xn)-one)/6._dp

c--   calculate propagator
      s34=s(3,4)
      s127=s(1,2)+s(1,ig)+s(2,ig)
      proptop=s34/cplx2((s34-zmass**2),zmass*zwidth)
      prop=proptop*s127/cplx2((s127-zmass**2),zmass*zwidth)

c     Deal with Higgs decay
      if (hdecaymode == 'tlta') then
         sH=s(5,6)+two*mtau**2
         hdecay=msqhtautau(sH)
      elseif (hdecaymode == 'bqba') then
         sH=s(5,6)+two*mb**2
         hdecay=msqhbb(sH)
c---  adjust for fixed H->bb BR if necessary
         if (FixBrHbb) then
            hdecay=hdecay*GamHbb/GamHbb0
         endif
      elseif (hdecaymode == 'gaga') then
         sH=s(5,6)
         hdecay=msqhgamgam(sH)
      elseif (hdecaymode == 'wpwm') then
         sH=s(5,6)+s(5,7)+s(5,8)+s(6,7)+s(6,8)+s(7,8)
         call hwwdecay(p,5,6,7,8,hdecay)
      elseif (hdecaymode == 'none') then
         sH=s(5,5)
         hdecay=one
      else
         write(6,*) 'Unimplemented process in qqb_higgs'
         stop
      endif
      if (hdecaymode /= 'none') then
         hdecay=hdecay/((sH-hmass**2)**2+(hmass*hwidth)**2)
      endif

      fac=8._dp*cf*xnsq*esq**2*gsq*hdecay
      if (ewscheme < 4) then
        fac=fac*gwsq*wmass**2/(one-xw)**2
      else
c-- note: in CMS, wmass**2 -> sqrt(wmass**4+(wmass*wwidth)**2)
        fac=fac*gwsq*sqrt(wmass**4+(wmass*wwidth)**2)/abs(cone-zxw)**2
      endif
      fac_top=two*cf*xn*hdecay*2._dp

c---- smalls cut for top loop piece stability (quick but dirty)
c      cutoff=1.E-3_dp
      call smalls(s,ig-2,*101)
      goto 102
 101  fac_top=zip
 102  continue
c      cutoff=cutoff_orig

      if(ig==7) then
      if ((useblha == 0).or.((blhafl(1) /= 0).and.(blhafl(2) /= 0))) then
      qqbZgLLtop=aveqq*fac_top*virt5_ZHtop(iqqbgLL,za,zb)
      qqbZgLRtop=aveqq*fac_top*virt5_ZHtop(iqqbgLR,za,zb)
      qqbZgRLtop=aveqq*fac_top*virt5_ZHtop(iqqbgLR,zb,za)
      qqbZgRRtop=aveqq*fac_top*virt5_ZHtop(iqqbgLL,zb,za)

      if (toponly .eqv. .false.) then
        qqbZgLL=aveqq*fac*virt5_VH(iqqbgLL,za,zb,.true.)/s34**2
        qqbZgLLax=-aveqq*fac*virt5ax/s34**2
        qqbZgLR=aveqq*fac*virt5_VH(iqqbgLR,za,zb,.true.)/s34**2
        qqbZgLRax=-aveqq*fac*virt5ax/s34**2
        qqbZgRL=aveqq*fac*virt5_VH(iqqbgLR,zb,za,.true.)/s34**2
        qqbZgRLax=aveqq*fac*virt5ax/s34**2
        qqbZgRR=aveqq*fac*virt5_VH(iqqbgLL,zb,za,.true.)/s34**2
        qqbZgRRax=aveqq*fac*virt5ax/s34**2
      else
        qqbZgLL=0._dp
        qqbZgLR=0._dp
        qqbZgRL=0._dp
        qqbZgRR=0._dp
        qqbZgLLax=0._dp
        qqbZgLRax=0._dp
        qqbZgRLax=0._dp
        qqbZgRRax=0._dp
      endif

      qbqZgLL=qqbZgRL
      qbqZgLR=qqbZgRR
      qbqZgRL=qqbZgLL
      qbqZgRR=qqbZgLR

      qbqZgLLax=-qqbZgRLax
      qbqZgLRax=-qqbZgRRax
      qbqZgRLax=-qqbZgLLax
      qbqZgRRax=-qqbZgLRax

      qbqZgLLtop=qqbZgRLtop
      qbqZgLRtop=qqbZgRRtop
      qbqZgRLtop=qqbZgLLtop
      qbqZgRRtop=qqbZgLRtop
      endif

      if ((useblha == 0).or.((blhafl(1) == 0).and.(blhafl(2) /= 0))) then
      gqZqLLtop=aveqg*fac_top*virt5_ZHtop(igqqLL,za,zb)
      gqZqLRtop=aveqg*fac_top*virt5_ZHtop(igqqLR,za,zb)
      gqZqRLtop=aveqg*fac_top*virt5_ZHtop(igqqLR,zb,za)
      gqZqRRtop=aveqg*fac_top*virt5_ZHtop(igqqLL,zb,za)

      if (toponly .eqv. .false.) then
        gqZqLL=aveqg*fac*virt5_VH(igqqLL,za,zb,.true.)/s34**2
        gqZqLLax=-aveqg*fac*virt5ax/s34**2
        gqZqLR=aveqg*fac*virt5_VH(igqqLR,za,zb,.true.)/s34**2
        gqZqLRax=-aveqg*fac*virt5ax/s34**2
        gqZqRL=aveqg*fac*virt5_VH(igqqLR,zb,za,.true.)/s34**2
        gqZqRLax=aveqg*fac*virt5ax/s34**2
        gqZqRR=aveqg*fac*virt5_VH(igqqLL,zb,za,.true.)/s34**2
        gqZqRRax=aveqg*fac*virt5ax/s34**2
      else
        gqZqLL=0._dp
        gqZqLR=0._dp
        gqZqRL=0._dp
        gqZqRR=0._dp
        gqZqLLax=0._dp
        gqZqLRax=0._dp
        gqZqRLax=0._dp
        gqZqRRax=0._dp
      endif

      gqbZqbRL=gqZqLL
      gqbZqbRR=gqZqLR
      gqbZqbLL=gqZqRL
      gqbZqbLR=gqZqRR

      gqbZqbRLax=-gqZqLLax
      gqbZqbRRax=-gqZqLRax
      gqbZqbLLax=-gqZqRLax
      gqbZqbLRax=-gqZqRRax

      gqbZqbRLtop=gqZqLLtop
      gqbZqbRRtop=gqZqLRtop
      gqbZqbLLtop=gqZqRLtop
      gqbZqbLRtop=gqZqRRtop
      endif

      if ((useblha == 0).or.((blhafl(1) /= 0).and.(blhafl(2) == 0))) then
      qgZqLLtop=aveqg*fac_top*virt5_ZHtop(iqgqLL,za,zb)
      qgZqLRtop=aveqg*fac_top*virt5_ZHtop(iqgqLR,za,zb)
      qgZqRLtop=aveqg*fac_top*virt5_ZHtop(iqgqLR,zb,za)
      qgZqRRtop=aveqg*fac_top*virt5_ZHtop(iqgqLL,zb,za)

      if (toponly .eqv. .false.) then
        qgZqLL=aveqg*fac*virt5_VH(iqgqLL,za,zb,.true.)/s34**2
        qgZqLLax=-aveqg*fac*virt5ax/s34**2
        qgZqLR=aveqg*fac*virt5_VH(iqgqLR,za,zb,.true.)/s34**2
        qgZqLRax=-aveqg*fac*virt5ax/s34**2
        qgZqRL=aveqg*fac*virt5_VH(iqgqLR,zb,za,.true.)/s34**2
        qgZqRLax=aveqg*fac*virt5ax/s34**2
        qgZqRR=aveqg*fac*virt5_VH(iqgqLL,zb,za,.true.)/s34**2
        qgZqRRax=aveqg*fac*virt5ax/s34**2
      else
        qgZqLL=0._dp
        qgZqLR=0._dp
        qgZqRL=0._dp
        qgZqRR=0._dp
        qgZqLLax=0._dp
        qgZqLRax=0._dp
        qgZqRLax=0._dp
        qgZqRRax=0._dp
      endif

      qbgZqbRL=qgZqLL
      qbgZqbRR=qgZqLR
      qbgZqbLL=qgZqRL
      qbgZqbLR=qgZqRR

      qbgZqbRLax=-qgZqLLax
      qbgZqbRRax=-qgZqLRax
      qbgZqbLLax=-qgZqRLax
      qbgZqbLRax=-qgZqRRax

      qbgZqbRLtop=qgZqLLtop
      qbgZqbRRtop=qgZqLRtop
      qbgZqbLLtop=qgZqRLtop
      qbgZqbLRtop=qgZqRRtop
      endif

      elseif(ig==9) then

         qqbZgLL=aveqq*fac_top*virt5_ZHtop(iqqbgLL9,za,zb)
         qqbZgLR=aveqq*fac_top*virt5_ZHtop(iqqbgLR9,za,zb)
         qqbZgRL=aveqq*fac_top*virt5_ZHtop(iqqbgLR9,zb,za)
         qqbZgRR=aveqq*fac_top*virt5_ZHtop(iqqbgLL9,zb,za)

         if (toponly .eqv. .false.) then
           qqbZgLL=aveqq*fac*virt5_VH(iqqbgLL9,za,zb,.true.)/s34**2
           qqbZgLLax=-aveqq*fac*virt5ax/s34**2
           qqbZgLR=aveqq*fac*virt5_VH(iqqbgLR9,za,zb,.true.)/s34**2
           qqbZgLRax=-aveqq*fac*virt5ax/s34**2
           qqbZgRL=aveqq*fac*virt5_VH(iqqbgLR9,zb,za,.true.)/s34**2
           qqbZgRLax=aveqq*fac*virt5ax/s34**2
           qqbZgRR=aveqq*fac*virt5_VH(iqqbgLL9,zb,za,.true.)/s34**2
           qqbZgRRax=aveqq*fac*virt5ax/s34**2
         else
           qqbZgLL=0._dp
           qqbZgLR=0._dp
           qqbZgRL=0._dp
           qqbZgRR=0._dp
           qqbZgLLax=0._dp
           qqbZgLRax=0._dp
           qqbZgRLax=0._dp
           qqbZgRRax=0._dp
         endif

         qbqZgLL=qqbZgRL
         qbqZgLR=qqbZgRR
         qbqZgRL=qqbZgLL
         qbqZgRR=qqbZgLR

         qbqZgLLax=-qqbZgRLax
         qbqZgLRax=-qqbZgRRax
         qbqZgRLax=-qqbZgLLax
         qbqZgRRax=-qqbZgLRax

         qbqZgLLtop=qqbZgRLtop
         qbqZgLRtop=qqbZgRRtop
         qbqZgRLtop=qqbZgLLtop
         qbqZgRRtop=qqbZgLRtop

         gqZqLL=aveqg*fac_top*virt5_ZHtop(igqqLL9,za,zb)
         gqZqLR=aveqg*fac_top*virt5_ZHtop(igqqLR9,za,zb)
         gqZqRL=aveqg*fac_top*virt5_ZHtop(igqqLR9,zb,za)
         gqZqRR=aveqg*fac_top*virt5_ZHtop(igqqLL9,zb,za)

         if (toponly .eqv. .false.) then
           gqZqLL=aveqg*fac*virt5_VH(igqqLL9,za,zb,.true.)/s34**2
           gqZqLLax=-aveqg*fac*virt5ax/s34**2
           gqZqLR=aveqg*fac*virt5_VH(igqqLR9,za,zb,.true.)/s34**2
           gqZqLRax=-aveqg*fac*virt5ax/s34**2
           gqZqRL=aveqg*fac*virt5_VH(igqqLR9,zb,za,.true.)/s34**2
           gqZqRLax=aveqg*fac*virt5ax/s34**2
           gqZqRR=aveqg*fac*virt5_VH(igqqLL9,zb,za,.true.)/s34**2
           gqZqRRax=aveqg*fac*virt5ax/s34**2
         else
           gqZqLL=0._dp
           gqZqLR=0._dp
           gqZqRL=0._dp
           gqZqRR=0._dp
           gqZqLLax=0._dp
           gqZqLRax=0._dp
           gqZqRLax=0._dp
           gqZqRRax=0._dp
         endif

         gqbZqbRL=gqZqLL
         gqbZqbRR=gqZqLR
         gqbZqbLL=gqZqRL
         gqbZqbLR=gqZqRR

         gqbZqbRLax=-gqZqLLax
         gqbZqbRRax=-gqZqLRax
         gqbZqbLLax=-gqZqRLax
         gqbZqbLRax=-gqZqRRax

         gqbZqbRLtop=gqZqLLtop
         gqbZqbRRtop=gqZqLRtop
         gqbZqbLLtop=gqZqRLtop
         gqbZqbLRtop=gqZqRRtop

         qgZqLL=aveqg*fac_top*virt5_ZHtop(iqgqLL9,za,zb)
         qgZqLR=aveqg*fac_top*virt5_ZHtop(iqgqLR9,za,zb)
         qgZqRL=aveqg*fac_top*virt5_ZHtop(iqgqLR9,zb,za)
         qgZqRR=aveqg*fac_top*virt5_ZHtop(iqgqLL9,zb,za)

         if (toponly .eqv. .false.) then
           qgZqLL=aveqg*fac*virt5_VH(iqgqLL9,za,zb,.true.)/s34**2
           qgZqLLax=-aveqg*fac*virt5ax/s34**2
           qgZqLR=aveqg*fac*virt5_VH(iqgqLR9,za,zb,.true.)/s34**2
           qgZqLRax=-aveqg*fac*virt5ax/s34**2
           qgZqRL=aveqg*fac*virt5_VH(iqgqLR9,zb,za,.true.)/s34**2
           qgZqRLax=aveqg*fac*virt5ax/s34**2
           qgZqRR=aveqg*fac*virt5_VH(iqgqLL9,zb,za,.true.)/s34**2
           qgZqRRax=aveqg*fac*virt5ax/s34**2
         else
           qgZqLL=0._dp
           qgZqLR=0._dp
           qgZqRL=0._dp
           qgZqRR=0._dp
           qgZqLLax=0._dp
           qgZqLRax=0._dp
           qgZqRLax=0._dp
           qgZqRRax=0._dp
         endif

         qbgZqbRL=qgZqLL
         qbgZqbRR=qgZqLR
         qbgZqbLL=qgZqRL
         qbgZqbLR=qgZqRR

         qbgZqbRLax=-qgZqLLax
         qbgZqbRRax=-qgZqLRax
         qbgZqbLLax=-qgZqRLax
         qbgZqbLRax=-qgZqRRax

         qbgZqbRLtop=qgZqLLtop
         qbgZqbRRtop=qgZqLRtop
         qbgZqbLLtop=qgZqRLtop
         qbgZqbLRtop=qgZqRRtop
      endif

      zsin2winv=zL(2)-zR(2)

      do j=-nflav,nflav
      do k=-nflav,nflav
      if( j  /=  0 .and. k  /=  0 .and. j  /=  -k) goto 19

      if     ((j == 0) .and. (k == 0)) then
          msq(j,k)=0._dp
      elseif ((j > 0) .and. (k < 0)) then
          msq(j,k)=+abs(zL(j)*zl1*prop)**2*qqbZgLL
     &             +abs(zR(j)*zr1*prop)**2*qqbZgRR
     &             +abs(zL(j)*zr1*prop)**2*qqbZgLR
     &             +abs(zR(j)*zl1*prop)**2*qqbZgRL
     &             -subuv*msq0(j,k)
     &             +real(zsin2winv*zl1*prop*(qqbZgLLax*conjg(zL(j)*zl1*prop)+qqbZgRLax*conjg(zR(j)*zl1*prop)),dp)
     &             +real(zsin2winv*zr1*prop*(qqbZgRRax*conjg(zR(j)*zr1*prop)+qqbZgLRax*conjg(zL(j)*zr1*prop)),dp)
     &             +real(conjg(zL(j)*zl1*prop)*(Q(j)*q1+zL(j)*zl1*proptop)*qqbZgLLtop,dp)
     &             +real(conjg(zR(j)*zr1*prop)*(Q(j)*q1+zR(j)*zr1*proptop)*qqbZgRRtop,dp)
     &             +real(conjg(zL(j)*zr1*prop)*(Q(j)*q1+zL(j)*zr1*proptop)*qqbZgLRtop,dp)
     &             +real(conjg(zR(j)*zl1*prop)*(Q(j)*q1+zR(j)*zl1*proptop)*qqbZgRLtop,dp)
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=+abs(zL(k)*zl1*prop)**2*qbqZgLL
     &             +abs(zR(k)*zr1*prop)**2*qbqZgRR
     &             +abs(zL(k)*zr1*prop)**2*qbqZgLR
     &             +abs(zR(k)*zl1*prop)**2*qbqZgRL
     &             -subuv*msq0(j,k)
     &             +real(zsin2winv*zl1*prop*(qbqZgLLax*conjg(zL(k)*zl1*prop)+qbqZgRLax*conjg(zR(k)*zl1*prop)),dp)
     &             +real(zsin2winv*zr1*prop*(qbqZgRRax*conjg(zR(k)*zr1*prop)+qbqZgLRax*conjg(zL(k)*zr1*prop)),dp)
     &             +real(conjg(zL(k)*zl1*prop)*(Q(k)*q1+zL(k)*zl1*proptop)*qbqZgLLtop,dp)
     &             +real(conjg(zR(k)*zr1*prop)*(Q(k)*q1+zR(k)*zr1*proptop)*qbqZgRRtop,dp)
     &             +real(conjg(zL(k)*zr1*prop)*(Q(k)*q1+zL(k)*zr1*proptop)*qbqZgLRtop,dp)
     &             +real(conjg(zR(k)*zl1*prop)*(Q(k)*q1+zR(k)*zl1*proptop)*qbqZgRLtop,dp)
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=+abs(zL(j)*zl1*prop)**2*qgZqLL
     &             +abs(zR(j)*zr1*prop)**2*qgZqRR
     &             +abs(zL(j)*zr1*prop)**2*qgZqLR
     &             +abs(zR(j)*zl1*prop)**2*qgZqRL
     &             -subuv*msq0(j,k)
     &             +real(zsin2winv*zl1*prop*(qgZqLLax*conjg(zL(j)*zl1*prop)+qgZqRLax*conjg(zR(j)*zl1*prop)),dp)
     &             +real(zsin2winv*zr1*prop*(qgZqRRax*conjg(zR(j)*zr1*prop)+qgZqLRax*conjg(zL(j)*zr1*prop)),dp)
     &             +real(conjg(zL(j)*zl1*prop)*(Q(j)*q1+zL(j)*zl1*proptop)*qgZqLLtop,dp)
     &             +real(conjg(zR(j)*zr1*prop)*(Q(j)*q1+zR(j)*zr1*proptop)*qgZqRRtop,dp)
     &             +real(conjg(zL(j)*zr1*prop)*(Q(j)*q1+zL(j)*zr1*proptop)*qgZqLRtop,dp)
     &             +real(conjg(zR(j)*zl1*prop)*(Q(j)*q1+zR(j)*zl1*proptop)*qgZqRLtop,dp)
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=+abs(zL(-j)*zl1*prop)**2*qbgZqbLL
     &             +abs(zR(-j)*zr1*prop)**2*qbgZqbRR
     &             +abs(zL(-j)*zr1*prop)**2*qbgZqbLR
     &             +abs(zR(-j)*zl1*prop)**2*qbgZqbRL
     &             -subuv*msq0(j,k)
     &             +real(zsin2winv*zl1*prop*(qbgZqbLLax*conjg(zL(-j)*zl1*prop)+qbgZqbRLax*conjg(zR(-j)*zl1*prop)),dp)
     &             +real(zsin2winv*zr1*prop*(qbgZqbRRax*conjg(zR(-j)*zr1*prop)+qbgZqbLRax*conjg(zL(-j)*zr1*prop)),dp)
     &             +real(conjg(zL(-j)*zl1*prop)*(Q(-j)*q1+zL(-j)*zl1*proptop)*qbgZqbLLtop,dp)
     &             +real(conjg(zR(-j)*zr1*prop)*(Q(-j)*q1+zR(-j)*zr1*proptop)*qbgZqbRRtop,dp)
     &             +real(conjg(zL(-j)*zr1*prop)*(Q(-j)*q1+zL(-j)*zr1*proptop)*qbgZqbLRtop,dp)
     &             +real(conjg(zR(-j)*zl1*prop)*(Q(-j)*q1+zR(-j)*zl1*proptop)*qbgZqbRLtop,dp)
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=+abs(zL(k)*zl1*prop)**2*gqZqLL
     &             +abs(zR(k)*zr1*prop)**2*gqZqRR
     &             +abs(zL(k)*zr1*prop)**2*gqZqLR
     &             +abs(zR(k)*zl1*prop)**2*gqZqRL
     &             -subuv*msq0(j,k)
     &             +real(zsin2winv*zl1*prop*(gqZqLLax*conjg(zL(k)*zl1*prop)+gqZqRLax*conjg(zR(k)*zl1*prop)),dp)
     &             +real(zsin2winv*zr1*prop*(gqZqRRax*conjg(zR(k)*zr1*prop)+gqZqLRax*conjg(zL(k)*zr1*prop)),dp)
     &             +real(conjg(zL(k)*zl1*prop)*(Q(k)*q1+zL(k)*zl1*proptop)*gqZqLLtop,dp)
     &             +real(conjg(zR(k)*zr1*prop)*(Q(k)*q1+zR(k)*zr1*proptop)*gqZqRRtop,dp)
     &             +real(conjg(zL(k)*zr1*prop)*(Q(k)*q1+zL(k)*zr1*proptop)*gqZqLRtop,dp)
     &             +real(conjg(zR(k)*zl1*prop)*(Q(k)*q1+zR(k)*zl1*proptop)*gqZqRLtop,dp)
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=+abs(zL(-k)*zl1*prop)**2*gqbZqbLL
     &             +abs(zR(-k)*zr1*prop)**2*gqbZqbRR
     &             +abs(zL(-k)*zr1*prop)**2*gqbZqbLR
     &             +abs(zR(-k)*zl1*prop)**2*gqbZqbRL
     &             -subuv*msq0(j,k)
     &             +real(zsin2winv*zl1*prop*(gqbZqbLLax*conjg(zL(-k)*zl1*prop)+gqbZqbRLax*conjg(zR(-k)*zl1*prop)),dp)
     &             +real(zsin2winv*zr1*prop*(gqbZqbRRax*conjg(zR(-k)*zr1*prop)+gqbZqbLRax*conjg(zL(-k)*zr1*prop)),dp)
     &             +real(conjg(zL(-k)*zl1*prop)*(Q(-k)*q1+zL(-k)*zl1*proptop)*gqbZqbLLtop,dp)
     &             +real(conjg(zR(-k)*zr1*prop)*(Q(-k)*q1+zR(-k)*zr1*proptop)*gqbZqbRRtop,dp)
     &             +real(conjg(zL(-k)*zr1*prop)*(Q(-k)*q1+zL(-k)*zr1*proptop)*gqbZqbLRtop,dp)
     &             +real(conjg(zR(-k)*zl1*prop)*(Q(-k)*q1+zR(-k)*zl1*proptop)*gqbZqbRLtop,dp)
      endif

   19 continue
      enddo
      enddo

      return
      end
