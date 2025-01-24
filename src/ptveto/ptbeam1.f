!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c#######################################################################
c#### next-to-leading order quark beam functions in qt scheme          #
c#### normalized to as/(4 pi) which has been extracted                 #
c#######################################################################
      subroutine ptbeam1(ih,zin,xb,Lperp,beam,ibeam)
        use LHAPDF
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'facscale.f'
      include 'tiny.f'
      include 'transitionlabels.f'
      include 'Lw.f'
      include 'Lnu.f'
      include 'scale.f'
      include 'energy.f'
      include 'distributions.f'
c     dmin=-2,dmax=1,rglr=-2,delt=-1,plus=0,lpls=1
      real(dp), intent(in) ::zin,xb,Lperp
      real(dp), intent(out) :: beam(-5:5)
      real(dp)::cI1(0:6,dmin:dmax,0:1),
     & fxq,fx(-5:5),fx0(-5:5),zb,jaco,LQ,
     & I1p(0:6,dmin:dmax,0:1),
     & I1p0(0:6,dmin:dmax),I1p1(0:6,dmin:dmax),
     & I1perp(0:6,dmin:dmax)
      integer,parameter::qqV=1,qbV=2,qqS=3

      integer :: ih,j,ibeam

      logical, parameter :: check = .false.

c--- changing variables so z integral ranges from xb to 1
      zb = (one - xb)*zin + xb
      jaco = abs(one - xb)

c alternative without rescaling
c      zb = zin
c      if (zb < xb) then
c        beam(:)=zip
c        return
c      endif

c--- catch numerical problems
      if (zb > one-tiny) then
        beam(:)=zip
        return
      endif

c calculate parton distribution
      call fdist(ih,xb,facscale,fx0,ibeam)     ! at xb
      call fdist(ih,xb/zb,facscale,fx,ibeam)   ! at xb/zb

c LQ = log(nu/w) = log(mu/w) - log(mu/nu) = log(mu/w) - Lnu
c where w = x*sqrts
      LQ = log(scale/(xb*sqrts))-Lnu
c calculate modification needed for beam function

      call ptI1(zb,Lperp,LQ,I1perp)

! Cross-check
      if (check) then
        call ptI1(zb, 0._dp,LQ,I1p0)
        call ptI1(zb, 1._dp,LQ,I1p1)
        I1p(:,:,0)=I1p0(:,:)
        I1p(:,:,1)=I1p1(:,:)-I1p(:,:,0)
! Lperp is log(mu/pt) but Lb is log(musq/bTsq)
        I1p(:,:,1)=I1p(:,:,1)/2._dp
        I1p(:,:,0)=I1p(:,:,0)
        Lw=LQ
        call calI1(zb,cI1)
        write(6,*) 'xb',xb
        write(6,*) 'zb',zb
        write(6,*) 'LQ',LQ
        write(6,*) 'Lnu',Lnu
        write(6,*) '============================================='
        write(6,*) 'cI1(gg,delt,1)',cI1(gg,delt,1)/I1p(gg,delt,1)
        write(6,*) 'cI1(gg,rglr,1)',cI1(gg,rglr,1)/I1p(gg,rglr,1)
        write(6,*) 'cI1(gg,plus,1)',cI1(gg,plus,1)/I1p(gg,plus,1)

        write(6,*) '============================================='
        write(6,*) 'cI1(gq,rglr,1)',cI1(gq,rglr,1)/I1p(gq,rglr,1)
        write(6,*) 'cI1(gq,rglr,0)',cI1(gq,rglr,0)/I1p(gq,rglr,0)

        write(6,*) '============================================='
        write(6,*) 'cI1(qg,rglr,1)',cI1(qg,rglr,1)/I1p(qg,rglr,1)
        write(6,*) 'cI1(qg,rglr,0)',cI1(qg,rglr,0)/I1p(qg,rglr,0)

        write(6,*) '============================================='
        write(6,*) 'cI1(qq,delt,1)',cI1(qqV,delt,1)/I1p(qq,delt,1)
        write(6,*) 'cI1(qq,rglr,1)',cI1(qqV,rglr,1)/I1p(qq,rglr,1)
        write(6,*) 'cI1(qq,rglr,0)',cI1(qqV,rglr,0)/I1p(qq,rglr,0)
        write(6,*) 'cI1(qq,plus,1)',cI1(qqV,plus,1)/I1p(qq,plus,1)
        pause
      endif


c calculate quark+antiquark sum
      fxq = zip
      do j = 1,nf
         fxq=fxq+fx(j)+fx(-j)
      enddo

      do j=-nf,nf
      if (j == 0) then
c---------------------------------------------
c----- Calculation of gluon beam function ----
c---------------------------------------------
! delta function piece
      beam(j)=I1perp(gg,delt)*fx0(0)
! plus piece
      beam(j)=beam(j)
     &+I1perp(gg,plus)*((fx(0)/zb-fx0(0))/(one-zb)*jaco
     & +fx0(0)*log(one-xb))
! log plus piece
      beam(j)=beam(j)
     & +I1perp(gg,lpls)*((fx(0)/zb-fx0(0))*log(one-zb)/(one-zb)*jaco
     & +fx0(0)*0.5_dp*log(one-xb)**2)
! no-distribution piece
      beam(j)=beam(j)
     & +(I1perp(gg,rglr)*fx(0)+I1perp(gq,rglr)*fxq)*jaco/zb
      else
c---------------------------------------------
c----- Calculation of quark beam function ----
c---------------------------------------------
c delta function piece
      beam(j)=I1perp(qq,delt)*fx0(j)
c plus piece
      beam(j)=beam(j)
     & +I1perp(qq,plus)*((fx(j)/zb-fx0(j))/(one-zb)*jaco
     & +fx0(j)*log(one-xb))
c log plus piece
      beam(j)=beam(j)
     & +I1perp(qq,lpls)*((fx(j)/zb-fx0(j))*log(one-zb)/(one-zb)*jaco
     & +fx0(j)*0.5_dp*log(one-xb)**2)
c no-distribution piece
      beam(j)=beam(j)
     & +(I1perp(qq,rglr)*fx(j)+I1perp(qg,rglr)*fx(0))*jaco/zb

      endif
      enddo

c alternative without rescaling
c      beam(:,:)=beam(:,:)/jaco

      return
      end


