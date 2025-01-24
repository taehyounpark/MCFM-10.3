!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c#######################################################################
c#### next-to-leading order quark beam functions in BNR scheme         #
c#### normalized to as/(4 pi) which has been extracted                 #
c#### beam(-5:5) contains the beam function                            #
c#######################################################################
      subroutine BNRptbeam1(ih,zin,xb,lnmuonpt,beam,ibeam)
      use LHAPDF
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'facscale.f'
      include 'tiny.f'
      include 'transitionlabels.f'
      include 'distributions.f'
c     dmin=-2,dmax=1,rglr=-2,delt=-1,plus=0,lpls=1
      real(dp), intent(in) ::zin,xb,lnmuonpt
      real(dp), intent(out) :: beam(-5:5)
      real(dp)::fxq,fx(-5:5),fx0(-5:5),zb,jaco,I1(0:6,dmin:dmax)

      integer :: ih,j,ibeam

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

c calculate modification needed for beam function in BNR scheme

      call I1bar(zb,lnmuonpt,I1)

c calculate quark+antiquark sum for gq element
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
      beam(j)=I1(gg,delt)*fx0(0)
! plus piece
      beam(j)=beam(j)
     &+I1(gg,plus)*((fx(0)/zb-fx0(0))/(one-zb)*jaco
     & +fx0(0)*log(one-xb))
! log plus piece
      beam(j)=beam(j)
     & +I1(gg,lpls)*((fx(0)/zb-fx0(0))*log(one-zb)/(one-zb)*jaco
     & +fx0(0)*0.5_dp*log(one-xb)**2)
! no-distribution piece
      beam(j)=beam(j)
     & +(I1(gg,rglr)*fx(0)+I1(gq,rglr)*fxq)*jaco/zb
      else
c---------------------------------------------
c----- Calculation of quark beam function ----
c---------------------------------------------
c delta function piece
      beam(j)=I1(qq,delt)*fx0(j)
c plus piece
      beam(j)=beam(j)
     & +I1(qq,plus)*((fx(j)/zb-fx0(j))/(one-zb)*jaco
     & +fx0(j)*log(one-xb))
c log plus piece
      beam(j)=beam(j)
     & +I1(qq,lpls)*((fx(j)/zb-fx0(j))*log(one-zb)/(one-zb)*jaco
     & +fx0(j)*0.5_dp*log(one-xb)**2)
c no-distribution piece
      beam(j)=beam(j)
     & +(I1(qq,rglr)*fx(j)+I1(qg,rglr)*fx(0))*jaco/zb

      endif
      enddo

c alternative without rescaling
c      beam(:,:)=beam(:,:)/jaco

      return
      end


