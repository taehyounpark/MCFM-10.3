!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c#######################################################################
c#### next-to-leading order quark beam functions in qt scheme          #
c#### normalized to as/(4 pi) which has been extracted                 #
c#######################################################################
      subroutine GLYbeam1(ih,zin,xb,beam,ibeam)
        use LHAPDF
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'facscale.f'
      include 'scale.f'
      include 'tiny.f'
      include 'transitionlabels.f'
      include 'distributions.f'
c     dmin=-2,dmax=1,rglr=-2,delt=-1,plus=0,lpls=1
      real(dp), intent(in) ::zin,xb
      real(dp), intent(out) :: beam(-5:5,0:2)
      real(dp)::cI1(0:6,dmin:dmax,0:2),
     & fxq,fx(-5:5),fx0(-5:5),zb,jaco
      integer :: ih,j,ibeam

c--- changing variables so z integral ranges from xb to 1
      zb = (one - xb)*zin + xb
      jaco = abs(one - xb)

c--- catch numerical problems
      if (zb > one-tiny) then
        beam(:,:)=zip
        return
      endif

c calculate parton distribution
      call fdist(ih,xb,facscale,fx0,ibeam)     ! at xb
      call fdist(ih,xb/zb,facscale,fx,ibeam)   ! at xb/zb

c calculate modification needed for beam function
      call GLYI1(zb,cI1)

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
c delta function piece
      beam(j,:)=cI1(gg,delt,:)*fx0(0)                    
c plus piece
      beam(j,:)=beam(j,:)
     & +cI1(gg,plus,:)*((fx(0)/zb-fx0(0))/(one-zb)*jaco
     & +fx0(0)*log(one-xb))
c log plus piece
      beam(j,:)=beam(j,:)
     & +cI1(gg,lpls,:)*((fx(0)/zb-fx0(0))*log(one-zb)/(one-zb)*jaco
     & +fx0(0)*0.5_dp*log(one-xb)**2)
c no-distribution piece
      beam(j,:)=beam(j,:)
     & +(cI1(gg,rglr,:)*fx(0)+cI1(gq,rglr,:)*fxq)*jaco/zb
      else
c---------------------------------------------
c----- Calculation of quark beam function ----
c---------------------------------------------
c delta function piece
      beam(j,:)=cI1(qq,delt,:)*fx0(j)
c plus piece
      beam(j,:)=beam(j,:)
     & +cI1(qq,plus,:)*((fx(j)/zb-fx0(j))/(one-zb)*jaco
     & +fx0(j)*log(one-xb))
c log plus piece
      beam(j,:)=beam(j,:)
     & +cI1(qq,lpls,:)*((fx(j)/zb-fx0(j))*log(one-zb)/(one-zb)*jaco
     & +fx0(j)*0.5_dp*log(one-xb)**2)
c no-distribution piece
      beam(j,:)=beam(j,:)
     & +(cI1(qq,rglr,:)*fx(j)+cI1(qg,rglr,:)*fx(0))*jaco/zb
      endif
      enddo
c      write(6,*)
c      j=-2
c      write(6,*) 'j',j
c      write(6,*) 'beam(j,2)',beam(j,2)
c      write(6,*) 'fx0(j)',fx0(j)
c      write(6,*) 'beam(j,2)/fx0(j)',beam(j,2)/fx0(j)
c      j=+1
c      write(6,*) 'j',j
c      write(6,*) 'beam(j,2)',beam(j,2)
c      write(6,*) 'fx0(j)',fx0(j)
c      write(6,*) 'beam(j,2)/fx0(j)',beam(j,2)/fx0(j)
      return
      end


