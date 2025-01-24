!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine wqq_sc(is,ic,ie,in,iq,iqb,p,msq)
      implicit none
      include 'types.f'
c Author: J. Campbell, April 2004.

c     s(-p1)+cbar(-p2) --> l(p3)+abar(p4)+q(p5)+qb(p6)
c     with c massive
c     multiplied by (((a+l)^2-M**2)/(a+l)^2)^2*g^4/gwsq^2/2
      include 'constants.f'
      include 'mxpart.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'zprods_decl.f'
      real(dp):: msq,p(mxpart,4),q(mxpart,4),mch
      real(dp):: s134,s156,se,sn,en,sq,sqb,qqb,cs,afac,prop,facqq,dot
      complex(dp):: ampmm,ampmp,amppm,amppp
      integer:: is,ic,ie,in,iq,iqb,nu,j

      se=dot(p,is,ie)
      sn=dot(p,is,in)
      en=dot(p,ie,in)
      sq=dot(p,is,iq)
      sqb=dot(p,is,iqb)
      qqb=dot(p,iq,iqb)

      cs=dot(p,ic,is)
      mch=sqrt(abs(p(ic,4)**2-p(ic,1)**2-p(ic,2)**2-p(ic,3)**2))

      s134=2._dp*(se+sn+en)-mch**2
      s156=2._dp*(sq+sqb+qqb)

c----define massless momenta
      afac=0.5_dp*mch**2/cs

      do nu=1,4
      do j=1,6
        if (j==ic) then
          q(j,nu)=p(ic,nu)-p(is,nu)*afac
        else
          q(j,nu)=p(j,nu)
        endif
      enddo
      enddo

      call spinoru(6,q,za,zb)

      ampmm =
     &  + s134**(-1) * ( 2*za(ie,ic)*za(iq,ic)*zb(is,in)*zb(iqb,ic) + 2
     &    *za(is,ie)*za(iq,ic)*zb(is,in)*zb(iqb,is)/za(is,ic)/zb(is,
     &    ic)*mch**2 + 2*za(iq,ie)*za(iq,ic)*zb(is,in)*zb(iq,iqb) - 2*
     &    za(iq,ie)*zb(is,in)*zb(iqb,is)/zb(is,ic)*mch**2 )
      ampmm = ampmm + s156**(-1) * (  - 2*za(ie,ic)*za(iq,is)*zb(is,in)
     &    *zb(iqb,is) - 2*za(ie,ic)*za(iq,iqb)*zb(iqb,is)*zb(iqb,in) )

      amppm =
     &  + s134**(-1) * ( 2*za(ie,ic)*za(iq,is)*zb(is,in)*zb(iqb,ic)
     &    /za(is,ic)*mch + 2*za(is,ie)*za(iq,is)*zb(is,in)*zb(iqb,is)
     &    /za(is,ic)**2/zb(is,ic)*mch**3 + 2*za(iq,ie)*za(iq,is)*zb(is,
     &    in)*zb(iq,iqb)/za(is,ic)*mch - 2*za(iq,ie)*zb(is,in)*zb(iqb,
     &    ic)*mch )
      amppm = amppm + s156**(-1) * ( 2*za(is,ie)*za(iq,is)*zb(is,in)*
     &    zb(iqb,is)/za(is,ic)*mch + 2*za(is,ie)*za(iq,iqb)*zb(iqb,is)*
     &    zb(iqb,in)/za(is,ic)*mch )

      ampmp =
     &  + s134**(-1) * ( 2*za(ie,ic)*za(iqb,ic)*zb(is,in)*zb(iq,ic) + 2
     &    *za(is,ie)*za(iqb,ic)*zb(is,in)*zb(iq,is)/za(is,ic)/zb(is,
     &    ic)*mch**2 - 2*za(iqb,ie)*za(iqb,ic)*zb(is,in)*zb(iq,iqb) - 2*
     &    za(iqb,ie)*zb(is,in)*zb(iq,is)/zb(is,ic)*mch**2 )
      ampmp = ampmp + s156**(-1) * ( 2*za(ie,ic)*za(iq,iqb)*zb(iq,is)*
     &    zb(iq,in) - 2*za(ie,ic)*za(iqb,is)*zb(is,in)*zb(iq,is) )

      amppp =
     &  + s134**(-1) * ( 2*za(ie,ic)*za(iqb,is)*zb(is,in)*zb(iq,ic)
     &    /za(is,ic)*mch + 2*za(is,ie)*za(iqb,is)*zb(is,in)*zb(iq,is)
     &    /za(is,ic)**2/zb(is,ic)*mch**3 - 2*za(iqb,ie)*za(iqb,is)*zb(
     &    is,in)*zb(iq,iqb)/za(is,ic)*mch - 2*za(iqb,ie)*zb(is,in)*zb(
     &    iq,ic)*mch )
      amppp = amppp + s156**(-1) * (  - 2*za(is,ie)*za(iq,iqb)*zb(iq,is
     &    )*zb(iq,in)/za(is,ic)*mch + 2*za(is,ie)*za(iqb,is)*zb(is,in)*
     &    zb(iq,is)/za(is,ic)*mch )

c-- insert the gluon propagator
      ampmm=ampmm/2._dp/qqb
      ampmp=ampmp/2._dp/qqb
      amppm=amppm/2._dp/qqb
      amppp=amppp/2._dp/qqb

c--- overall factor from sc_Wqq.frm
      prop=1._dp/((two*en-wmass**2)**2+wmass**2*wwidth**2)
      facqq=aveqq*xn*cf*(gsq*gwsq)**2/2._dp*prop

      msq=facqq*real(ampmm*conjg(ampmm)+ampmp*conjg(ampmp)
     &              +amppm*conjg(amppm)+amppp*conjg(amppp))

      return
      end
