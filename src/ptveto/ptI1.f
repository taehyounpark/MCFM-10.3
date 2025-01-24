!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine pTI1(z,Lperp,LQ,I1)
!      Lperp=log(mu/pt)
!      LQ=log(nu/Q)
      implicit none
!     Implementation of Eqs. 3.10-3.13 of 2207.07037v1
      include 'types.f'
      include 'constants.f'
      include 'transitionlabels.f'
      include 'distributions.f'
c     delt=-1,plus=0,lpls=1,rglr=2
      real(dp):: I1(0:6,dmin:dmax)
c      I1, index one transition flavor
c      I1, index2,distribution type
      real(dp):: tp0(0:6,dmin:dmax)
      real(dp):: z,omz,Lperp,LQ

c     initialize to zero
      I1(:,:)=0

!      tp0(qqV,plus)=2._dp;tp0(qqV,rglr)=-1._dp-z
!      tp0(qg,rglr)=z**2+omz**2
!      tp0(gq,rglr)=(1._dp+omz**2)/z
!      tp0(gg,plus)=1._dp; tp0(gg,rglr)=-1._dp+omz*opzsq/z

      call tildep0(+z,tp0)
      omz=1._dp-z

!qq   I1qqx2=2*CF*Lperp*4*LQ*delta(1-z)-4*CF*Lperp*P(q,q,z)+2*CF*(1-z)
      I1(qq,delt)=8*CF*Lperp*LQ
      I1(qq,rglr)=-4*CF*tp0(qq,rglr)*Lperp+2*CF*omz
      I1(qq,plus)=-4*CF*tp0(qq,plus)*Lperp

!gq   I1gqx1=-4*CF*P(g,q,z)*Lperp+2*z*CF;
      I1(gq,rglr)=-4*CF*tp0(gq,rglr)*Lperp+2*CF*z

!     nf?????
!qg   I1qgx1=-4*TF*Lperp*P(q,g,z)+4*TF*z*(1-z)
      I1(qg,rglr)=-4*Lperp*TF*tp0(qg,rglr)+4*TF*z*omz

!gg   I1ggx2=-8*CA*P(g,g,z)*Lperp
      I1(gg,delt)=8*CA*Lperp*LQ
      I1(gg,rglr)=-8*CA*tp0(gg,rglr)*Lperp
      I1(gg,plus)=-8*CA*tp0(gg,plus)*Lperp
      return
      end




