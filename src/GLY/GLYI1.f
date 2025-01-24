!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine GLYI1(z,I1)
c     Elaborated using the form results of 1403.6451v2
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'zeta.f'
      include 'transitionlabels.f'
      include 'distributions.f'
c     dmin=-2,dmax=1,rglr=-2,delt=-1,plus=0,lpls=1
c     -2=finite,-1=delta=-1,0=L0(z),1=L1(z)
      real(dp):: I1(0:6,dmin:dmax,0:2)
c      I1, index one transition flavor
c      I1, index2,distribution type
c      I1, index3,power of Lb
      real(dp):: tp0(0:6,dmin:dmax)
      real(dp):: z

c     initialize to zero
      I1(:,:,:)=0

      call tildep0(+z,tp0)

cgg   I1ggx2=+Lb^2*(+delta([1-z])*CA);I1ggx1=+Lb*(-4*P(g,g,z)*CA);I1ggx0=-delta([1-z])*CA*zeta2;
      I1(gg,delt,2)=CA
      I1(gg,rglr,1)=-4*CA*tp0(gg,rglr)
      I1(gg,plus,1)=-4*CA*tp0(gg,plus)
      I1(gg,delt,0)=-CA*zeta2

cqq   I1qqx2=+Lb^2*(+delta([1-z])*CF);I1qqx1=+Lb*(-2*P(q,q,z)*CF);I1qqx0=+2*[1-z]*CF-delta([1-z])*CF*zeta2;
      I1(qq,delt,2)=CF
      I1(qq,rglr,1)=-2*CF*tp0(qq,rglr)
      I1(qq,plus,1)=-2*CF*tp0(qq,plus)
      I1(qq,delt,0)=-CF*zeta2
      I1(qq,rglr,0)=+2*CF*(1._dp-z)

cgq   I1gqx1=+Lb*(-2*P(g,q,z)*CF);I1gqx0=+2*z*CF;
      I1(gq,rglr,1)=-2*CF*tp0(gq,rglr)
      I1(gq,plus,1)=-2*CF*tp0(gq,plus)
      I1(gq,rglr,0)=2*CF*z

cqg   I1qgx1=+Lb*(-2*P(q,g,z)*TF);I1qgx0=+2*TF-2*P(q,g,z)*TF;
      I1(qg,rglr,1)=-2*TF*tp0(qg,rglr)
      I1(qg,plus,1)=-2*TF*tp0(qg,plus)
      I1(qg,rglr,0)=2*TF-2*TF*tp0(qg,rglr)
      I1(qg,plus,0)=-2*TF*tp0(qg,plus)
      return
      end




