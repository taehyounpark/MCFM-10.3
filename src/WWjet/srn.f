!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine srn(p1,p2,p3,p4,p5,p6,p7,za,zb,vecm,Qid,swapen,unsym)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'nf.f'
      include 'zcouple_cms.f'
      include 'ewcharge.f'
      logical:: swapen
      integer:: p1,p2,p3,p4,p5,p6,p7,Qid
      real(dp):: s56,s3456,s356,s456,qq,qn,qe
      complex(dp)::unsym(2),PW,PZ,MZsq,MWsq,prop3456,
     & iza,izb,zba2,vqL,vqR,veL,vnL,vecm(mxpart,mxpart)
c     statement functions
      iza(p1,p2)=1._dp/za(p1,p2)
      izb(p1,p2)=1._dp/zb(p1,p2)
      zba2(p1,p2,p3,p4)=zb(p1,p2)*za(p2,p4)+zb(p1,p3)*za(p3,p4)
      PW(s56)=s56/(s56-MWsq)
      PZ(qq,qe,vqL,veL,prop3456)=qq*qe+vqL*veL*prop3456
c     end statement functions
      MWsq=cmplx(wmass**2,-wmass*wwidth,kind=dp)
      MZsq=cmplx(zmass**2,-zmass*zwidth,kind=dp)
      if (swapen) then
        qn=-1._dp
        qe=0._dp
        veL=zln
        vnL=zle
      else
        qe=-1._dp
        qn=0._dp
        vnL=zln
        veL=zle
      endif
      qq=Q(qid)
      vqL=zL(qid)
      vqR=zR(qid)
      s56=s(p5,p6)
      s356=s(p3,p5)+s(p3,p6)+s(p5,p6)
      s456=s(p4,p5)+s(p4,p6)+s(p5,p6)
      s3456=s(p3,p4)+s(p3,p5)+s(p3,p6)+s(p4,p5)+s(p4,p6)+s(p5,p6)
      prop3456=s3456/(s3456-MZsq)
      unsym(1)= + PZ(qq,qe,vqL,veL,prop3456)*PW(s56)*s3456**(-1)*s56**(-1)
     & *s456**(-1) * (  - 2*za(p1,p3)*zb(p4,p5)*iza(p1,p7)*izb(p1,p7)*
     &    zba2(p2,p4,p5,p6)*vecm(p1,p1) + 2*za(p1,p3)*zb(p4,p5)*iza(p2,
     &    p7)*izb(p2,p7)*zba2(p2,p4,p5,p6)*vecm(p2,p2) + 2*za(p1,p3)*
     &    zb(p4,p5)*iza(p2,p7)*izb(p2,p7)*zba2(p7,p4,p5,p6)*vecm(p7,p2)
     &     + 2*za(p3,p7)*zb(p4,p5)*iza(p1,p7)*izb(p1,p7)*zba2(p2,p4,p5,
     &    p6)*vecm(p1,p7) )
      unsym(1) = unsym(1) + PZ(qq,qn,vqL,vnL,prop3456)*PW(s56)*s3456**(-1)
     & *s56**(-1)*s356**(-1) * ( 2*za(p3,p6)*zb(p2,p4)*iza(p1,p7)*izb(
     &    p1,p7)*zba2(p5,p3,p6,p1)*vecm(p1,p1) + 2*za(p3,p6)*zb(p2,p4)*
     &    iza(p1,p7)*izb(p1,p7)*zba2(p5,p3,p6,p7)*vecm(p1,p7) - 2*za(p3
     &    ,p6)*zb(p2,p4)*iza(p2,p7)*izb(p2,p7)*zba2(p5,p3,p6,p1)*vecm(
     &    p2,p2) + 2*za(p3,p6)*zb(p4,p7)*iza(p2,p7)*izb(p2,p7)*zba2(p5,
     &    p3,p6,p1)*vecm(p7,p2) )
      unsym(2)= + PZ(qq,qe,vqR,veL,prop3456)*PW(s56)*s3456**(-1)*s56**(-1)
     & *s456**(-1) * (  - 2*za(p2,p3)*za(p4,p6)*zb(p1,p4)*zb(p4,p5)*
     &    iza(p1,p7)*izb(p1,p7)*vecm(p1,p1) + 2*za(p2,p3)*za(p4,p6)*zb(
     &    p4,p5)*zb(p4,p7)*iza(p1,p7)*izb(p1,p7)*vecm(p7,p1) - 2*za(p2,
     &    p3)*za(p5,p6)*zb(p1,p5)*zb(p4,p5)*iza(p1,p7)*izb(p1,p7)*vecm(
     &    p1,p1) + 2*za(p2,p3)*za(p5,p6)*zb(p4,p5)*zb(p5,p7)*iza(p1,p7)
     &    *izb(p1,p7)*vecm(p7,p1) + 2*za(p2,p3)*zb(p4,p5)*iza(p2,p7)*
     &    izb(p2,p7)*zba2(p1,p4,p5,p6)*vecm(p2,p2) - 2*za(p3,p7)*zb(p4,
     &    p5)*iza(p2,p7)*izb(p2,p7)*zba2(p1,p4,p5,p6)*vecm(p2,p7) )
      unsym(2) = unsym(2) + PZ(qq,qn,vqR,vnL,prop3456)*PW(s56)*s3456**(-1)
     & *s56**(-1)*s356**(-1) * ( 2*za(p2,p3)*za(p3,p6)*zb(p1,p4)*zb(p3,
     &    p5)*iza(p1,p7)*izb(p1,p7)*vecm(p1,p1) - 2*za(p2,p3)*za(p3,p6)
     &    *zb(p3,p5)*zb(p4,p7)*iza(p1,p7)*izb(p1,p7)*vecm(p7,p1) - 2*
     &    za(p2,p6)*za(p3,p6)*zb(p1,p4)*zb(p5,p6)*iza(p1,p7)*izb(p1,p7)
     &    *vecm(p1,p1) + 2*za(p2,p6)*za(p3,p6)*zb(p4,p7)*zb(p5,p6)*iza(
     &    p1,p7)*izb(p1,p7)*vecm(p7,p1) - 2*za(p3,p6)*zb(p1,p4)*iza(p2,
     &    p7)*izb(p2,p7)*zba2(p5,p3,p6,p2)*vecm(p2,p2) - 2*za(p3,p6)*
     &    zb(p1,p4)*iza(p2,p7)*izb(p2,p7)*zba2(p5,p3,p6,p7)*vecm(p2,p7)
     &     )
      return
      end
