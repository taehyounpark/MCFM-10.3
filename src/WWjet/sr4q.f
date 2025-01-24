!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine sr4q(p1,p2,p3,p4,p5,p6,p7,p8,Qid,swapen,unsym)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'nf.f'
      include 'zcouple_cms.f'
      include 'ewcharge.f'
      logical:: swapen
      integer:: p1,p2,p3,p4,p5,p6,p7,p8,Qid
      real(dp):: s56,s3456,s128,s278,s356,s456,s28,qq,qn,qe
      complex(dp):: unsym(2,2),PW,PZ,MZsq,MWsq,prop3456,
     & zba2,vqL,vqR,veL,vnL
c     statement functions
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
      s28=s(p2,p8)
      s56=s(p5,p6)
      s128=s(p1,p2)+s(p1,p8)+s(p2,p8)
      s278=s(p2,p7)+s(p2,p8)+s(p7,p8)
      s356=s(p3,p5)+s(p3,p6)+s(p5,p6)
      s456=s(p4,p5)+s(p4,p6)+s(p5,p6)
      s3456=s(p3,p4)+s(p3,p5)+s(p3,p6)+s(p4,p5)+s(p4,p6)+s(p5,p6)
      prop3456=s3456/(s3456-MZsq)
      unsym(1,1)= + PZ(qq,qe,vqL,veL,prop3456)*PW(s56)*s3456**(-1)*
     & s56**(-1)*s28**(-1)*s456**(-1) * ( 2*za(p1,p8)*za(p3,p7)*za(p4,
     &    p6)*zb(p1,p2)*zb(p1,p4)*zb(p4,p5)*s128**(-1) + 2*za(p1,p8)*
     &    za(p3,p7)*za(p5,p6)*zb(p1,p2)*zb(p1,p5)*zb(p4,p5)*s128**(-1)
     &     + 2*za(p2,p8)*za(p3,p7)*zb(p1,p2)*zb(p4,p5)*zba2(p2,p4,p5,p6
     &    )*s128**(-1) + 2*za(p3,p7)*za(p4,p6)*za(p7,p8)*zb(p1,p4)*zb(
     &    p2,p7)*zb(p4,p5)*s278**(-1) + 2*za(p3,p7)*za(p5,p6)*za(p7,p8)
     &    *zb(p1,p5)*zb(p2,p7)*zb(p4,p5)*s278**(-1) + 2*za(p3,p8)*za(p4
     &    ,p6)*za(p7,p8)*zb(p1,p4)*zb(p2,p8)*zb(p4,p5)*s278**(-1) + 2*
     &    za(p3,p8)*za(p5,p6)*za(p7,p8)*zb(p1,p5)*zb(p2,p8)*zb(p4,p5)*
     &    s278**(-1) )
      unsym(1,1) = unsym(1,1) + PZ(qq,qn,vqL,vnL,prop3456)*PW(s56)*
     & s3456**(-1)*s56**(-1)*s28**(-1)*s356**(-1) * (  - 2*za(p1,p8)*
     &    za(p3,p6)*za(p3,p7)*zb(p1,p2)*zb(p1,p4)*zb(p3,p5)*s128**(-1)
     &     + 2*za(p1,p8)*za(p3,p6)*za(p6,p7)*zb(p1,p2)*zb(p1,p4)*zb(p5,
     &    p6)*s128**(-1) + 2*za(p2,p8)*za(p3,p6)*zb(p1,p2)*zb(p2,p4)*
     &    zba2(p5,p3,p6,p7)*s128**(-1) - 2*za(p3,p6)*za(p3,p7)*za(p7,p8
     &    )*zb(p1,p4)*zb(p2,p7)*zb(p3,p5)*s278**(-1) - 2*za(p3,p6)*za(
     &    p3,p8)*za(p7,p8)*zb(p1,p4)*zb(p2,p8)*zb(p3,p5)*s278**(-1) + 2
     &    *za(p3,p6)*za(p6,p7)*za(p7,p8)*zb(p1,p4)*zb(p2,p7)*zb(p5,p6)*
     &    s278**(-1) + 2*za(p3,p6)*za(p6,p8)*za(p7,p8)*zb(p1,p4)*zb(p2,
     &    p8)*zb(p5,p6)*s278**(-1) )
      unsym(2,1)= + PZ(qq,qe,vqR,veL,prop3456)*PW(s56)*s3456**(-1)*
     & s56**(-1)*s28**(-1)*s456**(-1) * ( 2*za(p1,p3)*za(p1,p8)*za(p4,
     &    p6)*zb(p1,p2)*zb(p4,p5)*zb(p4,p7)*s128**(-1) + 2*za(p1,p3)*
     &    za(p1,p8)*za(p5,p6)*zb(p1,p2)*zb(p4,p5)*zb(p5,p7)*s128**(-1)
     &     - 2*za(p1,p3)*za(p2,p8)*za(p4,p6)*zb(p2,p4)*zb(p2,p7)*zb(p4,
     &    p5)*s278**(-1) - 2*za(p1,p3)*za(p2,p8)*za(p5,p6)*zb(p2,p5)*
     &    zb(p2,p7)*zb(p4,p5)*s278**(-1) + 2*za(p1,p3)*za(p4,p6)*za(p7,
     &    p8)*zb(p2,p7)*zb(p4,p5)*zb(p4,p7)*s278**(-1) + 2*za(p1,p3)*
     &    za(p5,p6)*za(p7,p8)*zb(p2,p7)*zb(p4,p5)*zb(p5,p7)*s278**(-1)
     &     - 2*za(p1,p8)*za(p3,p8)*zb(p2,p8)*zb(p4,p5)*zba2(p7,p4,p5,p6
     &    )*s128**(-1) )
      unsym(2,1) = unsym(2,1) + PZ(qq,qn,vqR,vnL,prop3456)*PW(s56)*
     & s3456**(-1)*s56**(-1)*s28**(-1)*s356**(-1) * (  - 2*za(p1,p3)*
     &    za(p1,p8)*za(p3,p6)*zb(p1,p2)*zb(p3,p5)*zb(p4,p7)*s128**(-1)
     &     + 2*za(p1,p3)*za(p2,p8)*za(p3,p6)*zb(p2,p4)*zb(p2,p7)*zb(p3,
     &    p5)*s278**(-1) - 2*za(p1,p3)*za(p3,p6)*za(p7,p8)*zb(p2,p7)*
     &    zb(p3,p5)*zb(p4,p7)*s278**(-1) + 2*za(p1,p6)*za(p1,p8)*za(p3,
     &    p6)*zb(p1,p2)*zb(p4,p7)*zb(p5,p6)*s128**(-1) - 2*za(p1,p6)*
     &    za(p2,p8)*za(p3,p6)*zb(p2,p4)*zb(p2,p7)*zb(p5,p6)*s278**(-1)
     &     + 2*za(p1,p6)*za(p3,p6)*za(p7,p8)*zb(p2,p7)*zb(p4,p7)*zb(p5,
     &    p6)*s278**(-1) + 2*za(p1,p8)*za(p3,p6)*zb(p2,p8)*zb(p4,p7)*
     &    zba2(p5,p3,p6,p8)*s128**(-1) )
      unsym(1,2)= + PZ(qq,qe,vqL,veL,prop3456)*PW(s56)*s3456**(-1)*
     & s56**(-1)*s28**(-1)*s456**(-1) * ( 2*za(p1,p2)*za(p3,p7)*za(p4,
     &    p6)*zb(p1,p4)*zb(p1,p8)*zb(p4,p5)*s128**(-1) + 2*za(p1,p2)*
     &    za(p3,p7)*za(p5,p6)*zb(p1,p5)*zb(p1,p8)*zb(p4,p5)*s128**(-1)
     &     - 2*za(p2,p3)*za(p2,p7)*za(p4,p6)*zb(p1,p4)*zb(p2,p8)*zb(p4,
     &    p5)*s278**(-1) - 2*za(p2,p3)*za(p2,p7)*za(p5,p6)*zb(p1,p5)*
     &    zb(p2,p8)*zb(p4,p5)*s278**(-1) + 2*za(p2,p7)*za(p3,p7)*za(p4,
     &    p6)*zb(p1,p4)*zb(p4,p5)*zb(p7,p8)*s278**(-1) + 2*za(p2,p7)*
     &    za(p3,p7)*za(p5,p6)*zb(p1,p5)*zb(p4,p5)*zb(p7,p8)*s278**(-1)
     &     - 2*za(p2,p8)*za(p3,p7)*zb(p1,p8)*zb(p4,p5)*zba2(p8,p4,p5,p6
     &    )*s128**(-1) )
      unsym(1,2) = unsym(1,2) + PZ(qq,qn,vqL,vnL,prop3456)*PW(s56)*
     & s3456**(-1)*s56**(-1)*s28**(-1)*s356**(-1) * (  - 2*za(p1,p2)*
     &    za(p3,p6)*za(p3,p7)*zb(p1,p4)*zb(p1,p8)*zb(p3,p5)*s128**(-1)
     &     + 2*za(p1,p2)*za(p3,p6)*za(p6,p7)*zb(p1,p4)*zb(p1,p8)*zb(p5,
     &    p6)*s128**(-1) + 2*za(p2,p3)*za(p2,p7)*za(p3,p6)*zb(p1,p4)*
     &    zb(p2,p8)*zb(p3,p5)*s278**(-1) - 2*za(p2,p6)*za(p2,p7)*za(p3,
     &    p6)*zb(p1,p4)*zb(p2,p8)*zb(p5,p6)*s278**(-1) - 2*za(p2,p7)*
     &    za(p3,p6)*za(p3,p7)*zb(p1,p4)*zb(p3,p5)*zb(p7,p8)*s278**(-1)
     &     + 2*za(p2,p7)*za(p3,p6)*za(p6,p7)*zb(p1,p4)*zb(p5,p6)*zb(p7,
     &    p8)*s278**(-1) + 2*za(p2,p8)*za(p3,p6)*zb(p1,p8)*zb(p4,p8)*
     &    zba2(p5,p3,p6,p7)*s128**(-1) )
      unsym(2,2)= + PZ(qq,qe,vqR,veL,prop3456)*PW(s56)*s3456**(-1)*
     & s56**(-1)*s28**(-1)*s456**(-1) * ( 2*za(p1,p2)*za(p1,p3)*za(p4,
     &    p6)*zb(p1,p8)*zb(p4,p5)*zb(p4,p7)*s128**(-1) + 2*za(p1,p2)*
     &    za(p1,p3)*za(p5,p6)*zb(p1,p8)*zb(p4,p5)*zb(p5,p7)*s128**(-1)
     &     - 2*za(p1,p2)*za(p2,p3)*zb(p2,p8)*zb(p4,p5)*zba2(p7,p4,p5,p6
     &    )*s128**(-1) + 2*za(p1,p3)*za(p2,p7)*za(p4,p6)*zb(p4,p5)*zb(
     &    p4,p7)*zb(p7,p8)*s278**(-1) + 2*za(p1,p3)*za(p2,p7)*za(p5,p6)
     &    *zb(p4,p5)*zb(p5,p7)*zb(p7,p8)*s278**(-1) + 2*za(p1,p3)*za(p2
     &    ,p8)*za(p4,p6)*zb(p4,p5)*zb(p4,p8)*zb(p7,p8)*s278**(-1) + 2*
     &    za(p1,p3)*za(p2,p8)*za(p5,p6)*zb(p4,p5)*zb(p5,p8)*zb(p7,p8)*
     &    s278**(-1) )
      unsym(2,2) = unsym(2,2) + PZ(qq,qn,vqR,vnL,prop3456)*PW(s56)*
     & s3456**(-1)*s56**(-1)*s28**(-1)*s356**(-1) * (  - 2*za(p1,p2)*
     &    za(p1,p3)*za(p3,p6)*zb(p1,p8)*zb(p3,p5)*zb(p4,p7)*s128**(-1)
     &     + 2*za(p1,p2)*za(p1,p6)*za(p3,p6)*zb(p1,p8)*zb(p4,p7)*zb(p5,
     &    p6)*s128**(-1) - 2*za(p1,p2)*za(p3,p6)*zb(p2,p8)*zb(p4,p7)*
     &    zba2(p5,p3,p6,p2)*s128**(-1) - 2*za(p1,p3)*za(p2,p7)*za(p3,p6
     &    )*zb(p3,p5)*zb(p4,p7)*zb(p7,p8)*s278**(-1) - 2*za(p1,p3)*za(
     &    p2,p8)*za(p3,p6)*zb(p3,p5)*zb(p4,p8)*zb(p7,p8)*s278**(-1) + 2
     &    *za(p1,p6)*za(p2,p7)*za(p3,p6)*zb(p4,p7)*zb(p5,p6)*zb(p7,p8)*
     &    s278**(-1) + 2*za(p1,p6)*za(p2,p8)*za(p3,p6)*zb(p4,p8)*zb(p5,
     &    p6)*zb(p7,p8)*s278**(-1) )
      return
      end
