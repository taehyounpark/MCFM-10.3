!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a6WpgamxnWrad(p1,p2,p3,p4,p5,p6,za,zb,b22,b21)
      implicit none
c---- Matrix element for Wgamma radiation from W decay
c  u(-p1) +dbar(-p2) --> ve(p3)+e^+(p4)+gam(p5)+g(p6)

c               5     3-----<--4
c               \      /
c           gam \    / W
c                 \  /
c                  \/
c                     |W
c     2 ----<-------|-------------1
c                0
c                0
c                0
c             jtype=3

      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'ewcharge.f'
      include 'masses.f'
      integer::p1,p2,p3,p4,p5,p6
      real(dp):: s3,s16,s26,s34,s345,Qdiff
      complex(dp):: iza,izb,zab2,b22,b21,
     & Lsm1,L0,L1,prop34,prop345
c      amplitude b(h5,h6)

      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p1,p3)
      zab2(p1,p2,p3,p4)=+za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      iza(p1,p2)=cone/za(p1,p2)
      izb(p1,p2)=cone/zb(p1,p2)

      Qdiff=Q(2)-Q(1)
      s16=s(p1,p6)
      s26=s(p2,p6)
      s34=s(p3,p4)
      s345=s3(p3,p4,p5)
      prop34=cmplx(s34-wmass**2,wmass*wwidth,kind=dp)
      prop345=cmplx(s345-wmass**2,wmass*wwidth,kind=dp)

      b22=Qdiff/prop345*iza(p2,p5)*iza(p1,p6)*iza(p2,p6)*(
     & -(za(p2,p1)*za(p2,p3)*L0(-s26,-s345)
     & *(zb(p1,p5)*zb(p4,p3)*za(p2,p3)
     & +0.5_dp*zb(p1,p6)*zb(p4,p5)*za(p2,p6)*(1._dp+s(p3,p4)/s345))
     & -zb(p4,p3)*zab2(p2,p3,p4,p5)*za(p2,p3)**2
     & *Lsm1(-s16,-s345,-s26,-s345)
     & +0.5_dp*za(p2,p1)**2*zab2(p3,p2,p6,p1)*L1(-s26,-s345)/s345
     & *(zb(p1,p5)*zb(p4,p3)*za(p2,p3)+zb(p4,p5)*zb(p1,p6)*za(p2,p6)))
     & /prop34)

      b21=Qdiff/prop345*iza(p2,p5)*izb(p1,p6)*izb(p2,p6)*(
     &  +zb(p1,p2)*za(p2,p3)*L0(-s16,-s345)/prop34
     &  *(za(p2,p3)*zb(p3,p4)*zb(p1,p5)
     &  -0.5_dp*za(p2,p6)*zb(p1,p6)*zb(p4,p5)*(1._dp+s(p3,p4)/s345))
     &  +zab2(p3,p2,p6,p1)*Lsm1(-s26,-s345,-s16,-s345)/prop34
     &  *(zb(p1,p4)*zab2(p2,p3,p4,p5)+za(p2,p5)*zb(p1,p5)*zb(p4,p5))
     &  +0.5_dp*(za(p2,p3)*zb(p1,p2))**2*zb(p3,p4)*zab2(p2,p1,p6,p5)
     &  *L1(-s16,-s345)/s345/prop34)

      return
      end
