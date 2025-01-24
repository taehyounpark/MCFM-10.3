!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a6Wpgam(p1,p2,p3,p4,p5,p6,za,zb,b)
      implicit none
c---- Matrix element for W^+gamma radiation from line 12
c  u(-p1) +dbar(-p2) --> ve(p3)+e^+(p4)+gam(p5)+g(p6)

c                                                5     3-----<--4
c                                                \      /
c                                           gam \    / W
c                                                  \  /
c            5         3-----<--4                   \/
c            |gam            |W                        |W
c   2 ----<--|-------------|----1      2 ----<-------|-------------1
c                 0                               0
c                 0                               0
c                 0                               0
c              jtype=1                         jtype=3

c        3-----<-- 4       5
c            |W            |gam
c   2 ----<--|-------------|----1
c                 0
c                 0
c                 0
c                jtype=2


c Overall factor of  i_*gs^2*e^2*gw^4/8*(T^A)_{i2,i1}
c     Tr{T^A T^B} = delta^{AB}
c     including QED type propagators
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'ewcharge.f'
      include 'masses.f'
      integer::p1,p2,p3,p4,p5,p6
      real(dp):: s3,s134,s234,s34,s345,Qdiff
      complex(dp):: zab2,b(2,2),prop34,prop345
c      amplitude b(h5,h6)

      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p1,p3)
      zab2(p1,p2,p3,p4)=+za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      Qdiff=Q(2)-Q(1)
      s134=s3(p1,p3,p4)
      s234=s3(p2,p3,p4)
      s34=s(p3,p4)
      s345=s3(p3,p4,p5)
      prop34=cmplx(s34-wmass**2,wmass*wwidth,kind=dp)
      prop345=cmplx(s345-wmass**2,wmass*wwidth,kind=dp)

      b(2,2)=(Q(2)*za(p2,p3)**2
     & *(za(p1,p5)*zab2(p2,p4,p3,p5)-za(p1,p6)*zab2(p2,p1,p5,p6))
     & *zb(p3,p4))
     & /(za(p1,p5)*za(p1,p6)*za(p2,p5)*za(p2,p6)*prop34*s234)
     & +(Qdiff*za(p2,p3)**2
     & *(zab2(p2,p1,p6,p3)*zb(p4,p5)
     & -zab2(p2,p1,p6,p4)*zb(p3,p5)))
     & /(za(p1,p6)*za(p2,p5)*za(p2,p6)*prop34*prop345)
     & -(Qdiff*za(p2,p3)
     & *(zab2(p2,p4,p3,p5)*za(p2,p5)
     & -zab2(p2,p1,p6,p4)*za(p2,p4)))
     & /(za(p1,p6)*za(p2,p5)*za(p2,p6)*za(p4,p5)*prop345)

      b(2,1)=(-(Q(2)*(zab2(p2,p1,p5,p4)*zab2(p3,p2,p6,p1)*s234
     & +zb(p1,p5)*za(p2,p3)*za(p2,p5)*zb(p2,p6)
     & *zab2(p6,p2,p3,p4)))
     & /(za(p1,p5)*zb(p1,p6)*za(p2,p5)*zb(p2,p6)*prop34*s234))
     & -(Q(1)*zb(p1,p4)*za(p2,p6)*zab2(p3,p1,p4,p5))
     & /(za(p2,p5)*zb(p2,p6)*prop34*s134)
     & +(Qdiff*(zb(p1,p4)*za(p2,p3)*zab2(p4,p2,p6,p1)*zb(p4,p5)
     & +zb(p1,p5)*zb(p1,p6)*za(p2,p6)*za(p3,p5)*zb(p4,p5)
     & +zb(p1,p3)*za(p2,p3)*zab2(p3,p2,p6,p1)*zb(p4,p5)
     & -zb(p1,p4)*zab2(p2,p4,p3,p5)*zab2(p3,p2,p6,p1)))
     & /(zb(p1,p6)*za(p2,p5)*zb(p2,p6)*prop34*prop345)
     & -(Qdiff*zab2(p2,p4,p5,p1)*zab2(p3,p2,p6,p1))
     & /(zb(p1,p6)*za(p2,p5)*zb(p2,p6)*za(p4,p5)*prop345)

      b(1,2)=(-(Q(2)*zb(p1,p6)*za(p2,p3)*zab2(p5,p2,p3,p4))
     & /(zb(p1,p5)*za(p1,p6)*prop34*s234))
     & -(Q(1)*(zab2(p2,p1,p6,p4)*zab2(p3,p2,p5,p1)*s134
     & +zb(p1,p4)*zb(p1,p5)*za(p1,p6)*za(p2,p5)
     & *zab2(p3,p1,p4,p6)))
     & /(zb(p1,p5)*za(p1,p6)*zb(p2,p5)*za(p2,p6)*prop34*s134)
     & +(Qdiff*(zab2(p2,p1,p6,p4)*za(p2,p3)*zab2(p5,p4,p3,p1)
     & -zb(p1,p6)*za(p2,p5)*za(p2,p6)*za(p3,p5)*zb(p4,p5)
     & -zb(p1,p4)*zab2(p2,p1,p6,p4)*za(p2,p4)*za(p3,p5)
     & -zb(p1,p4)*zab2(p2,p1,p6,p3)*za(p2,p3)*za(p3,p5)))
     & /(zb(p1,p5)*za(p1,p6)*za(p2,p6)*prop34*prop345)
     & -(Qdiff*zb(p1,p4)*zab2(p2,p1,p6,p4)*za(p2,p3))
     & /(zb(p1,p5)*za(p1,p6)*za(p2,p6)*zb(p4,p5)*prop345)

      b(1,1)=(Q(1)*zb(p1,p4)**2*za(p3,p4)
     & *(zb(p2,p6)*zab2(p6,p2,p5,p1)
     & -zb(p2,p5)*zab2(p5,p4,p3,p1)))
     & /(zb(p1,p5)*zb(p1,p6)*zb(p2,p5)*zb(p2,p6)*prop34*s134)
     & +(Qdiff*zb(p1,p4)**2
     & *(zab2(p3,p2,p6,p1)*za(p4,p5)
     & -za(p3,p5)*zab2(p4,p2,p6,p1)))
     & /(zb(p1,p5)*zb(p1,p6)*zb(p2,p6)*prop34*prop345)
     & -(Qdiff*zb(p1,p4)**2*zab2(p3,p2,p6,p1))
     & /(zb(p1,p5)*zb(p1,p6)*zb(p2,p6)*zb(p4,p5)*prop345)

      return
      end
