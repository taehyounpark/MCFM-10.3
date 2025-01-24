!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a6Wgam(p1,p2,p3,p4,p5,p6,za,zb,b)
      implicit none
c---- Matrix element for Wgamma radiation from line 12
c  d(-p1) +ubar(-p2) --> e^-(p3)+ve~(p4)+gam(p5)+g(p6)

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
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'ewcharge.f'
      include 'masses.f'
      integer::p1,p2,p3,p4,p5,p6
      real(dp):: s3,s134,s234,s34,s345,Qdiff
      complex(dp):: iza,izb,zab2,b(2,2),prop34,prop345
c      amplitude b(h5,h6)

      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p1,p3)
      zab2(p1,p2,p3,p4)=+za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      iza(p1,p2)=cone/za(p1,p2)
      izb(p1,p2)=cone/zb(p1,p2)

      Qdiff=Q(2)-Q(1)
      s134=s3(p1,p3,p4)
      s234=s3(p2,p3,p4)
      s34=s(p3,p4)
      s345=s3(p3,p4,p5)
      prop34=cmplx(s34-wmass**2,wmass*wwidth,kind=dp)
      prop345=cmplx(s345-wmass**2,wmass*wwidth,kind=dp)

      b(1,1)= + Qdiff*prop34**(-1)*prop345**(-1) * (  - 1.D0/2.D0*za(p3
     &    ,p5)*zb(p1,p2)*zb(p1,p4)*izb(p1,p5)*izb(p1,p6)*izb(p2,p6)*
     &    zab2(p2,p1,p6,p1) + za(p3,p5)*zb(p1,p3)*zb(p1,p4)*izb(p1,p5)*
     &    izb(p1,p6)*izb(p2,p6)*zab2(p3,p2,p6,p1) + za(p3,p5)*zb(p1,p4)
     &    **2*izb(p1,p5)*izb(p1,p6)*izb(p2,p6)*zab2(p4,p2,p6,p1) - 1.D0/
     &    2.D0*za(p3,p5)*zb(p1,p4)*izb(p1,p5)*izb(p2,p6)*zab2(p6,p1,p2,
     &    p1) + 1.D0/2.D0*zb(p1,p4)*izb(p1,p5)*izb(p1,p6)*izb(p2,p6)*
     &    zab2(p3,p2,p6,p1)*zab2(p5,p2,p6,p1) - 1.D0/2.D0*zb(p1,p4)*
     &    izb(p1,p5)*izb(p1,p6)*izb(p2,p6)*zab2(p3,p2,p6,p1)*zab2(p5,p3
     &    ,p4,p1) )
      b(1,1) = b(1,1) + Qdiff*prop345**(-1) * ( zb(p1,p3)*zb(p1,p4)*
     &    izb(p1,p5)*izb(p1,p6)*izb(p2,p6)*izb(p3,p5)*zab2(p3,p2,p6,p1)
     &     + zb(p1,p4)*izb(p1,p6)*izb(p2,p6)*izb(p3,p5)*zab2(p5,p2,p6,
     &    p1) )
      b(1,1) = b(1,1) + Q(2)*prop34**(-1) * (  - zb(p1,p4)*izb(p1,p5)*
     &    izb(p1,p6)*izb(p2,p5)*zab2(p3,p1,p4,p1)*zab2(p6,p2,p5,p1)*
     &    s134**(-1) - zb(p1,p4)*izb(p1,p5)*izb(p1,p6)*izb(p2,p6)*zab2(
     &    p3,p1,p4,p1)*zab2(p5,p2,p6,p1)*s134**(-1) )
      b(1,2)= + Qdiff*prop34**(-1)*prop345**(-1) * ( 1.D0/2.D0*za(p1,p2
     &    )*za(p3,p5)*zb(p1,p4)*zb(p1,p6)*iza(p1,p6)*izb(p1,p5) + 1.D0/
     &    2.D0*za(p1,p2)*za(p3,p5)*zb(p1,p4)*iza(p1,p6)*iza(p2,p6)*izb(
     &    p1,p5)*zab2(p2,p1,p6,p1) + za(p2,p3)*za(p3,p5)*zb(p1,p4)*iza(
     &    p1,p6)*iza(p2,p6)*izb(p1,p5)*zab2(p2,p1,p6,p3) + 1.D0/2.D0*
     &    za(p2,p3)*iza(p1,p6)*iza(p2,p6)*izb(p1,p5)*zab2(p2,p1,p6,p4)*
     &    zab2(p5,p2,p6,p1) - 1.D0/2.D0*za(p2,p3)*iza(p1,p6)*iza(p2,p6)
     &    *izb(p1,p5)*zab2(p2,p1,p6,p4)*zab2(p5,p3,p4,p1) + za(p2,p4)*
     &    za(p3,p5)*zb(p1,p4)*iza(p1,p6)*iza(p2,p6)*izb(p1,p5)*zab2(p2,
     &    p1,p6,p4) - za(p2,p5)*za(p3,p5)*zb(p4,p5)*iza(p1,p6)*iza(p2,
     &    p6)*izb(p1,p5)*zab2(p2,p1,p6,p1) )
      b(1,2) = b(1,2) + Qdiff*prop345**(-1) * (  - iza(p1,p6)*iza(p2,p6
     &    )*izb(p1,p5)*izb(p3,p5)*zab2(p2,p1,p6,p4)*zab2(p2,p3,p5,p1) )
      b(1,2) = b(1,2) + Q(1)*prop34**(-1) * ( za(p2,p3)*iza(p1,p6)*iza(
     &    p2,p6)*izb(p1,p5)*zab2(p2,p1,p6,p1)*zab2(p5,p2,p3,p4)*
     &    s234**(-1) )
      b(1,2) = b(1,2) + Q(2)*prop34**(-1) * (  - za(p2,p5)*zb(p1,p4)*
     &    iza(p2,p6)*izb(p2,p5)*zab2(p3,p1,p4,p6)*s134**(-1) - iza(p1,
     &    p6)*iza(p2,p6)*izb(p1,p5)*izb(p2,p5)*zab2(p2,p1,p6,p4)*zab2(
     &    p3,p2,p5,p1) )
      b(2,1)= + Qdiff*prop34**(-1)*prop345**(-1) * ( 1.D0/2.D0*za(p2,p3
     &    )*zb(p1,p2)*zb(p4,p5)*iza(p2,p5)*izb(p1,p6)*izb(p2,p6)*zab2(
     &    p2,p1,p6,p1) - za(p2,p3)*zb(p1,p4)*zb(p4,p5)*iza(p2,p5)*izb(
     &    p1,p6)*izb(p2,p6)*zab2(p4,p2,p6,p1) + 1.D0/2.D0*za(p2,p3)*zb(
     &    p4,p5)*iza(p2,p5)*izb(p2,p6)*zab2(p6,p1,p2,p1) - za(p2,p6)*
     &    za(p3,p5)*zb(p1,p5)*zb(p4,p5)*iza(p2,p5)*izb(p2,p6) - 1.D0/2.D
     &    0*zb(p1,p4)*iza(p2,p5)*izb(p1,p6)*izb(p2,p6)*zab2(p2,p1,p6,p5
     &    )*zab2(p3,p2,p6,p1) + 1.D0/2.D0*zb(p1,p4)*iza(p2,p5)*izb(p1,
     &    p6)*izb(p2,p6)*zab2(p2,p3,p4,p5)*zab2(p3,p2,p6,p1) + zb(p1,p5
     &    )*zb(p4,p5)*izb(p1,p6)*izb(p2,p6)*zab2(p3,p2,p6,p1) + zb(p4,
     &    p5)*iza(p2,p5)*izb(p1,p6)*izb(p2,p6)*zab2(p2,p3,p5,p1)*zab2(
     &    p3,p2,p6,p1) )
      b(2,1) = b(2,1) + Qdiff*prop345**(-1) * (  - za(p2,p3)*zb(p1,p4)*
     &    iza(p2,p5)*iza(p3,p5)*izb(p1,p6)*izb(p2,p6)*zab2(p3,p2,p6,p1)
     &     )
      b(2,1) = b(2,1) + Q(1)*prop34**(-1) * ( za(p2,p3)*iza(p1,p5)*iza(
     &    p2,p5)*izb(p1,p6)*zab2(p2,p1,p5,p1)*zab2(p6,p2,p3,p4)*
     &    s234**(-1) - iza(p1,p5)*iza(p2,p5)*izb(p1,p6)*izb(p2,p6)*
     &    zab2(p2,p1,p5,p4)*zab2(p3,p2,p6,p1) )
      b(2,1) = b(2,1) + Q(2)*prop34**(-1) * (  - za(p2,p6)*zb(p1,p4)*
     &    iza(p2,p5)*izb(p2,p6)*zab2(p3,p1,p4,p5)*s134**(-1) )
      b(2,2)= + Qdiff*prop34**(-1)*prop345**(-1) * (  - 1.D0/2.D0*za(p1
     &    ,p2)*za(p2,p3)*zb(p1,p6)*zb(p4,p5)*iza(p1,p6)*iza(p2,p5) - 1.D
     &    0/2.D0*za(p1,p2)*za(p2,p3)*zb(p4,p5)*iza(p1,p6)*iza(p2,p5)*
     &    iza(p2,p6)*zab2(p2,p1,p6,p1) - za(p2,p3)**2*zb(p4,p5)*iza(p1,
     &    p6)*iza(p2,p5)*iza(p2,p6)*zab2(p2,p1,p6,p3) - za(p2,p3)*za(p2
     &    ,p4)*zb(p4,p5)*iza(p1,p6)*iza(p2,p5)*iza(p2,p6)*zab2(p2,p1,p6
     &    ,p4) - 1.D0/2.D0*za(p2,p3)*iza(p1,p6)*iza(p2,p5)*iza(p2,p6)*
     &    zab2(p2,p1,p6,p4)*zab2(p2,p1,p6,p5) + 1.D0/2.D0*za(p2,p3)*
     &    iza(p1,p6)*iza(p2,p5)*iza(p2,p6)*zab2(p2,p1,p6,p4)*zab2(p2,p3
     &    ,p4,p5) )
      b(2,2) = b(2,2) + Qdiff*prop345**(-1) * (  - za(p2,p3)**2*iza(p1,
     &    p6)*iza(p2,p5)*iza(p2,p6)*iza(p3,p5)*zab2(p2,p1,p6,p4) )
      b(2,2) = b(2,2) + Q(1)*prop34**(-1) * (  - za(p2,p3)**2*zb(p3,p4)
     &    *iza(p1,p5)*iza(p2,p5)*iza(p2,p6)*zab2(p2,p1,p5,p6)*
     &    s234**(-1) - za(p2,p3)**2*zb(p3,p4)*iza(p1,p6)*iza(p2,p5)*
     &    iza(p2,p6)*zab2(p2,p1,p6,p5)*s234**(-1) )

      return
      end
