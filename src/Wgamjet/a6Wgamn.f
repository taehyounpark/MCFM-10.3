!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a6Wgamn(p1,p2,p3,p4,p5,p6,p,n,za,zb,b)
      implicit none
c---- Matrix element for Wgamma radiation from line 12
c  d(-p1) +ubar(-p2) --> e^-(p3)+ve~^+(p4)+gam(p5)+g(p6)

c                                                5     3-----<--4
c                                                \      /
c                                           gamma \    / W
c                                                  \  /
c            5         3-----<--4                   \/
c            |gamma        |W                        |W
c   2 ----<--|-------------|----1      2 ----<-------|-------------1
c                 0                               0
c                 0                               0
c                 0                               0
c              jtype=1                         jtype=3

c        3-----<-- 4       5
c            |W            |gamma
c   2 ----<--|-------------|----1
c                 0
c                 0
c                 0
c                jtype=2


c Overall factor of  4*gs^2*e^2*gw^4*(T^A)_{i2,i1}
c     Tr{T^A T^B} = delta^{AB}
c     including QED type propagators
      include 'types.f'
      include 'cplx.h'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'ewcharge.f'
      include 'masses.f'
      include 'nwz.f'
      integer::p1,p2,p3,p4,p5,p6,i,j
      real(dp):: p(mxpart,4),n(4),s3,s134,s234,s34,s345,s345ms34
      complex(dp):: iza,izb,zab2,b(2),prop34,prop345,vecm(mxpart,mxpart)

      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p1,p3)
      zab2(p1,p2,p3,p4)=+za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      iza(p1,p2)=cone/za(p1,p2)
      izb(p1,p2)=cone/zb(p1,p2)

      s134=s3(p1,p3,p4)
      s234=s3(p2,p3,p4)
      s34=s(p3,p4)
      s345=s3(p3,p4,p5)
      prop34=s34/cplx2(s34-wmass**2,wmass*wwidth)
      prop345=s345/cplx2(s345-wmass**2,wmass*wwidth)
      s345ms34=s345-s34

      do i=1,6
      do j=i,6
      call ndveccur(i,j,n,p,vecm)
      enddo
      enddo
c For nwz=+1 we call with za and zb interchanged so we need [a|n|b> = <b|n|a]
      if (nwz == +1) vecm=transpose(vecm)

      b(1)= + Q(1)*s345**(-1)*prop345 * (  - 2*za(p2,p3)*zb(p1,p3)*zb(
     &    p1,p4)*iza(p1,p6)*izb(p1,p5)*izb(p1,p6)*izb(p3,p5)*vecm(p1,p1
     &    ) + 2*za(p2,p3)*zb(p1,p3)*zb(p1,p4)*iza(p2,p6)*izb(p1,p5)*
     &    izb(p2,p6)*izb(p3,p5)*vecm(p2,p2) + 2*za(p2,p3)*zb(p1,p3)*zb(
     &    p4,p6)*iza(p1,p6)*izb(p1,p5)*izb(p1,p6)*izb(p3,p5)*vecm(p6,p1
     &    ) + za(p2,p3)*zb(p1,p4)*iza(p1,p6)*izb(p1,p5)*izb(p1,p6)*
     &    zab2(p5,p2,p6,p1)*vecm(p1,p1)*s345ms34**(-1) - za(p2,p3)*zb(
     &    p1,p4)*iza(p1,p6)*izb(p1,p5)*izb(p1,p6)*zab2(p5,p3,p4,p1)*
     &    vecm(p1,p1)*s345ms34**(-1) - za(p2,p3)*zb(p1,p4)*iza(p2,p6)*
     &    izb(p1,p5)*izb(p2,p6)*zab2(p5,p2,p6,p1)*vecm(p2,p2)*
     &    s345ms34**(-1) + za(p2,p3)*zb(p1,p4)*iza(p2,p6)*izb(p1,p5)*
     &    izb(p2,p6)*zab2(p5,p3,p4,p1)*vecm(p2,p2)*s345ms34**(-1) - za(
     &    p2,p3)*zb(p4,p6)*iza(p1,p6)*izb(p1,p5)*izb(p1,p6)*zab2(p5,p2,
     &    p6,p1)*vecm(p6,p1)*s345ms34**(-1) + za(p2,p3)*zb(p4,p6)*iza(
     &    p1,p6)*izb(p1,p5)*izb(p1,p6)*zab2(p5,p3,p4,p1)*vecm(p6,p1)*
     &    s345ms34**(-1) )
      b(1) = b(1) + Q(1)*s345**(-1)*prop345 * ( 2*za(p2,p5)*za(p3,p5)*
     &    zb(p4,p5)*iza(p1,p6)*izb(p1,p5)*vecm(p6,p1)*s345ms34**(-1) -
     &    2*za(p2,p5)*zb(p1,p4)*iza(p1,p6)*izb(p1,p6)*izb(p3,p5)*vecm(
     &    p1,p1) + 2*za(p2,p5)*zb(p1,p4)*iza(p2,p6)*izb(p2,p6)*izb(p3,
     &    p5)*vecm(p2,p2) + 2*za(p2,p5)*zb(p4,p6)*iza(p1,p6)*izb(p1,p6)
     &    *izb(p3,p5)*vecm(p6,p1) - za(p3,p5)*zb(p1,p4)*iza(p1,p6)*izb(
     &    p1,p5)*izb(p1,p6)*zab2(p2,p1,p6,p1)*vecm(p1,p1)*
     &    s345ms34**(-1) - za(p3,p5)*zb(p1,p4)*iza(p1,p6)*izb(p1,p5)*
     &    izb(p1,p6)*zab2(p2,p1,p6,p6)*vecm(p6,p1)*s345ms34**(-1) - 2*
     &    za(p3,p5)*zb(p1,p4)*iza(p1,p6)*izb(p1,p5)*izb(p1,p6)*zab2(p2,
     &    p3,p4,p1)*vecm(p1,p1)*s345ms34**(-1) - 2*za(p3,p5)*zb(p1,p4)*
     &    iza(p1,p6)*izb(p1,p5)*izb(p1,p6)*zab2(p2,p3,p4,p6)*vecm(p6,p1
     &    )*s345ms34**(-1) + za(p3,p5)*zb(p1,p4)*iza(p2,p6)*izb(p1,p5)*
     &    izb(p2,p6)*zab2(p2,p1,p6,p1)*vecm(p2,p2)*s345ms34**(-1) + 2*
     &    za(p3,p5)*zb(p1,p4)*iza(p2,p6)*izb(p1,p5)*izb(p2,p6)*zab2(p2,
     &    p3,p4,p1)*vecm(p2,p2)*s345ms34**(-1) )
      b(1) = b(1) + Q(1)*s345**(-1)*prop345 * ( za(p3,p5)*zb(p1,p4)*
     &    iza(p2,p6)*izb(p1,p5)*izb(p2,p6)*zab2(p6,p1,p2,p1)*vecm(p2,p6
     &    )*s345ms34**(-1) + 2*za(p3,p5)*zb(p1,p4)*iza(p2,p6)*izb(p1,p5
     &    )*izb(p2,p6)*zab2(p6,p3,p4,p1)*vecm(p2,p6)*s345ms34**(-1) - 2
     &    *za(p3,p6)*zb(p1,p3)*zb(p1,p4)*iza(p2,p6)*izb(p1,p5)*izb(p2,
     &    p6)*izb(p3,p5)*vecm(p2,p6) + za(p3,p6)*zb(p1,p4)*iza(p2,p6)*
     &    izb(p1,p5)*izb(p2,p6)*zab2(p5,p2,p6,p1)*vecm(p2,p6)*
     &    s345ms34**(-1) - za(p3,p6)*zb(p1,p4)*iza(p2,p6)*izb(p1,p5)*
     &    izb(p2,p6)*zab2(p5,p3,p4,p1)*vecm(p2,p6)*s345ms34**(-1) - 2*
     &    za(p5,p6)*zb(p1,p4)*iza(p2,p6)*izb(p2,p6)*izb(p3,p5)*vecm(p2,
     &    p6) )
      b(1) = b(1) + Q(1)*s34**(-1)*prop34 * (  - za(p2,p3)*zb(p1,p4)*
     &    iza(p1,p6)*izb(p1,p5)*izb(p1,p6)*zab2(p5,p2,p6,p1)*vecm(p1,p1
     &    )*s345ms34**(-1) + za(p2,p3)*zb(p1,p4)*iza(p1,p6)*izb(p1,p5)*
     &    izb(p1,p6)*zab2(p5,p3,p4,p1)*vecm(p1,p1)*s345ms34**(-1) + za(
     &    p2,p3)*zb(p1,p4)*iza(p2,p6)*izb(p1,p5)*izb(p2,p6)*zab2(p5,p2,
     &    p6,p1)*vecm(p2,p2)*s345ms34**(-1) - za(p2,p3)*zb(p1,p4)*iza(
     &    p2,p6)*izb(p1,p5)*izb(p2,p6)*zab2(p5,p3,p4,p1)*vecm(p2,p2)*
     &    s345ms34**(-1) + za(p2,p3)*zb(p4,p6)*iza(p1,p6)*izb(p1,p5)*
     &    izb(p1,p6)*zab2(p5,p2,p6,p1)*vecm(p6,p1)*s345ms34**(-1) - za(
     &    p2,p3)*zb(p4,p6)*iza(p1,p6)*izb(p1,p5)*izb(p1,p6)*zab2(p5,p3,
     &    p4,p1)*vecm(p6,p1)*s345ms34**(-1) - 2*za(p2,p3)*iza(p1,p6)*
     &    izb(p1,p5)*zab2(p5,p2,p3,p4)*vecm(p6,p1)*s234**(-1) - 2*za(p2
     &    ,p5)*za(p3,p5)*zb(p4,p5)*iza(p1,p6)*izb(p1,p5)*vecm(p6,p1)*
     &    s345ms34**(-1) + za(p3,p5)*zb(p1,p4)*iza(p1,p6)*izb(p1,p5)*
     &    izb(p1,p6)*zab2(p2,p1,p6,p1)*vecm(p1,p1)*s345ms34**(-1) + za(
     &    p3,p5)*zb(p1,p4)*iza(p1,p6)*izb(p1,p5)*izb(p1,p6)*zab2(p2,p1,
     &    p6,p6)*vecm(p6,p1)*s345ms34**(-1) )
      b(1) = b(1) + Q(1)*s34**(-1)*prop34 * ( 2*za(p3,p5)*zb(p1,p4)*
     &    iza(p1,p6)*izb(p1,p5)*izb(p1,p6)*zab2(p2,p3,p4,p1)*vecm(p1,p1
     &    )*s345ms34**(-1) + 2*za(p3,p5)*zb(p1,p4)*iza(p1,p6)*izb(p1,p5
     &    )*izb(p1,p6)*zab2(p2,p3,p4,p6)*vecm(p6,p1)*s345ms34**(-1) -
     &    za(p3,p5)*zb(p1,p4)*iza(p2,p6)*izb(p1,p5)*izb(p2,p6)*zab2(p2,
     &    p1,p6,p1)*vecm(p2,p2)*s345ms34**(-1) - 2*za(p3,p5)*zb(p1,p4)*
     &    iza(p2,p6)*izb(p1,p5)*izb(p2,p6)*zab2(p2,p3,p4,p1)*vecm(p2,p2
     &    )*s345ms34**(-1) - za(p3,p5)*zb(p1,p4)*iza(p2,p6)*izb(p1,p5)*
     &    izb(p2,p6)*zab2(p6,p1,p2,p1)*vecm(p2,p6)*s345ms34**(-1) - 2*
     &    za(p3,p5)*zb(p1,p4)*iza(p2,p6)*izb(p1,p5)*izb(p2,p6)*zab2(p6,
     &    p3,p4,p1)*vecm(p2,p6)*s345ms34**(-1) - za(p3,p6)*zb(p1,p4)*
     &    iza(p2,p6)*izb(p1,p5)*izb(p2,p6)*zab2(p5,p2,p6,p1)*vecm(p2,p6
     &    )*s345ms34**(-1) + za(p3,p6)*zb(p1,p4)*iza(p2,p6)*izb(p1,p5)*
     &    izb(p2,p6)*zab2(p5,p3,p4,p1)*vecm(p2,p6)*s345ms34**(-1) )
      b(1) = b(1) + Q(2)*s345**(-1)*prop345 * ( 2*za(p2,p3)*zb(p1,p3)*
     &    zb(p1,p4)*iza(p1,p6)*izb(p1,p5)*izb(p1,p6)*izb(p3,p5)*vecm(p1
     &    ,p1) - 2*za(p2,p3)*zb(p1,p3)*zb(p1,p4)*iza(p2,p6)*izb(p1,p5)*
     &    izb(p2,p6)*izb(p3,p5)*vecm(p2,p2) - 2*za(p2,p3)*zb(p1,p3)*zb(
     &    p4,p6)*iza(p1,p6)*izb(p1,p5)*izb(p1,p6)*izb(p3,p5)*vecm(p6,p1
     &    ) - za(p2,p3)*zb(p1,p4)*iza(p1,p6)*izb(p1,p5)*izb(p1,p6)*
     &    zab2(p5,p2,p6,p1)*vecm(p1,p1)*s345ms34**(-1) + za(p2,p3)*zb(
     &    p1,p4)*iza(p1,p6)*izb(p1,p5)*izb(p1,p6)*zab2(p5,p3,p4,p1)*
     &    vecm(p1,p1)*s345ms34**(-1) + za(p2,p3)*zb(p1,p4)*iza(p2,p6)*
     &    izb(p1,p5)*izb(p2,p6)*zab2(p5,p2,p6,p1)*vecm(p2,p2)*
     &    s345ms34**(-1) - za(p2,p3)*zb(p1,p4)*iza(p2,p6)*izb(p1,p5)*
     &    izb(p2,p6)*zab2(p5,p3,p4,p1)*vecm(p2,p2)*s345ms34**(-1) + za(
     &    p2,p3)*zb(p4,p6)*iza(p1,p6)*izb(p1,p5)*izb(p1,p6)*zab2(p5,p2,
     &    p6,p1)*vecm(p6,p1)*s345ms34**(-1) - za(p2,p3)*zb(p4,p6)*iza(
     &    p1,p6)*izb(p1,p5)*izb(p1,p6)*zab2(p5,p3,p4,p1)*vecm(p6,p1)*
     &    s345ms34**(-1) )
      b(1) = b(1) + Q(2)*s345**(-1)*prop345 * (  - 2*za(p2,p5)*za(p3,p5
     &    )*zb(p4,p5)*iza(p1,p6)*izb(p1,p5)*vecm(p6,p1)*s345ms34**(-1)
     &     + 2*za(p2,p5)*zb(p1,p4)*iza(p1,p6)*izb(p1,p6)*izb(p3,p5)*
     &    vecm(p1,p1) - 2*za(p2,p5)*zb(p1,p4)*iza(p2,p6)*izb(p2,p6)*
     &    izb(p3,p5)*vecm(p2,p2) - 2*za(p2,p5)*zb(p4,p6)*iza(p1,p6)*
     &    izb(p1,p6)*izb(p3,p5)*vecm(p6,p1) + za(p3,p5)*zb(p1,p4)*iza(
     &    p1,p6)*izb(p1,p5)*izb(p1,p6)*zab2(p2,p1,p6,p1)*vecm(p1,p1)*
     &    s345ms34**(-1) + za(p3,p5)*zb(p1,p4)*iza(p1,p6)*izb(p1,p5)*
     &    izb(p1,p6)*zab2(p2,p1,p6,p6)*vecm(p6,p1)*s345ms34**(-1) + 2*
     &    za(p3,p5)*zb(p1,p4)*iza(p1,p6)*izb(p1,p5)*izb(p1,p6)*zab2(p2,
     &    p3,p4,p1)*vecm(p1,p1)*s345ms34**(-1) + 2*za(p3,p5)*zb(p1,p4)*
     &    iza(p1,p6)*izb(p1,p5)*izb(p1,p6)*zab2(p2,p3,p4,p6)*vecm(p6,p1
     &    )*s345ms34**(-1) - za(p3,p5)*zb(p1,p4)*iza(p2,p6)*izb(p1,p5)*
     &    izb(p2,p6)*zab2(p2,p1,p6,p1)*vecm(p2,p2)*s345ms34**(-1) - 2*
     &    za(p3,p5)*zb(p1,p4)*iza(p2,p6)*izb(p1,p5)*izb(p2,p6)*zab2(p2,
     &    p3,p4,p1)*vecm(p2,p2)*s345ms34**(-1) )
      b(1) = b(1) + Q(2)*s345**(-1)*prop345 * (  - za(p3,p5)*zb(p1,p4)*
     &    iza(p2,p6)*izb(p1,p5)*izb(p2,p6)*zab2(p6,p1,p2,p1)*vecm(p2,p6
     &    )*s345ms34**(-1) - 2*za(p3,p5)*zb(p1,p4)*iza(p2,p6)*izb(p1,p5
     &    )*izb(p2,p6)*zab2(p6,p3,p4,p1)*vecm(p2,p6)*s345ms34**(-1) + 2
     &    *za(p3,p6)*zb(p1,p3)*zb(p1,p4)*iza(p2,p6)*izb(p1,p5)*izb(p2,
     &    p6)*izb(p3,p5)*vecm(p2,p6) - za(p3,p6)*zb(p1,p4)*iza(p2,p6)*
     &    izb(p1,p5)*izb(p2,p6)*zab2(p5,p2,p6,p1)*vecm(p2,p6)*
     &    s345ms34**(-1) + za(p3,p6)*zb(p1,p4)*iza(p2,p6)*izb(p1,p5)*
     &    izb(p2,p6)*zab2(p5,p3,p4,p1)*vecm(p2,p6)*s345ms34**(-1) + 2*
     &    za(p5,p6)*zb(p1,p4)*iza(p2,p6)*izb(p2,p6)*izb(p3,p5)*vecm(p2,
     &    p6) )
      b(1) = b(1) + Q(2)*s34**(-1)*prop34 * ( 2*za(p1,p3)*zb(p1,p2)*zb(
     &    p1,p4)*izb(p1,p5)*izb(p2,p5)*vecm(p2,p1)*s134**(-1) + 2*za(p1
     &    ,p3)*zb(p1,p4)*zb(p1,p5)*izb(p1,p5)*izb(p2,p5)*vecm(p5,p1)*
     &    s134**(-1) + za(p2,p3)*zb(p1,p4)*iza(p1,p6)*izb(p1,p5)*izb(p1
     &    ,p6)*zab2(p5,p2,p6,p1)*vecm(p1,p1)*s345ms34**(-1) - za(p2,p3)
     &    *zb(p1,p4)*iza(p1,p6)*izb(p1,p5)*izb(p1,p6)*zab2(p5,p3,p4,p1)
     &    *vecm(p1,p1)*s345ms34**(-1) - za(p2,p3)*zb(p1,p4)*iza(p2,p6)*
     &    izb(p1,p5)*izb(p2,p6)*zab2(p5,p2,p6,p1)*vecm(p2,p2)*
     &    s345ms34**(-1) + za(p2,p3)*zb(p1,p4)*iza(p2,p6)*izb(p1,p5)*
     &    izb(p2,p6)*zab2(p5,p3,p4,p1)*vecm(p2,p2)*s345ms34**(-1) - za(
     &    p2,p3)*zb(p4,p6)*iza(p1,p6)*izb(p1,p5)*izb(p1,p6)*zab2(p5,p2,
     &    p6,p1)*vecm(p6,p1)*s345ms34**(-1) + za(p2,p3)*zb(p4,p6)*iza(
     &    p1,p6)*izb(p1,p5)*izb(p1,p6)*zab2(p5,p3,p4,p1)*vecm(p6,p1)*
     &    s345ms34**(-1) + 2*za(p2,p5)*za(p3,p5)*zb(p4,p5)*iza(p1,p6)*
     &    izb(p1,p5)*vecm(p6,p1)*s345ms34**(-1) + 2*za(p2,p5)*zb(p1,p4)
     &    *iza(p2,p6)*izb(p1,p5)*izb(p2,p6)*zab2(p3,p1,p4,p1)*vecm(p2,
     &    p2)*s134**(-1) )
      b(1) = b(1) + Q(2)*s34**(-1)*prop34 * (  - za(p3,p5)*zb(p1,p4)*
     &    iza(p1,p6)*izb(p1,p5)*izb(p1,p6)*zab2(p2,p1,p6,p1)*vecm(p1,p1
     &    )*s345ms34**(-1) - za(p3,p5)*zb(p1,p4)*iza(p1,p6)*izb(p1,p5)*
     &    izb(p1,p6)*zab2(p2,p1,p6,p6)*vecm(p6,p1)*s345ms34**(-1) - 2*
     &    za(p3,p5)*zb(p1,p4)*iza(p1,p6)*izb(p1,p5)*izb(p1,p6)*zab2(p2,
     &    p3,p4,p1)*vecm(p1,p1)*s345ms34**(-1) - 2*za(p3,p5)*zb(p1,p4)*
     &    iza(p1,p6)*izb(p1,p5)*izb(p1,p6)*zab2(p2,p3,p4,p6)*vecm(p6,p1
     &    )*s345ms34**(-1) + za(p3,p5)*zb(p1,p4)*iza(p2,p6)*izb(p1,p5)*
     &    izb(p2,p6)*zab2(p2,p1,p6,p1)*vecm(p2,p2)*s345ms34**(-1) + 2*
     &    za(p3,p5)*zb(p1,p4)*iza(p2,p6)*izb(p1,p5)*izb(p2,p6)*zab2(p2,
     &    p3,p4,p1)*vecm(p2,p2)*s345ms34**(-1) + za(p3,p5)*zb(p1,p4)*
     &    iza(p2,p6)*izb(p1,p5)*izb(p2,p6)*zab2(p6,p1,p2,p1)*vecm(p2,p6
     &    )*s345ms34**(-1) + 2*za(p3,p5)*zb(p1,p4)*iza(p2,p6)*izb(p1,p5
     &    )*izb(p2,p6)*zab2(p6,p3,p4,p1)*vecm(p2,p6)*s345ms34**(-1) +
     &    za(p3,p6)*zb(p1,p4)*iza(p2,p6)*izb(p1,p5)*izb(p2,p6)*zab2(p5,
     &    p2,p6,p1)*vecm(p2,p6)*s345ms34**(-1) )
      b(1) = b(1) + Q(2)*s34**(-1)*prop34 * (  - za(p3,p6)*zb(p1,p4)*
     &    iza(p2,p6)*izb(p1,p5)*izb(p2,p6)*zab2(p5,p3,p4,p1)*vecm(p2,p6
     &    )*s345ms34**(-1) + 2*za(p4,p3)*zb(p1,p2)*zb(p1,p4)*izb(p1,p5)
     &    *izb(p2,p5)*vecm(p2,p4)*s134**(-1) + 2*za(p4,p3)*zb(p1,p4)*
     &    zb(p1,p5)*izb(p1,p5)*izb(p2,p5)*vecm(p5,p4)*s134**(-1) - 2*
     &    za(p5,p6)*zb(p1,p4)*iza(p2,p6)*izb(p1,p5)*izb(p2,p6)*zab2(p3,
     &    p1,p4,p1)*vecm(p2,p6)*s134**(-1) - 2*zb(p1,p4)*iza(p1,p6)*
     &    izb(p1,p5)*izb(p1,p6)*izb(p2,p5)*zab2(p3,p2,p5,p1)*vecm(p1,p1
     &    ) + 2*zb(p4,p6)*iza(p1,p6)*izb(p1,p5)*izb(p1,p6)*izb(p2,p5)*
     &    zab2(p3,p2,p5,p1)*vecm(p6,p1) )
      b(2)= + Q(1)*s345**(-1)*prop345 * ( 2*za(p2,p3)**2*zb(p1,p4)*iza(
     &    p1,p6)*iza(p2,p5)*iza(p3,p5)*izb(p1,p6)*vecm(p1,p1) - 2*za(p2
     &    ,p3)**2*zb(p1,p4)*iza(p2,p5)*iza(p2,p6)*iza(p3,p5)*izb(p2,p6)
     &    *vecm(p2,p2) - 2*za(p2,p3)**2*zb(p4,p6)*iza(p1,p6)*iza(p2,p5)
     &    *iza(p3,p5)*izb(p1,p6)*vecm(p6,p1) + 2*za(p2,p3)*za(p3,p6)*
     &    zb(p1,p4)*iza(p2,p5)*iza(p2,p6)*iza(p3,p5)*izb(p2,p6)*vecm(p2
     &    ,p6) - za(p2,p3)*zb(p1,p4)*iza(p1,p6)*iza(p2,p5)*izb(p1,p6)*
     &    zab2(p2,p1,p6,p5)*vecm(p1,p1)*s345ms34**(-1) + za(p2,p3)*zb(
     &    p1,p4)*iza(p1,p6)*iza(p2,p5)*izb(p1,p6)*zab2(p2,p3,p4,p5)*
     &    vecm(p1,p1)*s345ms34**(-1) + za(p2,p3)*zb(p1,p4)*iza(p2,p5)*
     &    iza(p2,p6)*izb(p2,p6)*zab2(p2,p1,p6,p5)*vecm(p2,p2)*
     &    s345ms34**(-1) - za(p2,p3)*zb(p1,p4)*iza(p2,p5)*iza(p2,p6)*
     &    izb(p2,p6)*zab2(p2,p3,p4,p5)*vecm(p2,p2)*s345ms34**(-1) + za(
     &    p2,p3)*zb(p4,p5)*iza(p1,p6)*iza(p2,p5)*izb(p1,p6)*zab2(p2,p1,
     &    p6,p1)*vecm(p1,p1)*s345ms34**(-1) + za(p2,p3)*zb(p4,p5)*iza(
     &    p1,p6)*iza(p2,p5)*izb(p1,p6)*zab2(p2,p1,p6,p6)*vecm(p6,p1)*
     &    s345ms34**(-1) )
      b(2) = b(2) + Q(1)*s345**(-1)*prop345 * ( 2*za(p2,p3)*zb(p4,p5)*
     &    iza(p1,p6)*iza(p2,p5)*izb(p1,p6)*zab2(p2,p3,p4,p1)*vecm(p1,p1
     &    )*s345ms34**(-1) + 2*za(p2,p3)*zb(p4,p5)*iza(p1,p6)*iza(p2,p5
     &    )*izb(p1,p6)*zab2(p2,p3,p4,p6)*vecm(p6,p1)*s345ms34**(-1) -
     &    za(p2,p3)*zb(p4,p5)*iza(p2,p5)*iza(p2,p6)*izb(p2,p6)*zab2(p2,
     &    p1,p6,p1)*vecm(p2,p2)*s345ms34**(-1) - 2*za(p2,p3)*zb(p4,p5)*
     &    iza(p2,p5)*iza(p2,p6)*izb(p2,p6)*zab2(p2,p3,p4,p1)*vecm(p2,p2
     &    )*s345ms34**(-1) - za(p2,p3)*zb(p4,p5)*iza(p2,p5)*iza(p2,p6)*
     &    izb(p2,p6)*zab2(p6,p1,p2,p1)*vecm(p2,p6)*s345ms34**(-1) - 2*
     &    za(p2,p3)*zb(p4,p5)*iza(p2,p5)*iza(p2,p6)*izb(p2,p6)*zab2(p6,
     &    p3,p4,p1)*vecm(p2,p6)*s345ms34**(-1) + za(p2,p3)*zb(p4,p6)*
     &    iza(p1,p6)*iza(p2,p5)*izb(p1,p6)*zab2(p2,p1,p6,p5)*vecm(p6,p1
     &    )*s345ms34**(-1) - za(p2,p3)*zb(p4,p6)*iza(p1,p6)*iza(p2,p5)*
     &    izb(p1,p6)*zab2(p2,p3,p4,p5)*vecm(p6,p1)*s345ms34**(-1) + 2*
     &    za(p3,p5)*zb(p1,p5)*zb(p4,p5)*iza(p2,p5)*izb(p2,p6)*vecm(p2,
     &    p6)*s345ms34**(-1) )
      b(2) = b(2) + Q(1)*s345**(-1)*prop345 * (  - za(p3,p6)*zb(p1,p4)*
     &    iza(p2,p5)*iza(p2,p6)*izb(p2,p6)*zab2(p2,p1,p6,p5)*vecm(p2,p6
     &    )*s345ms34**(-1) + za(p3,p6)*zb(p1,p4)*iza(p2,p5)*iza(p2,p6)*
     &    izb(p2,p6)*zab2(p2,p3,p4,p5)*vecm(p2,p6)*s345ms34**(-1) )
      b(2) = b(2) + Q(1)*s34**(-1)*prop34 * (  - 2*za(p1,p2)*za(p2,p3)*
     &    zb(p4,p2)*iza(p1,p5)*iza(p2,p5)*vecm(p2,p1)*s234**(-1) - 2*
     &    za(p1,p2)*za(p2,p3)*zb(p4,p3)*iza(p1,p5)*iza(p2,p5)*vecm(p3,
     &    p1)*s234**(-1) - 2*za(p2,p3)*za(p5,p2)*zb(p4,p2)*iza(p1,p5)*
     &    iza(p2,p5)*vecm(p2,p5)*s234**(-1) - 2*za(p2,p3)*za(p5,p2)*zb(
     &    p4,p3)*iza(p1,p5)*iza(p2,p5)*vecm(p3,p5)*s234**(-1) + za(p2,
     &    p3)*zb(p1,p4)*iza(p1,p6)*iza(p2,p5)*izb(p1,p6)*zab2(p2,p1,p6,
     &    p5)*vecm(p1,p1)*s345ms34**(-1) - za(p2,p3)*zb(p1,p4)*iza(p1,
     &    p6)*iza(p2,p5)*izb(p1,p6)*zab2(p2,p3,p4,p5)*vecm(p1,p1)*
     &    s345ms34**(-1) - za(p2,p3)*zb(p1,p4)*iza(p2,p5)*iza(p2,p6)*
     &    izb(p2,p6)*zab2(p2,p1,p6,p5)*vecm(p2,p2)*s345ms34**(-1) + za(
     &    p2,p3)*zb(p1,p4)*iza(p2,p5)*iza(p2,p6)*izb(p2,p6)*zab2(p2,p3,
     &    p4,p5)*vecm(p2,p2)*s345ms34**(-1) - 2*za(p2,p3)*zb(p1,p5)*
     &    iza(p1,p6)*iza(p2,p5)*izb(p1,p6)*zab2(p2,p2,p3,p4)*vecm(p1,p1
     &    )*s234**(-1) - za(p2,p3)*zb(p4,p5)*iza(p1,p6)*iza(p2,p5)*izb(
     &    p1,p6)*zab2(p2,p1,p6,p1)*vecm(p1,p1)*s345ms34**(-1) )
      b(2) = b(2) + Q(1)*s34**(-1)*prop34 * (  - za(p2,p3)*zb(p4,p5)*
     &    iza(p1,p6)*iza(p2,p5)*izb(p1,p6)*zab2(p2,p1,p6,p6)*vecm(p6,p1
     &    )*s345ms34**(-1) - 2*za(p2,p3)*zb(p4,p5)*iza(p1,p6)*iza(p2,p5
     &    )*izb(p1,p6)*zab2(p2,p3,p4,p1)*vecm(p1,p1)*s345ms34**(-1) - 2
     &    *za(p2,p3)*zb(p4,p5)*iza(p1,p6)*iza(p2,p5)*izb(p1,p6)*zab2(p2
     &    ,p3,p4,p6)*vecm(p6,p1)*s345ms34**(-1) + za(p2,p3)*zb(p4,p5)*
     &    iza(p2,p5)*iza(p2,p6)*izb(p2,p6)*zab2(p2,p1,p6,p1)*vecm(p2,p2
     &    )*s345ms34**(-1) + 2*za(p2,p3)*zb(p4,p5)*iza(p2,p5)*iza(p2,p6
     &    )*izb(p2,p6)*zab2(p2,p3,p4,p1)*vecm(p2,p2)*s345ms34**(-1) +
     &    za(p2,p3)*zb(p4,p5)*iza(p2,p5)*iza(p2,p6)*izb(p2,p6)*zab2(p6,
     &    p1,p2,p1)*vecm(p2,p6)*s345ms34**(-1) + 2*za(p2,p3)*zb(p4,p5)*
     &    iza(p2,p5)*iza(p2,p6)*izb(p2,p6)*zab2(p6,p3,p4,p1)*vecm(p2,p6
     &    )*s345ms34**(-1) - za(p2,p3)*zb(p4,p6)*iza(p1,p6)*iza(p2,p5)*
     &    izb(p1,p6)*zab2(p2,p1,p6,p5)*vecm(p6,p1)*s345ms34**(-1) + za(
     &    p2,p3)*zb(p4,p6)*iza(p1,p6)*iza(p2,p5)*izb(p1,p6)*zab2(p2,p3,
     &    p4,p5)*vecm(p6,p1)*s345ms34**(-1) )
      b(2) = b(2) + Q(1)*s34**(-1)*prop34 * ( 2*za(p2,p3)*zb(p5,p6)*
     &    iza(p1,p6)*iza(p2,p5)*izb(p1,p6)*zab2(p2,p2,p3,p4)*vecm(p6,p1
     &    )*s234**(-1) + 2*za(p2,p3)*iza(p1,p5)*iza(p2,p5)*iza(p2,p6)*
     &    izb(p2,p6)*zab2(p2,p1,p5,p4)*vecm(p2,p2) - 2*za(p3,p5)*zb(p1,
     &    p5)*zb(p4,p5)*iza(p2,p5)*izb(p2,p6)*vecm(p2,p6)*
     &    s345ms34**(-1) + za(p3,p6)*zb(p1,p4)*iza(p2,p5)*iza(p2,p6)*
     &    izb(p2,p6)*zab2(p2,p1,p6,p5)*vecm(p2,p6)*s345ms34**(-1) - za(
     &    p3,p6)*zb(p1,p4)*iza(p2,p5)*iza(p2,p6)*izb(p2,p6)*zab2(p2,p3,
     &    p4,p5)*vecm(p2,p6)*s345ms34**(-1) - 2*za(p3,p6)*iza(p1,p5)*
     &    iza(p2,p5)*iza(p2,p6)*izb(p2,p6)*zab2(p2,p1,p5,p4)*vecm(p2,p6
     &    ) )
      b(2) = b(2) + Q(2)*s345**(-1)*prop345 * (  - 2*za(p2,p3)**2*zb(p1
     &    ,p4)*iza(p1,p6)*iza(p2,p5)*iza(p3,p5)*izb(p1,p6)*vecm(p1,p1)
     &     + 2*za(p2,p3)**2*zb(p1,p4)*iza(p2,p5)*iza(p2,p6)*iza(p3,p5)*
     &    izb(p2,p6)*vecm(p2,p2) + 2*za(p2,p3)**2*zb(p4,p6)*iza(p1,p6)*
     &    iza(p2,p5)*iza(p3,p5)*izb(p1,p6)*vecm(p6,p1) - 2*za(p2,p3)*
     &    za(p3,p6)*zb(p1,p4)*iza(p2,p5)*iza(p2,p6)*iza(p3,p5)*izb(p2,
     &    p6)*vecm(p2,p6) + za(p2,p3)*zb(p1,p4)*iza(p1,p6)*iza(p2,p5)*
     &    izb(p1,p6)*zab2(p2,p1,p6,p5)*vecm(p1,p1)*s345ms34**(-1) - za(
     &    p2,p3)*zb(p1,p4)*iza(p1,p6)*iza(p2,p5)*izb(p1,p6)*zab2(p2,p3,
     &    p4,p5)*vecm(p1,p1)*s345ms34**(-1) - za(p2,p3)*zb(p1,p4)*iza(
     &    p2,p5)*iza(p2,p6)*izb(p2,p6)*zab2(p2,p1,p6,p5)*vecm(p2,p2)*
     &    s345ms34**(-1) + za(p2,p3)*zb(p1,p4)*iza(p2,p5)*iza(p2,p6)*
     &    izb(p2,p6)*zab2(p2,p3,p4,p5)*vecm(p2,p2)*s345ms34**(-1) - za(
     &    p2,p3)*zb(p4,p5)*iza(p1,p6)*iza(p2,p5)*izb(p1,p6)*zab2(p2,p1,
     &    p6,p1)*vecm(p1,p1)*s345ms34**(-1) - za(p2,p3)*zb(p4,p5)*iza(
     &    p1,p6)*iza(p2,p5)*izb(p1,p6)*zab2(p2,p1,p6,p6)*vecm(p6,p1)*
     &    s345ms34**(-1) )
      b(2) = b(2) + Q(2)*s345**(-1)*prop345 * (  - 2*za(p2,p3)*zb(p4,p5
     &    )*iza(p1,p6)*iza(p2,p5)*izb(p1,p6)*zab2(p2,p3,p4,p1)*vecm(p1,
     &    p1)*s345ms34**(-1) - 2*za(p2,p3)*zb(p4,p5)*iza(p1,p6)*iza(p2,
     &    p5)*izb(p1,p6)*zab2(p2,p3,p4,p6)*vecm(p6,p1)*s345ms34**(-1)
     &     + za(p2,p3)*zb(p4,p5)*iza(p2,p5)*iza(p2,p6)*izb(p2,p6)*zab2(
     &    p2,p1,p6,p1)*vecm(p2,p2)*s345ms34**(-1) + 2*za(p2,p3)*zb(p4,
     &    p5)*iza(p2,p5)*iza(p2,p6)*izb(p2,p6)*zab2(p2,p3,p4,p1)*vecm(
     &    p2,p2)*s345ms34**(-1) + za(p2,p3)*zb(p4,p5)*iza(p2,p5)*iza(p2
     &    ,p6)*izb(p2,p6)*zab2(p6,p1,p2,p1)*vecm(p2,p6)*s345ms34**(-1)
     &     + 2*za(p2,p3)*zb(p4,p5)*iza(p2,p5)*iza(p2,p6)*izb(p2,p6)*
     &    zab2(p6,p3,p4,p1)*vecm(p2,p6)*s345ms34**(-1) - za(p2,p3)*zb(
     &    p4,p6)*iza(p1,p6)*iza(p2,p5)*izb(p1,p6)*zab2(p2,p1,p6,p5)*
     &    vecm(p6,p1)*s345ms34**(-1) + za(p2,p3)*zb(p4,p6)*iza(p1,p6)*
     &    iza(p2,p5)*izb(p1,p6)*zab2(p2,p3,p4,p5)*vecm(p6,p1)*
     &    s345ms34**(-1) - 2*za(p3,p5)*zb(p1,p5)*zb(p4,p5)*iza(p2,p5)*
     &    izb(p2,p6)*vecm(p2,p6)*s345ms34**(-1) )
      b(2) = b(2) + Q(2)*s345**(-1)*prop345 * ( za(p3,p6)*zb(p1,p4)*
     &    iza(p2,p5)*iza(p2,p6)*izb(p2,p6)*zab2(p2,p1,p6,p5)*vecm(p2,p6
     &    )*s345ms34**(-1) - za(p3,p6)*zb(p1,p4)*iza(p2,p5)*iza(p2,p6)*
     &    izb(p2,p6)*zab2(p2,p3,p4,p5)*vecm(p2,p6)*s345ms34**(-1) )
      b(2) = b(2) + Q(2)*s34**(-1)*prop34 * (  - za(p2,p3)*zb(p1,p4)*
     &    iza(p1,p6)*iza(p2,p5)*izb(p1,p6)*zab2(p2,p1,p6,p5)*vecm(p1,p1
     &    )*s345ms34**(-1) + za(p2,p3)*zb(p1,p4)*iza(p1,p6)*iza(p2,p5)*
     &    izb(p1,p6)*zab2(p2,p3,p4,p5)*vecm(p1,p1)*s345ms34**(-1) + za(
     &    p2,p3)*zb(p1,p4)*iza(p2,p5)*iza(p2,p6)*izb(p2,p6)*zab2(p2,p1,
     &    p6,p5)*vecm(p2,p2)*s345ms34**(-1) - za(p2,p3)*zb(p1,p4)*iza(
     &    p2,p5)*iza(p2,p6)*izb(p2,p6)*zab2(p2,p3,p4,p5)*vecm(p2,p2)*
     &    s345ms34**(-1) + za(p2,p3)*zb(p4,p5)*iza(p1,p6)*iza(p2,p5)*
     &    izb(p1,p6)*zab2(p2,p1,p6,p1)*vecm(p1,p1)*s345ms34**(-1) + za(
     &    p2,p3)*zb(p4,p5)*iza(p1,p6)*iza(p2,p5)*izb(p1,p6)*zab2(p2,p1,
     &    p6,p6)*vecm(p6,p1)*s345ms34**(-1) + 2*za(p2,p3)*zb(p4,p5)*
     &    iza(p1,p6)*iza(p2,p5)*izb(p1,p6)*zab2(p2,p3,p4,p1)*vecm(p1,p1
     &    )*s345ms34**(-1) + 2*za(p2,p3)*zb(p4,p5)*iza(p1,p6)*iza(p2,p5
     &    )*izb(p1,p6)*zab2(p2,p3,p4,p6)*vecm(p6,p1)*s345ms34**(-1) -
     &    za(p2,p3)*zb(p4,p5)*iza(p2,p5)*iza(p2,p6)*izb(p2,p6)*zab2(p2,
     &    p1,p6,p1)*vecm(p2,p2)*s345ms34**(-1) )
      b(2) = b(2) + Q(2)*s34**(-1)*prop34 * (  - 2*za(p2,p3)*zb(p4,p5)*
     &    iza(p2,p5)*iza(p2,p6)*izb(p2,p6)*zab2(p2,p3,p4,p1)*vecm(p2,p2
     &    )*s345ms34**(-1) - za(p2,p3)*zb(p4,p5)*iza(p2,p5)*iza(p2,p6)*
     &    izb(p2,p6)*zab2(p6,p1,p2,p1)*vecm(p2,p6)*s345ms34**(-1) - 2*
     &    za(p2,p3)*zb(p4,p5)*iza(p2,p5)*iza(p2,p6)*izb(p2,p6)*zab2(p6,
     &    p3,p4,p1)*vecm(p2,p6)*s345ms34**(-1) + za(p2,p3)*zb(p4,p6)*
     &    iza(p1,p6)*iza(p2,p5)*izb(p1,p6)*zab2(p2,p1,p6,p5)*vecm(p6,p1
     &    )*s345ms34**(-1) - za(p2,p3)*zb(p4,p6)*iza(p1,p6)*iza(p2,p5)*
     &    izb(p1,p6)*zab2(p2,p3,p4,p5)*vecm(p6,p1)*s345ms34**(-1) + 2*
     &    za(p3,p5)*zb(p1,p5)*zb(p4,p5)*iza(p2,p5)*izb(p2,p6)*vecm(p2,
     &    p6)*s345ms34**(-1) - za(p3,p6)*zb(p1,p4)*iza(p2,p5)*iza(p2,p6
     &    )*izb(p2,p6)*zab2(p2,p1,p6,p5)*vecm(p2,p6)*s345ms34**(-1) +
     &    za(p3,p6)*zb(p1,p4)*iza(p2,p5)*iza(p2,p6)*izb(p2,p6)*zab2(p2,
     &    p3,p4,p5)*vecm(p2,p6)*s345ms34**(-1) + 2*zb(p1,p4)*iza(p2,p5)
     &    *izb(p2,p6)*zab2(p3,p1,p4,p5)*vecm(p2,p6)*s134**(-1) )

      return
      end
