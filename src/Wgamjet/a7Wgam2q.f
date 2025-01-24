!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a7Wgam2q(p1,p2,p3,p4,p5,p6,p7,za,zb,a,b)
      implicit none
c---- Matrix element for Wgam radiation from line 16
c  d(-p1) +c(-p2) --> e^-(p3)+ve~^+(p4)+gamma(p5)+u(p6)+c(p7)

c                                                5     3-----<--4
c                                                \      /
c                                             gam \    / W
c                                                  \  /
c            5         3-----<--4                   \/
c            |gam          |W                        |W
c   6 ----<--|-------------|----1      6 ----<-------|-------------1
c                 0                               0
c                 0                               0
c                 0                               0
c   7 -----<--------------------2      7 -----<--------------------2
c              jtype=1                         jtype=3

c        3-----<-- 4       5
c            |W            |gam
c   6 ----<--|-------------|----1
c                 0
c                 0
c                 0
c   7 -----<--------------------2
c                jtype=2
c  + Singly resonant diagram

c Overall factor of  i_*gs^2*e*rt2*gw^2*(T^A)_{i7,i1}*(T^A)_{i8,i2}
c     Tr{T^A T^B} = delta^{AB}
c     including QED type propagators
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'nf.f'
      include 'ewcharge.f'
      integer::p1,p2,p3,p4,p5,p6,p7
      real(dp):: s3,s134,s345,s257,s267,s34,s346,s127,s27
      complex(dp):: zab2,iza,izb,a(2,2),b(2,2),prop34,prop345,
     & zaa22,zbb22
c     amplitude a(h5,h27)
      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p1,p3)
      iza(p1,p2)=cone/za(p1,p2)
      izb(p1,p2)=cone/zb(p1,p2)
      zab2(p1,p2,p3,p4)=+za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zaa22(p1,p2,p3,p4,p5,p6)=
     & +zab2(p1,p2,p3,p4)*za(p4,p6)+zab2(p1,p2,p3,p5)*za(p5,p6)
      zbb22(p1,p2,p3,p4,p5,p6)=
     & +zb(p1,p2)*zab2(p2,p4,p5,p6)+zb(p1,p3)*zab2(p3,p4,p5,p6)

      s127=s3(p1,p2,p7)
      s134=s3(p1,p3,p4)
      s267=s3(p2,p6,p7)
      s345=s3(p3,p4,p5)
      s346=s3(p3,p4,p6)
      s257=s3(p2,p5,p7)
      s34=s(p3,p4)
      s27=s(p2,p7)

      prop34=cmplx(s34-wmass**2,wmass*wwidth,kind=dp)
      prop345=cmplx(s345-wmass**2,wmass*wwidth,kind=dp)
      a(1,1)= + Q(1)*s27**(-1) * ( za(p3,p5)*za(p3,p6)*zb(p1,p2)*zb(p1,
     &    p4)*izb(p1,p5)*zab2(p7,p1,p2,p3)*prop34**(-1)*prop345**(-1)*
     &    s127**(-1) + za(p3,p5)*za(p4,p6)*zb(p1,p2)*zb(p1,p4)*izb(p1,
     &    p5)*zab2(p7,p1,p2,p4)*prop34**(-1)*prop345**(-1)*s127**(-1)
     &     - za(p3,p5)*za(p5,p6)*zb(p1,p2)*zb(p1,p4)*izb(p1,p5)*zab2(p7
     &    ,p1,p2,p5)*prop34**(-1)*prop345**(-1)*s127**(-1) - 2*za(p3,p5
     &    )*za(p5,p6)*zb(p1,p2)*zb(p4,p5)*izb(p1,p5)*zab2(p7,p1,p2,p1)*
     &    prop34**(-1)*prop345**(-1)*s127**(-1) + za(p3,p5)*za(p6,p7)*
     &    zb(p1,p3)*zb(p1,p4)*izb(p1,p5)*zab2(p3,p6,p7,p2)*prop34**(-1)
     &    *prop345**(-1)*s267**(-1) + za(p3,p5)*za(p6,p7)*zb(p1,p4)**2*
     &    izb(p1,p5)*zab2(p4,p6,p7,p2)*prop34**(-1)*prop345**(-1)*
     &    s267**(-1) - za(p3,p5)*za(p6,p7)*zb(p1,p4)*zab2(p5,p6,p7,p2)*
     &    prop34**(-1)*prop345**(-1)*s267**(-1) + 2*za(p3,p6)*zb(p1,p2)
     &    *zb(p1,p3)*izb(p1,p5)*izb(p3,p5)*zab2(p7,p1,p2,p4)*
     &    prop345**(-1)*s127**(-1) - 2*za(p3,p6)*zb(p1,p2)*izb(p1,p5)*
     &    zab2(p5,p3,p4,p1)*zab2(p7,p1,p2,p4)*prop34**(-1)*
     &    prop345**(-1)*s127**(-1) )
      a(1,1) = a(1,1) + Q(1)*s27**(-1) * (  - 2*za(p3,p6)*zb(p1,p2)*
     &    izb(p1,p5)*zab2(p5,p3,p6,p4)*zab2(p7,p1,p2,p1)*prop34**(-1)*
     &    s127**(-1)*s346**(-1) + 2*za(p5,p6)*zb(p1,p2)*izb(p3,p5)*
     &    zab2(p7,p1,p2,p4)*prop345**(-1)*s127**(-1) + 2*za(p6,p7)*zb(
     &    p1,p3)*zb(p1,p4)*izb(p1,p5)*izb(p3,p5)*zab2(p3,p6,p7,p2)*
     &    prop345**(-1)*s267**(-1) - 2*za(p6,p7)*zb(p1,p4)*izb(p1,p5)*
     &    zab2(p3,p6,p7,p2)*zab2(p5,p3,p4,p1)*prop34**(-1)*
     &    prop345**(-1)*s267**(-1) + 2*za(p6,p7)*zb(p1,p4)*izb(p3,p5)*
     &    zab2(p5,p6,p7,p2)*prop345**(-1)*s267**(-1) )
      a(1,1) = a(1,1) + Q(2)*s27**(-1) * (  - za(p3,p5)*za(p3,p6)*zb(p1
     &    ,p2)*zb(p1,p4)*izb(p1,p5)*zab2(p7,p1,p2,p3)*prop34**(-1)*
     &    prop345**(-1)*s127**(-1) - za(p3,p5)*za(p4,p6)*zb(p1,p2)*zb(
     &    p1,p4)*izb(p1,p5)*zab2(p7,p1,p2,p4)*prop34**(-1)*
     &    prop345**(-1)*s127**(-1) + za(p3,p5)*za(p5,p6)*zb(p1,p2)*zb(
     &    p1,p4)*izb(p1,p5)*zab2(p7,p1,p2,p5)*prop34**(-1)*
     &    prop345**(-1)*s127**(-1) + 2*za(p3,p5)*za(p5,p6)*zb(p1,p2)*
     &    zb(p4,p5)*izb(p1,p5)*zab2(p7,p1,p2,p1)*prop34**(-1)*
     &    prop345**(-1)*s127**(-1) - za(p3,p5)*za(p6,p7)*zb(p1,p3)*zb(
     &    p1,p4)*izb(p1,p5)*zab2(p3,p6,p7,p2)*prop34**(-1)*
     &    prop345**(-1)*s267**(-1) - za(p3,p5)*za(p6,p7)*zb(p1,p4)**2*
     &    izb(p1,p5)*zab2(p4,p6,p7,p2)*prop34**(-1)*prop345**(-1)*
     &    s267**(-1) + za(p3,p5)*za(p6,p7)*zb(p1,p4)*zab2(p5,p6,p7,p2)*
     &    prop34**(-1)*prop345**(-1)*s267**(-1) - 2*za(p3,p6)*zb(p1,p2)
     &    *zb(p1,p3)*izb(p1,p5)*izb(p3,p5)*zab2(p7,p1,p2,p4)*
     &    prop345**(-1)*s127**(-1) )
      a(1,1) = a(1,1) + Q(2)*s27**(-1) * ( 2*za(p3,p6)*zb(p1,p2)*izb(p1
     &    ,p5)*zab2(p5,p3,p4,p1)*zab2(p7,p1,p2,p4)*prop34**(-1)*
     &    prop345**(-1)*s127**(-1) - 2*za(p5,p6)*zb(p1,p2)*izb(p3,p5)*
     &    zab2(p7,p1,p2,p4)*prop345**(-1)*s127**(-1) - 2*za(p6,p7)*zb(
     &    p1,p3)*zb(p1,p4)*izb(p1,p5)*izb(p3,p5)*zab2(p3,p6,p7,p2)*
     &    prop345**(-1)*s267**(-1) + 2*za(p6,p7)*zb(p1,p4)*izb(p1,p5)*
     &    zab2(p3,p1,p4,p1)*zab2(p5,p6,p7,p2)*prop34**(-1)*s267**(-1)*
     &    s134**(-1) + 2*za(p6,p7)*zb(p1,p4)*izb(p1,p5)*zab2(p3,p6,p7,
     &    p2)*zab2(p5,p3,p4,p1)*prop34**(-1)*prop345**(-1)*s267**(-1)
     &     - 2*za(p6,p7)*zb(p1,p4)*izb(p3,p5)*zab2(p5,p6,p7,p2)*
     &    prop345**(-1)*s267**(-1) + 2*zb(p1,p2)*izb(p1,p5)*izb(p5,p6)*
     &    zab2(p3,p5,p6,p1)*zab2(p7,p1,p2,p4)*prop34**(-1)*s127**(-1)
     &     + 2*zb(p1,p4)*izb(p1,p5)*izb(p5,p6)*zab2(p3,p1,p4,p2)*zab2(
     &    p7,p5,p6,p1)*prop34**(-1)*s134**(-1) )
      a(1,2)= + Q(1)*s27**(-1) * (  - 2*za(p3,p5)*za(p6,p7)*zb(p1,p5)*
     &    zb(p4,p5)*iza(p5,p6)*zab2(p6,p2,p7,p2)*prop34**(-1)*
     &    prop345**(-1)*s267**(-1) - za(p3,p6)**2*zb(p1,p2)*zb(p4,p5)*
     &    iza(p5,p6)*zab2(p7,p1,p2,p3)*prop34**(-1)*prop345**(-1)*
     &    s127**(-1) - 2*za(p3,p6)**2*zb(p1,p2)*iza(p3,p5)*iza(p5,p6)*
     &    zab2(p7,p1,p2,p4)*prop345**(-1)*s127**(-1) - za(p3,p6)*za(p4,
     &    p6)*zb(p1,p2)*zb(p4,p5)*iza(p5,p6)*zab2(p7,p1,p2,p4)*
     &    prop34**(-1)*prop345**(-1)*s127**(-1) - 2*za(p3,p6)*za(p6,p3)
     &    *zb(p1,p2)*zb(p3,p4)*iza(p5,p6)*zab2(p7,p1,p2,p5)*
     &    prop34**(-1)*s127**(-1)*s346**(-1) - za(p3,p6)*za(p6,p7)*zb(
     &    p1,p3)*zb(p4,p5)*iza(p5,p6)*zab2(p3,p6,p7,p2)*prop34**(-1)*
     &    prop345**(-1)*s267**(-1) - za(p3,p6)*za(p6,p7)*zb(p1,p4)*zb(
     &    p4,p5)*iza(p5,p6)*zab2(p4,p6,p7,p2)*prop34**(-1)*
     &    prop345**(-1)*s267**(-1) - 2*za(p3,p6)*za(p6,p7)*zb(p1,p4)*
     &    iza(p3,p5)*iza(p5,p6)*zab2(p3,p6,p7,p2)*prop345**(-1)*
     &    s267**(-1) )
      a(1,2) = a(1,2) + Q(1)*s27**(-1) * ( za(p3,p6)*za(p6,p7)*zb(p1,p5
     &    )*zb(p4,p5)*iza(p5,p6)*zab2(p5,p6,p7,p2)*prop34**(-1)*
     &    prop345**(-1)*s267**(-1) + za(p3,p6)*zb(p1,p2)*zb(p4,p5)*
     &    zab2(p7,p1,p2,p5)*prop34**(-1)*prop345**(-1)*s127**(-1) - 2*
     &    za(p3,p6)*zb(p1,p2)*iza(p5,p6)*zab2(p6,p3,p4,p5)*zab2(p7,p1,
     &    p2,p4)*prop34**(-1)*prop345**(-1)*s127**(-1) + 2*za(p3,p6)*
     &    iza(p1,p5)*iza(p5,p6)*zab2(p6,p1,p5,p2)*zab2(p7,p3,p6,p4)*
     &    prop34**(-1)*s346**(-1) - 2*za(p6,p7)*zb(p1,p4)*iza(p5,p6)*
     &    zab2(p3,p6,p7,p2)*zab2(p6,p3,p4,p5)*prop34**(-1)*
     &    prop345**(-1)*s267**(-1) - 2*za(p6,p7)*iza(p1,p5)*iza(p5,p6)*
     &    zab2(p3,p6,p7,p2)*zab2(p6,p1,p5,p4)*prop34**(-1)*s267**(-1) )
      a(1,2) = a(1,2) + Q(2)*s27**(-1) * ( 2*za(p3,p5)*za(p6,p7)*zb(p1,
     &    p5)*zb(p4,p5)*iza(p5,p6)*zab2(p6,p2,p7,p2)*prop34**(-1)*
     &    prop345**(-1)*s267**(-1) + za(p3,p6)**2*zb(p1,p2)*zb(p4,p5)*
     &    iza(p5,p6)*zab2(p7,p1,p2,p3)*prop34**(-1)*prop345**(-1)*
     &    s127**(-1) + 2*za(p3,p6)**2*zb(p1,p2)*iza(p3,p5)*iza(p5,p6)*
     &    zab2(p7,p1,p2,p4)*prop345**(-1)*s127**(-1) + za(p3,p6)*za(p4,
     &    p6)*zb(p1,p2)*zb(p4,p5)*iza(p5,p6)*zab2(p7,p1,p2,p4)*
     &    prop34**(-1)*prop345**(-1)*s127**(-1) + za(p3,p6)*za(p6,p7)*
     &    zb(p1,p3)*zb(p4,p5)*iza(p5,p6)*zab2(p3,p6,p7,p2)*prop34**(-1)
     &    *prop345**(-1)*s267**(-1) + za(p3,p6)*za(p6,p7)*zb(p1,p4)*zb(
     &    p4,p5)*iza(p5,p6)*zab2(p4,p6,p7,p2)*prop34**(-1)*
     &    prop345**(-1)*s267**(-1) + 2*za(p3,p6)*za(p6,p7)*zb(p1,p4)*
     &    iza(p3,p5)*iza(p5,p6)*zab2(p3,p6,p7,p2)*prop345**(-1)*
     &    s267**(-1) - za(p3,p6)*za(p6,p7)*zb(p1,p5)*zb(p4,p5)*iza(p5,
     &    p6)*zab2(p5,p6,p7,p2)*prop34**(-1)*prop345**(-1)*s267**(-1)
     &     - za(p3,p6)*zb(p1,p2)*zb(p4,p5)*zab2(p7,p1,p2,p5)*
     &    prop34**(-1)*prop345**(-1)*s127**(-1) )
      a(1,2) = a(1,2) + Q(2)*s27**(-1) * ( 2*za(p3,p6)*zb(p1,p2)*iza(p5
     &    ,p6)*zab2(p6,p3,p4,p5)*zab2(p7,p1,p2,p4)*prop34**(-1)*
     &    prop345**(-1)*s127**(-1) + 2*za(p6,p7)*zb(p1,p4)*iza(p5,p6)*
     &    zab2(p3,p1,p4,p5)*zab2(p6,p2,p7,p2)*prop34**(-1)*s267**(-1)*
     &    s134**(-1) + 2*za(p6,p7)*zb(p1,p4)*iza(p5,p6)*zab2(p3,p6,p7,
     &    p2)*zab2(p6,p3,p4,p5)*prop34**(-1)*prop345**(-1)*s267**(-1) )
      a(2,1)= + Q(1)*s27**(-1) * (  - za(p2,p6)*za(p3,p5)*zb(p1,p3)*zb(
     &    p1,p4)*izb(p1,p5)*zab2(p3,p2,p6,p7)*prop34**(-1)*
     &    prop345**(-1)*s267**(-1) - za(p2,p6)*za(p3,p5)*zb(p1,p4)**2*
     &    izb(p1,p5)*zab2(p4,p2,p6,p7)*prop34**(-1)*prop345**(-1)*
     &    s267**(-1) + za(p2,p6)*za(p3,p5)*zb(p1,p4)*zab2(p5,p2,p6,p7)*
     &    prop34**(-1)*prop345**(-1)*s267**(-1) - 2*za(p2,p6)*zb(p1,p3)
     &    *zb(p1,p4)*izb(p1,p5)*izb(p3,p5)*zab2(p3,p2,p6,p7)*
     &    prop345**(-1)*s267**(-1) + 2*za(p2,p6)*zb(p1,p4)*izb(p1,p5)*
     &    zab2(p3,p2,p6,p7)*zab2(p5,p3,p4,p1)*prop34**(-1)*
     &    prop345**(-1)*s267**(-1) - 2*za(p2,p6)*zb(p1,p4)*izb(p3,p5)*
     &    zab2(p5,p2,p6,p7)*prop345**(-1)*s267**(-1) + za(p3,p5)*za(p3,
     &    p6)*zb(p1,p4)*zb(p1,p7)*izb(p1,p5)*zab2(p2,p1,p7,p3)*
     &    prop34**(-1)*prop345**(-1)*s127**(-1) + za(p3,p5)*za(p4,p6)*
     &    zb(p1,p4)*zb(p1,p7)*izb(p1,p5)*zab2(p2,p1,p7,p4)*prop34**(-1)
     &    *prop345**(-1)*s127**(-1) - za(p3,p5)*za(p5,p6)*zb(p1,p4)*zb(
     &    p1,p7)*izb(p1,p5)*zab2(p2,p1,p7,p5)*prop34**(-1)*
     &    prop345**(-1)*s127**(-1) )
      a(2,1) = a(2,1) + Q(1)*s27**(-1) * (  - 2*za(p3,p5)*za(p5,p6)*zb(
     &    p1,p7)*zb(p4,p5)*izb(p1,p5)*zab2(p2,p1,p7,p1)*prop34**(-1)*
     &    prop345**(-1)*s127**(-1) + 2*za(p3,p6)*zb(p1,p3)*zb(p1,p7)*
     &    izb(p1,p5)*izb(p3,p5)*zab2(p2,p1,p7,p4)*prop345**(-1)*
     &    s127**(-1) - 2*za(p3,p6)*zb(p1,p7)*izb(p1,p5)*zab2(p2,p1,p7,
     &    p1)*zab2(p5,p3,p6,p4)*prop34**(-1)*s127**(-1)*s346**(-1) - 2*
     &    za(p3,p6)*zb(p1,p7)*izb(p1,p5)*zab2(p2,p1,p7,p4)*zab2(p5,p3,
     &    p4,p1)*prop34**(-1)*prop345**(-1)*s127**(-1) + 2*za(p5,p6)*
     &    zb(p1,p7)*izb(p3,p5)*zab2(p2,p1,p7,p4)*prop345**(-1)*
     &    s127**(-1) )
      a(2,1) = a(2,1) + Q(2)*s27**(-1) * ( za(p2,p6)*za(p3,p5)*zb(p1,p3
     &    )*zb(p1,p4)*izb(p1,p5)*zab2(p3,p2,p6,p7)*prop34**(-1)*
     &    prop345**(-1)*s267**(-1) + za(p2,p6)*za(p3,p5)*zb(p1,p4)**2*
     &    izb(p1,p5)*zab2(p4,p2,p6,p7)*prop34**(-1)*prop345**(-1)*
     &    s267**(-1) - za(p2,p6)*za(p3,p5)*zb(p1,p4)*zab2(p5,p2,p6,p7)*
     &    prop34**(-1)*prop345**(-1)*s267**(-1) + 2*za(p2,p6)*zb(p1,p3)
     &    *zb(p1,p4)*izb(p1,p5)*izb(p3,p5)*zab2(p3,p2,p6,p7)*
     &    prop345**(-1)*s267**(-1) - 2*za(p2,p6)*zb(p1,p4)*izb(p1,p5)*
     &    zab2(p3,p1,p4,p1)*zab2(p5,p2,p6,p7)*prop34**(-1)*s267**(-1)*
     &    s134**(-1) - 2*za(p2,p6)*zb(p1,p4)*izb(p1,p5)*zab2(p3,p2,p6,
     &    p7)*zab2(p5,p3,p4,p1)*prop34**(-1)*prop345**(-1)*s267**(-1)
     &     + 2*za(p2,p6)*zb(p1,p4)*izb(p3,p5)*zab2(p5,p2,p6,p7)*
     &    prop345**(-1)*s267**(-1) - za(p3,p5)*za(p3,p6)*zb(p1,p4)*zb(
     &    p1,p7)*izb(p1,p5)*zab2(p2,p1,p7,p3)*prop34**(-1)*
     &    prop345**(-1)*s127**(-1) - za(p3,p5)*za(p4,p6)*zb(p1,p4)*zb(
     &    p1,p7)*izb(p1,p5)*zab2(p2,p1,p7,p4)*prop34**(-1)*
     &    prop345**(-1)*s127**(-1) )
      a(2,1) = a(2,1) + Q(2)*s27**(-1) * ( za(p3,p5)*za(p5,p6)*zb(p1,p4
     &    )*zb(p1,p7)*izb(p1,p5)*zab2(p2,p1,p7,p5)*prop34**(-1)*
     &    prop345**(-1)*s127**(-1) + 2*za(p3,p5)*za(p5,p6)*zb(p1,p7)*
     &    zb(p4,p5)*izb(p1,p5)*zab2(p2,p1,p7,p1)*prop34**(-1)*
     &    prop345**(-1)*s127**(-1) - 2*za(p3,p6)*zb(p1,p3)*zb(p1,p7)*
     &    izb(p1,p5)*izb(p3,p5)*zab2(p2,p1,p7,p4)*prop345**(-1)*
     &    s127**(-1) + 2*za(p3,p6)*zb(p1,p7)*izb(p1,p5)*zab2(p2,p1,p7,
     &    p4)*zab2(p5,p3,p4,p1)*prop34**(-1)*prop345**(-1)*s127**(-1)
     &     - 2*za(p5,p6)*zb(p1,p7)*izb(p3,p5)*zab2(p2,p1,p7,p4)*
     &    prop345**(-1)*s127**(-1) + 2*zb(p1,p4)*izb(p1,p5)*izb(p5,p6)*
     &    zab2(p2,p5,p6,p1)*zab2(p3,p1,p4,p7)*prop34**(-1)*s134**(-1)
     &     + 2*zb(p1,p7)*izb(p1,p5)*izb(p5,p6)*zab2(p2,p1,p7,p4)*zab2(
     &    p3,p5,p6,p1)*prop34**(-1)*s127**(-1) )
      a(2,2)= + Q(1)*s27**(-1) * ( 2*za(p2,p6)*za(p3,p5)*za(p6,p2)*zb(
     &    p1,p5)*zb(p2,p7)*zb(p4,p5)*iza(p5,p6)*prop34**(-1)*
     &    prop345**(-1)*s267**(-1) + za(p2,p6)*za(p3,p6)*zb(p1,p3)*zb(
     &    p4,p5)*iza(p5,p6)*zab2(p3,p2,p6,p7)*prop34**(-1)*
     &    prop345**(-1)*s267**(-1) + za(p2,p6)*za(p3,p6)*zb(p1,p4)*zb(
     &    p4,p5)*iza(p5,p6)*zab2(p4,p2,p6,p7)*prop34**(-1)*
     &    prop345**(-1)*s267**(-1) + 2*za(p2,p6)*za(p3,p6)*zb(p1,p4)*
     &    iza(p3,p5)*iza(p5,p6)*zab2(p3,p2,p6,p7)*prop345**(-1)*
     &    s267**(-1) - za(p2,p6)*za(p3,p6)*zb(p1,p5)*zb(p4,p5)*iza(p5,
     &    p6)*zab2(p5,p2,p6,p7)*prop34**(-1)*prop345**(-1)*s267**(-1)
     &     + 2*za(p2,p6)*zb(p1,p4)*iza(p5,p6)*zab2(p3,p2,p6,p7)*zab2(p6
     &    ,p3,p4,p5)*prop34**(-1)*prop345**(-1)*s267**(-1) + 2*za(p2,p6
     &    )*iza(p1,p5)*iza(p5,p6)*zab2(p3,p2,p6,p7)*zab2(p6,p1,p5,p4)*
     &    prop34**(-1)*s267**(-1) - za(p3,p6)**2*zb(p1,p7)*zb(p4,p5)*
     &    iza(p5,p6)*zab2(p2,p1,p7,p3)*prop34**(-1)*prop345**(-1)*
     &    s127**(-1) )
      a(2,2) = a(2,2) + Q(1)*s27**(-1) * (  - 2*za(p3,p6)**2*zb(p1,p7)*
     &    iza(p3,p5)*iza(p5,p6)*zab2(p2,p1,p7,p4)*prop345**(-1)*
     &    s127**(-1) - za(p3,p6)*za(p4,p6)*zb(p1,p7)*zb(p4,p5)*iza(p5,
     &    p6)*zab2(p2,p1,p7,p4)*prop34**(-1)*prop345**(-1)*s127**(-1)
     &     - 2*za(p3,p6)*za(p6,p3)*zb(p1,p7)*zb(p3,p4)*iza(p5,p6)*zab2(
     &    p2,p1,p7,p5)*prop34**(-1)*s127**(-1)*s346**(-1) + za(p3,p6)*
     &    zb(p1,p7)*zb(p4,p5)*zab2(p2,p1,p7,p5)*prop34**(-1)*
     &    prop345**(-1)*s127**(-1) - 2*za(p3,p6)*zb(p1,p7)*iza(p5,p6)*
     &    zab2(p2,p1,p7,p4)*zab2(p6,p3,p4,p5)*prop34**(-1)*
     &    prop345**(-1)*s127**(-1) + 2*za(p3,p6)*iza(p1,p5)*iza(p5,p6)*
     &    zab2(p2,p3,p6,p4)*zab2(p6,p1,p5,p7)*prop34**(-1)*s346**(-1) )
      a(2,2) = a(2,2) + Q(2)*s27**(-1) * (  - 2*za(p2,p6)*za(p3,p5)*za(
     &    p6,p2)*zb(p1,p5)*zb(p2,p7)*zb(p4,p5)*iza(p5,p6)*prop34**(-1)*
     &    prop345**(-1)*s267**(-1) - za(p2,p6)*za(p3,p6)*zb(p1,p3)*zb(
     &    p4,p5)*iza(p5,p6)*zab2(p3,p2,p6,p7)*prop34**(-1)*
     &    prop345**(-1)*s267**(-1) - za(p2,p6)*za(p3,p6)*zb(p1,p4)*zb(
     &    p4,p5)*iza(p5,p6)*zab2(p4,p2,p6,p7)*prop34**(-1)*
     &    prop345**(-1)*s267**(-1) - 2*za(p2,p6)*za(p3,p6)*zb(p1,p4)*
     &    iza(p3,p5)*iza(p5,p6)*zab2(p3,p2,p6,p7)*prop345**(-1)*
     &    s267**(-1) + za(p2,p6)*za(p3,p6)*zb(p1,p5)*zb(p4,p5)*iza(p5,
     &    p6)*zab2(p5,p2,p6,p7)*prop34**(-1)*prop345**(-1)*s267**(-1)
     &     - 2*za(p2,p6)*za(p6,p2)*zb(p1,p4)*zb(p2,p7)*iza(p5,p6)*zab2(
     &    p3,p1,p4,p5)*prop34**(-1)*s267**(-1)*s134**(-1) - 2*za(p2,p6)
     &    *zb(p1,p4)*iza(p5,p6)*zab2(p3,p2,p6,p7)*zab2(p6,p3,p4,p5)*
     &    prop34**(-1)*prop345**(-1)*s267**(-1) + za(p3,p6)**2*zb(p1,p7
     &    )*zb(p4,p5)*iza(p5,p6)*zab2(p2,p1,p7,p3)*prop34**(-1)*
     &    prop345**(-1)*s127**(-1) )
      a(2,2) = a(2,2) + Q(2)*s27**(-1) * ( 2*za(p3,p6)**2*zb(p1,p7)*
     &    iza(p3,p5)*iza(p5,p6)*zab2(p2,p1,p7,p4)*prop345**(-1)*
     &    s127**(-1) + za(p3,p6)*za(p4,p6)*zb(p1,p7)*zb(p4,p5)*iza(p5,
     &    p6)*zab2(p2,p1,p7,p4)*prop34**(-1)*prop345**(-1)*s127**(-1)
     &     - za(p3,p6)*zb(p1,p7)*zb(p4,p5)*zab2(p2,p1,p7,p5)*
     &    prop34**(-1)*prop345**(-1)*s127**(-1) + 2*za(p3,p6)*zb(p1,p7)
     &    *iza(p5,p6)*zab2(p2,p1,p7,p4)*zab2(p6,p3,p4,p5)*prop34**(-1)*
     &    prop345**(-1)*s127**(-1) )

c---- Matrix element for W radiation from line16
c---- and gamma radiation line 27
c---- Four diagrams, before and after on both lines
c     d(-p1) +c(-p2) --> mu^-(p3)+vmu~(p4)+gamma(p5)+u(p6)+c(p7)


c           3-----<--4
c               |W
c      6 ----<--|-------------------1
c                    0
c                    0
c                    0
c      7 -----<---------------|-----2
c                             |
c                      gamma  5

c  Overall factor of  i_*rt2*gs^2*e*gw^2*(T^A)_{i7,i1}*(T^A)_{i8,i2}
c     Tr{T^A T^B} = delta^{AB}
c     including QED type propagators

      b(1,1)=2._dp*(+za(p3,p6)*zb(p1,p2)*zbb22(p2,p5,p7,p3,p6,p4)/s346
     & -zb(p1,p4)*zab2(p3,p1,p4,p2)*zab2(p6,p5,p7,p2)/s134)
     & /(zb(p2,p5)*zb(p5,p7)*s257*prop34)

      b(2,1)=2._dp*(-za(p3,p6)*zb(p1,p7)*zbb22(p7,p2,p5,p3,p6,p4)/s346
     & +zb(p1,p4)*zab2(p3,p1,p4,p7)*zab2(p6,p2,p5,p7)/s134)
     & /(zb(p2,p5)*zb(p5,p7)*s257*prop34)

      b(1,2)=-2._dp*(za(p6,p7)*zb(p1,p4)*zaa22(p3,p1,p4,p2,p5,p7)/s134
     & +za(p3,p6)*zab2(p7,p2,p5,p1)*zab2(p7,p3,p6,p4)/s346)
     & /(za(p2,p5)*za(p5,p7)*s257*prop34)

      b(2,2)=2._dp*(-za(p2,p6)*zb(p1,p4)*zaa22(p3,p1,p4,p5,p7,p2)/s134
     & +za(p3,p6)*zab2(p2,p3,p6,p4)*zab2(p2,p5,p7,p1)/s346)
     & /(za(p2,p5)*za(p5,p7)*s257*prop34)
      return
      end
