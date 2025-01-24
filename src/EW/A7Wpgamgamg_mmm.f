!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine A7Wpgamgamg_mmm(p1,p2,p3,p4,p5,p6,p7,b5,b6,b7,za,zb,b)
      implicit none
c---- Matrix element for Wgamgamma radiation from line 12
c  u(-p1) +dbar(-p2) --> ve(p3)+e^+(p4)+gamgam(p5)+g(p6)
c  including radiation from decay electron

c Overall factor of  i_*e^2*gw^2*gs
c     including QED type propagators
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      integer::p1,p2,p3,p4,p5,p6,p7,b5,b6,b7
      real(dp):: s3,s34,s134,s234,s127,s156,s157,
     & s167,s256,s257,s267,s345,s346,s456,Qsum
      complex(dp):: izb,zab2,zab3,b,Mwsq,
     & prop34,prop127,prop1257,prop1267
      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p1,p3)
      zab2(p1,p2,p3,p4)=+za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zab3(p1,p2,p3,p4,p5)=
     & za(p1,p2)*zb(p2,p5)+za(p1,p3)*zb(p3,p5)+za(p1,p4)*zb(p4,p5)
      izb(p1,p2)=cone/zb(p1,p2)

      Qsum=1._dp/3._dp
      s134=s3(p1,p3,p4)
      s127=s3(p1,p2,p7)
      s156=s3(p1,p5,p6)
      s157=s3(p1,p5,p7)
      s167=s3(p1,p6,p7)
      s234=s3(p2,p3,p4)
      s256=s3(p2,p5,p6)
      s257=s3(p2,p5,p7)
      s267=s3(p2,p6,p7)
      s345=s3(p3,p4,p5)
      s346=s3(p3,p4,p6)
      s456=s3(p4,p5,p6)
      s34=s(p3,p4)
      Mwsq=cmplx(wmass**2,-wmass*wwidth,kind=dp)
      prop34=cmplx(s34,kind=dp)-Mwsq
      prop127=cmplx(s127,kind=dp)-Mwsq
      prop1257=cmplx(s346,kind=dp)-Mwsq
      prop1267=cmplx(s345,kind=dp)-Mwsq

      b= + prop34**(-1)*prop127**(-1)*prop1267**(-1)*Mwsq**(-1) * ( za(
     &    p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s345**2 + za(p2,p6)*za(p3,p5
     &    )*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*s(p1,p3)*s345 - 2*za(p2,p6)*za(p3,p5)*zb(p1,b6
     &    )*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7
     &    ,b7)*s(p1,p3)*s(p3,p4) - za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,
     &    b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(
     &    p1,p3)*s(p3,p5) - za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(
     &    p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p1,p3)*
     &    s(p4,p5) + za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*
     &    izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p1,p4)*s345 - 2
     &    *za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)
     &    *izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p1,p4)*s(p3,p4) - za(p2,
     &    p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5
     &    ,b5)*izb(p6,b6)*izb(p7,b7)*s(p1,p4)*s(p3,p5) )
      b = b + prop34**(-1)*prop127**(-1)*prop1267**(-1)*Mwsq**(-1) * (
     &     - za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p1,p4)*s(p4,p5) + za(
     &    p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p1,p5)*s345 - 2*za(p2,p6)*
     &    za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)
     &    *izb(p6,b6)*izb(p7,b7)*s(p1,p5)*s(p3,p4) - za(p2,p6)*za(p3,p5
     &    )*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*s(p1,p5)*s(p3,p5) - za(p2,p6)*za(p3,p5)*zb(p1,
     &    b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*s(p1,p5)*s(p4,p5) + za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(
     &    p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    s(p2,p3)*s345 - 2*za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(
     &    p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p2,p3)*
     &    s(p3,p4) - za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*
     &    izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p2,p3)*s(p3,p5)
     &     )
      b = b + prop34**(-1)*prop127**(-1)*prop1267**(-1)*Mwsq**(-1) * (
     &     - za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p2,p3)*s(p4,p5) + za(
     &    p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p2,p4)*s345 - 2*za(p2,p6)*
     &    za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)
     &    *izb(p6,b6)*izb(p7,b7)*s(p2,p4)*s(p3,p4) - za(p2,p6)*za(p3,p5
     &    )*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*s(p2,p4)*s(p3,p5) - za(p2,p6)*za(p3,p5)*zb(p1,
     &    b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*s(p2,p4)*s(p4,p5) + za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(
     &    p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    s(p2,p5)*s345 - 2*za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(
     &    p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p2,p5)*
     &    s(p3,p4) - za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*
     &    izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p2,p5)*s(p3,p5)
     &     )
      b = b + prop34**(-1)*prop127**(-1)*prop1267**(-1)*Mwsq**(-1) * (
     &     - za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p2,p5)*s(p4,p5) - 2*
     &    za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p3,p4)*s345 - 2*za(p2,p6)*
     &    za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)
     &    *izb(p6,b6)*izb(p7,b7)*s(p3,p4)*s(p3,p7) - 2*za(p2,p6)*za(p3,
     &    p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*s(p3,p4)*s(p4,p7) - 2*za(p2,p6)*za(p3,p5)*
     &    zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6
     &    )*izb(p7,b7)*s(p3,p4)*s(p5,p7) - za(p2,p6)*za(p3,p5)*zb(p1,b6
     &    )*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7
     &    ,b7)*s(p3,p5)*s345 - za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*
     &    zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p3,p5
     &    )*s(p3,p7) - za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5
     &    )*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p3,p5)*s(p4,
     &    p7) )
      b = b + prop34**(-1)*prop127**(-1)*prop1267**(-1)*Mwsq**(-1) * (
     &     - za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p3,p5)*s(p5,p7) + za(
     &    p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p3,p7)*s345 - za(p2,p6)*
     &    za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)
     &    *izb(p6,b6)*izb(p7,b7)*s(p3,p7)*s(p4,p5) - za(p2,p6)*za(p3,p5
     &    )*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*s(p4,p5)*s345 - za(p2,p6)*za(p3,p5)*zb(p1,b6)*
     &    zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*s(p4,p5)*s(p4,p7) - za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,
     &    b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(
     &    p4,p5)*s(p5,p7) + za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(
     &    p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p4,p7)*
     &    s345 + za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(
     &    p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p5,p7)*s345 )
      b = b + prop34**(-1)*prop127**(-1)*prop1267**(-1)*Mwsq**(-1) * (
     &    za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6
     &    )*izb(p7,b7)*zab2(p6,p2,p7,b7)*s345**2 + za(p3,p5)*zb(p1,b6)*
     &    zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6
     &    ,p2,p7,b7)*s(p1,p3)*s345 - 2*za(p3,p5)*zb(p1,b6)*zb(p4,b5)*
     &    izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)
     &    *s(p1,p3)*s(p3,p4) - za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)
     &    *izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p1,p3)*
     &    s(p3,p5) - za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5
     &    )*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p1,p3)*s(p4,p5)
     &     + za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6
     &    ,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p1,p4)*s345 - 2*za(p3,p5)
     &    *zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*zab2(p6,p2,p7,b7)*s(p1,p4)*s(p3,p4) - za(p3,p5)*zb(p1,b6)
     &    *zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(
     &    p6,p2,p7,b7)*s(p1,p4)*s(p3,p5) )
      b = b + prop34**(-1)*prop127**(-1)*prop1267**(-1)*Mwsq**(-1) * (
     &     - za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6
     &    ,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p1,p4)*s(p4,p5) + za(p3,
     &    p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*zab2(p6,p2,p7,b7)*s(p1,p5)*s345 - 2*za(p3,p5)*zb(p1,b6
     &    )*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(
     &    p6,p2,p7,b7)*s(p1,p5)*s(p3,p4) - za(p3,p5)*zb(p1,b6)*zb(p4,b5
     &    )*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,
     &    b7)*s(p1,p5)*s(p3,p5) - za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p1,
     &    p5)*s(p4,p5) + za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(
     &    p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p2,p3)*s345
     &     - 2*za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p2,p3)*s(p3,p4) - za(p3
     &    ,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p2,p3)*s(p3,p5) )
      b = b + prop34**(-1)*prop127**(-1)*prop1267**(-1)*Mwsq**(-1) * (
     &     - za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6
     &    ,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p2,p3)*s(p4,p5) + za(p3,
     &    p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*zab2(p6,p2,p7,b7)*s(p2,p4)*s345 - 2*za(p3,p5)*zb(p1,b6
     &    )*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(
     &    p6,p2,p7,b7)*s(p2,p4)*s(p3,p4) - za(p3,p5)*zb(p1,b6)*zb(p4,b5
     &    )*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,
     &    b7)*s(p2,p4)*s(p3,p5) - za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p2,
     &    p4)*s(p4,p5) + za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(
     &    p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p2,p5)*s345
     &     - 2*za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p2,p5)*s(p3,p4) - za(p3
     &    ,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p2,p5)*s(p3,p5) )
      b = b + prop34**(-1)*prop127**(-1)*prop1267**(-1)*Mwsq**(-1) * (
     &     - za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6
     &    ,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p2,p5)*s(p4,p5) - 2*za(p3
     &    ,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p3,p4)*s345 - 2*za(p3,p5)*zb(
     &    p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p6,p2,p7,b7)*s(p3,p4)*s(p3,p7) - 2*za(p3,p5)*zb(p1,b6)*
     &    zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6
     &    ,p2,p7,b7)*s(p3,p4)*s(p4,p7) - 2*za(p3,p5)*zb(p1,b6)*zb(p4,b5
     &    )*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,
     &    b7)*s(p3,p4)*s(p5,p7) - za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p3,
     &    p5)*s345 - za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5
     &    )*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p3,p5)*s(p3,p7)
     &     - za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6
     &    ,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p3,p5)*s(p4,p7) )
      b = b + prop34**(-1)*prop127**(-1)*prop1267**(-1)*Mwsq**(-1) * (
     &     - za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6
     &    ,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p3,p5)*s(p5,p7) + za(p3,
     &    p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*zab2(p6,p2,p7,b7)*s(p3,p7)*s345 - za(p3,p5)*zb(p1,b6)*
     &    zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6
     &    ,p2,p7,b7)*s(p3,p7)*s(p4,p5) - za(p3,p5)*zb(p1,b6)*zb(p4,b5)*
     &    izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)
     &    *s(p4,p5)*s345 - za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p4,p5)*
     &    s(p4,p7) - za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5
     &    )*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p4,p5)*s(p5,p7)
     &     + za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6
     &    ,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p4,p7)*s345 + za(p3,p5)*
     &    zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*zab2(p6,p2,p7,b7)*s(p5,p7)*s345 )
      b = b + prop34**(-1)*prop127**(-1)*prop1267**(-1) * (  - 2*za(p2,
     &    p3)*za(p3,p5)*zb(p1,p3)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5
     &    ,b5)*izb(p6,b6)*izb(p7,b7)*zab3(p6,p1,p2,p7,b6) + 2*za(p2,p3)
     &    *za(p3,p5)*zb(p1,p3)*zb(p2,b7)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5
     &    )*izb(p6,b6)*izb(p7,b7)*zab3(p6,p1,p2,p7,b6) - 2*za(p2,p3)*
     &    zb(p1,p4)*zb(p1,b7)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*zab2(p5,p3,p4,b5)*zab3(p6,p1,p2,p7,b6) - 2*za(p2,p4)*za(
     &    p3,p5)*zb(p1,p4)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab3(p6,p1,p2,p7,b6) + 2*za(p2,p4)*za(
     &    p3,p5)*zb(p1,p4)*zb(p2,b7)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab3(p6,p1,p2,p7,b6) - 2*za(p2,p5)*za(
     &    p3,p5)*zb(p1,b5)*zb(p1,b7)*zb(p4,p5)*izb(p1,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab3(p6,p1,p2,p7,b6) + 2*za(p2,p6)*za(
     &    p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,p5)*izb(p1,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab3(p5,p1,p2,p7,b5) - za(p2,p6)*za(p3,
     &    p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*s345 )
      b = b + prop34**(-1)*prop127**(-1)*prop1267**(-1) * ( za(p2,p6)*
     &    za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)
     &    *izb(p6,b6)*izb(p7,b7)*s(p1,p3) + za(p2,p6)*za(p3,p5)*zb(p1,
     &    b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*s(p1,p4) - za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(
     &    p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p1,p5)
     &     + za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p2,p3) + za(p2,p6)*za(
     &    p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*s(p2,p4) - za(p2,p6)*za(p3,p5)*zb(p1,b6
     &    )*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7
     &    ,b7)*s(p2,p5) + 2*za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(
     &    p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p3,p4)
     &     + za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p3,p5) + za(p2,p6)*za(
     &    p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*s(p3,p7) )
      b = b + prop34**(-1)*prop127**(-1)*prop1267**(-1) * ( za(p2,p6)*
     &    za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)
     &    *izb(p6,b6)*izb(p7,b7)*s(p4,p5) + za(p2,p6)*za(p3,p5)*zb(p1,
     &    b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*s(p4,p7) - za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(
     &    p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p5,p7)
     &     - 2*za(p2,p6)*zb(p1,b6)*zb(p1,b7)*izb(p1,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*zab2(p3,p5,p6,p4)*zab2(p5,p3,p4,b5) + 2*za(
     &    p2,p7)*za(p3,p5)*za(p5,p6)*zb(p1,b7)*zb(p4,p5)*zb(b5,b6)*izb(
     &    p5,b5)*izb(p6,b6)*izb(p7,b7) - 2*za(p2,p7)*za(p3,p5)*zb(p1,b7
     &    )*zb(p4,b5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p3,p4,b6
     &    ) + 2*za(p2,p7)*za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*zab2(p5,p3,p4,b5) - 2*za(p3,p5)*za(p3,p7)*
     &    zb(p1,p3)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*zab3(p6,
     &    p1,p2,p7,b6) - 2*za(p3,p5)*za(p4,p7)*zb(p1,p4)*zb(p4,b5)*izb(
     &    p2,p7)*izb(p5,b5)*izb(p6,b6)*zab3(p6,p1,p2,p7,b6) )
      b = b + prop34**(-1)*prop127**(-1)*prop1267**(-1) * (  - 2*za(p3,
     &    p5)*za(p5,p6)*zb(p1,p6)*zb(p4,p5)*zb(b5,b6)*izb(p2,p7)*izb(p5
     &    ,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7) + 2*za(p3,p5)*
     &    za(p5,p6)*zb(p1,b7)*zb(p4,p5)*zb(b5,b6)*izb(p1,p7)*izb(p5,b5)
     &    *izb(p6,b6)*izb(p7,b7)*zab2(p2,p6,p7,p1) + 2*za(p3,p5)*zb(p1,
     &    p6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p6,p2,p7,b7)*zab2(p6,p3,p4,b6) - 2*za(p3,p5)*zb(p1,b5)*
     &    zb(p4,p5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5
     &    ,p2,p7,b7)*zab3(p6,p1,p2,p7,b6) + 2*za(p3,p5)*zb(p1,b6)*zb(p4
     &    ,p5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,
     &    p7,b7)*zab3(p5,p1,p2,p7,b5) - za(p3,p5)*zb(p1,b6)*zb(p4,b5)*
     &    izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)
     &    *s345 + za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p1,p3) + za(p3,p5)*
     &    zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*zab2(p6,p2,p7,b7)*s(p1,p4) )
      b = b + prop34**(-1)*prop127**(-1)*prop1267**(-1) * (  - za(p3,p5
     &    )*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7
     &    ,b7)*zab2(p6,p2,p7,b7)*s(p1,p5) + za(p3,p5)*zb(p1,b6)*zb(p4,
     &    b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7
     &    ,b7)*s(p2,p3) + za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(
     &    p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p2,p4) - za(
     &    p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p2,p5) + 2*za(p3,p5)*zb(p1,b6)
     &    *zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(
     &    p6,p2,p7,b7)*s(p3,p4) + za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p3,
     &    p5) + za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p3,p7) + za(p3,p5)*
     &    zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*zab2(p6,p2,p7,b7)*s(p4,p5) + za(p3,p5)*zb(p1,b6)*zb(p4,b5
     &    )*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,
     &    b7)*s(p4,p7) )
      b = b + prop34**(-1)*prop127**(-1)*prop1267**(-1) * (  - za(p3,p5
     &    )*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7
     &    ,b7)*zab2(p6,p2,p7,b7)*s(p5,p7) - 2*za(p3,p5)*zb(p1,b7)*zb(p4
     &    ,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p6,
     &    p7,p1)*zab2(p6,p3,p4,b6) - 2*za(p3,p6)*zb(p1,p6)*zb(p4,b6)*
     &    izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p3,p4,b5)
     &    *zab2(p6,p2,p7,b7) + 2*za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p1,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p6,p7,p1)*zab2(
     &    p5,p3,p4,b5) - 2*zb(p1,p4)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p3,p2,p7,b7)*zab2(p5,p3,p4,b5)*zab3(p6,p1,p2,
     &    p7,b6) - 2*zb(p1,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*zab2(p3,p5,p6,p4)*zab2(p5,p3,p4,b5)*zab2(p6,p2,p7,b7) )
      b = b + prop34**(-1)*prop127**(-1)*prop1257**(-1)*Mwsq**(-1) * (
     &    za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s346**2 + za(p2,p5)*za(p3,p6
     &    )*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*s(p1,p3)*s346 - 2*za(p2,p5)*za(p3,p6)*zb(p1,b5
     &    )*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7
     &    ,b7)*s(p1,p3)*s(p3,p4) - za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,
     &    b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(
     &    p1,p3)*s(p3,p6) - za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(
     &    p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p1,p3)*
     &    s(p4,p6) + za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*
     &    izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p1,p4)*s346 - 2
     &    *za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)
     &    *izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p1,p4)*s(p3,p4) - za(p2,
     &    p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5
     &    ,b5)*izb(p6,b6)*izb(p7,b7)*s(p1,p4)*s(p3,p6) )
      b = b + prop34**(-1)*prop127**(-1)*prop1257**(-1)*Mwsq**(-1) * (
     &     - za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p1,p4)*s(p4,p6) + za(
     &    p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p1,p6)*s346 - 2*za(p2,p5)*
     &    za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)
     &    *izb(p6,b6)*izb(p7,b7)*s(p1,p6)*s(p3,p4) - za(p2,p5)*za(p3,p6
     &    )*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*s(p1,p6)*s(p3,p6) - za(p2,p5)*za(p3,p6)*zb(p1,
     &    b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*s(p1,p6)*s(p4,p6) + za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(
     &    p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    s(p2,p3)*s346 - 2*za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(
     &    p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p2,p3)*
     &    s(p3,p4) - za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*
     &    izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p2,p3)*s(p3,p6)
     &     )
      b = b + prop34**(-1)*prop127**(-1)*prop1257**(-1)*Mwsq**(-1) * (
     &     - za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p2,p3)*s(p4,p6) + za(
     &    p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p2,p4)*s346 - 2*za(p2,p5)*
     &    za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)
     &    *izb(p6,b6)*izb(p7,b7)*s(p2,p4)*s(p3,p4) - za(p2,p5)*za(p3,p6
     &    )*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*s(p2,p4)*s(p3,p6) - za(p2,p5)*za(p3,p6)*zb(p1,
     &    b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*s(p2,p4)*s(p4,p6) + za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(
     &    p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    s(p2,p6)*s346 - 2*za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(
     &    p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p2,p6)*
     &    s(p3,p4) - za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*
     &    izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p2,p6)*s(p3,p6)
     &     )
      b = b + prop34**(-1)*prop127**(-1)*prop1257**(-1)*Mwsq**(-1) * (
     &     - za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p2,p6)*s(p4,p6) - 2*
     &    za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p3,p4)*s346 - 2*za(p2,p5)*
     &    za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)
     &    *izb(p6,b6)*izb(p7,b7)*s(p3,p4)*s(p3,p7) - 2*za(p2,p5)*za(p3,
     &    p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*s(p3,p4)*s(p4,p7) - 2*za(p2,p5)*za(p3,p6)*
     &    zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6
     &    )*izb(p7,b7)*s(p3,p4)*s(p6,p7) - za(p2,p5)*za(p3,p6)*zb(p1,b5
     &    )*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7
     &    ,b7)*s(p3,p6)*s346 - za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*
     &    zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p3,p6
     &    )*s(p3,p7) - za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6
     &    )*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p3,p6)*s(p4,
     &    p7) )
      b = b + prop34**(-1)*prop127**(-1)*prop1257**(-1)*Mwsq**(-1) * (
     &     - za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p3,p6)*s(p6,p7) + za(
     &    p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p3,p7)*s346 - za(p2,p5)*
     &    za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)
     &    *izb(p6,b6)*izb(p7,b7)*s(p3,p7)*s(p4,p6) - za(p2,p5)*za(p3,p6
     &    )*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*s(p4,p6)*s346 - za(p2,p5)*za(p3,p6)*zb(p1,b5)*
     &    zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*s(p4,p6)*s(p4,p7) - za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,
     &    b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(
     &    p4,p6)*s(p6,p7) + za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(
     &    p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p4,p7)*
     &    s346 + za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(
     &    p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p6,p7)*s346 )
      b = b + prop34**(-1)*prop127**(-1)*prop1257**(-1)*Mwsq**(-1) * (
     &     - 1.D0/2.D0*za(p2,p7)*za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p5,
     &    b5)*izb(p6,b6)*izb(p7,b7)*zab3(p5,p1,p2,p7,b5)*s346 + za(p2,
     &    p7)*za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*zab3(p5,p1,p2,p7,b5)*s(p3,p4) + 1.D0/2.D0*za(p2,p7)*
     &    za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7
     &    )*zab3(p5,p1,p2,p7,b5)*s(p3,p6) + 1.D0/2.D0*za(p2,p7)*za(p3,
     &    p6)*zb(p1,b7)*zb(p4,b6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab3(p5,p1,p2,p7,b5)*s(p4,p6) - 1.D0/2.D0*za(p2,p7)*za(p3,p6)
     &    *zb(p1,b7)*zb(p4,b6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab3(p5
     &    ,p3,p4,p6,b5)*s346 + za(p2,p7)*za(p3,p6)*zb(p1,b7)*zb(p4,b6)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab3(p5,p3,p4,p6,b5)*s(p3,p4
     &    ) + 1.D0/2.D0*za(p2,p7)*za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p5,
     &    b5)*izb(p6,b6)*izb(p7,b7)*zab3(p5,p3,p4,p6,b5)*s(p3,p6) + 1.D0
     &    /2.D0*za(p2,p7)*za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*zab3(p5,p3,p4,p6,b5)*s(p4,p6) )
      b = b + prop34**(-1)*prop127**(-1)*prop1257**(-1)*Mwsq**(-1) * (
     &    1.D0/2.D0*za(p3,p6)*zb(p1,p5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)
     &    *izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*zab3(p5,p1,p2,p7,b5)
     &    *s346 - za(p3,p6)*zb(p1,p5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*zab3(p5,p1,p2,p7,b5)*
     &    s(p3,p4) - 1.D0/2.D0*za(p3,p6)*zb(p1,p5)*zb(p4,b6)*izb(p2,p7)
     &    *izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*zab3(p5,
     &    p1,p2,p7,b5)*s(p3,p6) - 1.D0/2.D0*za(p3,p6)*zb(p1,p5)*zb(p4,
     &    b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7
     &    ,b7)*zab3(p5,p1,p2,p7,b5)*s(p4,p6) + 1.D0/2.D0*za(p3,p6)*zb(
     &    p1,p5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p5,p2,p7,b7)*zab3(p5,p3,p4,p6,b5)*s346 - za(p3,p6)*zb(p1
     &    ,p5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p5,p2,p7,b7)*zab3(p5,p3,p4,p6,b5)*s(p3,p4) - 1.D0/2.D0*
     &    za(p3,p6)*zb(p1,p5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6
     &    )*izb(p7,b7)*zab2(p5,p2,p7,b7)*zab3(p5,p3,p4,p6,b5)*s(p3,p6)
     &     )
      b = b + prop34**(-1)*prop127**(-1)*prop1257**(-1)*Mwsq**(-1) * (
     &     - 1.D0/2.D0*za(p3,p6)*zb(p1,p5)*zb(p4,b6)*izb(p2,p7)*izb(p5,
     &    b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*zab3(p5,p3,p4,p6,
     &    b5)*s(p4,p6) + za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(
     &    p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s346**2 + za(
     &    p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p1,p3)*s346 - 2*za(p3,p6)*zb(
     &    p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p5,p2,p7,b7)*s(p1,p3)*s(p3,p4) - za(p3,p6)*zb(p1,b5)*zb(
     &    p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2
     &    ,p7,b7)*s(p1,p3)*s(p3,p6) - za(p3,p6)*zb(p1,b5)*zb(p4,b6)*
     &    izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)
     &    *s(p1,p3)*s(p4,p6) + za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)
     &    *izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p1,p4)*
     &    s346 - 2*za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p1,p4)*s(p3,p4) )
      b = b + prop34**(-1)*prop127**(-1)*prop1257**(-1)*Mwsq**(-1) * (
     &     - za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6
     &    ,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p1,p4)*s(p3,p6) - za(p3,
     &    p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*zab2(p5,p2,p7,b7)*s(p1,p4)*s(p4,p6) + za(p3,p6)*zb(p1,
     &    b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p5,p2,p7,b7)*s(p1,p6)*s346 - 2*za(p3,p6)*zb(p1,b5)*zb(p4
     &    ,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,
     &    p7,b7)*s(p1,p6)*s(p3,p4) - za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(
     &    p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(
     &    p1,p6)*s(p3,p6) - za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p1,p6)*
     &    s(p4,p6) + za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5
     &    )*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p2,p3)*s346 - 2*
     &    za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6
     &    )*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p2,p3)*s(p3,p4) )
      b = b + prop34**(-1)*prop127**(-1)*prop1257**(-1)*Mwsq**(-1) * (
     &     - za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6
     &    ,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p2,p3)*s(p3,p6) - za(p3,
     &    p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*zab2(p5,p2,p7,b7)*s(p2,p3)*s(p4,p6) + za(p3,p6)*zb(p1,
     &    b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p5,p2,p7,b7)*s(p2,p4)*s346 - 2*za(p3,p6)*zb(p1,b5)*zb(p4
     &    ,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,
     &    p7,b7)*s(p2,p4)*s(p3,p4) - za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(
     &    p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(
     &    p2,p4)*s(p3,p6) - za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p2,p4)*
     &    s(p4,p6) + za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5
     &    )*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p2,p6)*s346 - 2*
     &    za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6
     &    )*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p2,p6)*s(p3,p4) )
      b = b + prop34**(-1)*prop127**(-1)*prop1257**(-1)*Mwsq**(-1) * (
     &     - za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6
     &    ,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p2,p6)*s(p3,p6) - za(p3,
     &    p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*zab2(p5,p2,p7,b7)*s(p2,p6)*s(p4,p6) - 2*za(p3,p6)*zb(
     &    p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p5,p2,p7,b7)*s(p3,p4)*s346 - 2*za(p3,p6)*zb(p1,b5)*zb(p4
     &    ,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,
     &    p7,b7)*s(p3,p4)*s(p3,p7) - 2*za(p3,p6)*zb(p1,b5)*zb(p4,b6)*
     &    izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)
     &    *s(p3,p4)*s(p4,p7) - 2*za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p3,
     &    p4)*s(p6,p7) - za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(
     &    p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p3,p6)*s346
     &     - za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6
     &    ,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p3,p6)*s(p3,p7) )
      b = b + prop34**(-1)*prop127**(-1)*prop1257**(-1)*Mwsq**(-1) * (
     &     - za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6
     &    ,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p3,p6)*s(p4,p7) - za(p3,
     &    p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*zab2(p5,p2,p7,b7)*s(p3,p6)*s(p6,p7) + za(p3,p6)*zb(p1,
     &    b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p5,p2,p7,b7)*s(p3,p7)*s346 - za(p3,p6)*zb(p1,b5)*zb(p4,
     &    b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7
     &    ,b7)*s(p3,p7)*s(p4,p6) - za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2
     &    ,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p4,
     &    p6)*s346 - za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5
     &    )*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p4,p6)*s(p4,p7)
     &     - za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6
     &    ,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p4,p6)*s(p6,p7) + za(p3,
     &    p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*zab2(p5,p2,p7,b7)*s(p4,p7)*s346 )
      b = b + prop34**(-1)*prop127**(-1)*prop1257**(-1)*Mwsq**(-1) * (
     &    za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6
     &    )*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p6,p7)*s346 - 1.D0/2.D0*za(
     &    p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p2,p5,p7,p1)*zab3(p5,p1,p2,p7,b5)*s346 + za(
     &    p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p2,p5,p7,p1)*zab3(p5,p1,p2,p7,b5)*s(p3,p4) +
     &    1.D0/2.D0*za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)
     &    *izb(p6,b6)*izb(p7,b7)*zab2(p2,p5,p7,p1)*zab3(p5,p1,p2,p7,b5)
     &    *s(p3,p6) + 1.D0/2.D0*za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7
     &    )*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p5,p7,p1)*zab3(p5,
     &    p1,p2,p7,b5)*s(p4,p6) - 1.D0/2.D0*za(p3,p6)*zb(p1,b7)*zb(p4,
     &    b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p5,p7
     &    ,p1)*zab3(p5,p3,p4,p6,b5)*s346 + za(p3,p6)*zb(p1,b7)*zb(p4,b6
     &    )*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p5,p7,
     &    p1)*zab3(p5,p3,p4,p6,b5)*s(p3,p4) )
      b = b + prop34**(-1)*prop127**(-1)*prop1257**(-1)*Mwsq**(-1) * (
     &    1.D0/2.D0*za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)
     &    *izb(p6,b6)*izb(p7,b7)*zab2(p2,p5,p7,p1)*zab3(p5,p3,p4,p6,b5)
     &    *s(p3,p6) + 1.D0/2.D0*za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7
     &    )*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p5,p7,p1)*zab3(p5,
     &    p3,p4,p6,b5)*s(p4,p6) )
      b = b + prop34**(-1)*prop127**(-1)*prop1257**(-1) * (  - za(p2,p3
     &    )*za(p3,p6)*zb(p1,p3)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,
     &    b5)*izb(p6,b6)*izb(p7,b7)*zab3(p5,p1,p2,p7,b5) + za(p2,p3)*
     &    za(p3,p6)*zb(p1,p3)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)
     &    *izb(p6,b6)*izb(p7,b7)*zab3(p5,p3,p4,p6,b5) + za(p2,p3)*za(p3
     &    ,p6)*zb(p1,p3)*zb(p2,b7)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*zab3(p5,p1,p2,p7,b5) - za(p2,p3)*za(p3,p6)*
     &    zb(p1,p3)*zb(p2,b7)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6
     &    )*izb(p7,b7)*zab3(p5,p3,p4,p6,b5) - za(p2,p3)*zb(p1,p4)*zb(p1
     &    ,b7)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p3,
     &    p4,b6)*zab3(p5,p1,p2,p7,b5) + za(p2,p3)*zb(p1,p4)*zb(p1,b7)*
     &    izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p3,p4,b6)
     &    *zab3(p5,p3,p4,p6,b5) - za(p2,p4)*za(p3,p6)*zb(p1,p4)*zb(p1,
     &    b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab3(p5,p1,p2,p7,b5) + za(p2,p4)*za(p3,p6)*zb(p1,p4)*zb(p1,b7
     &    )*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab3(
     &    p5,p3,p4,p6,b5) )
      b = b + prop34**(-1)*prop127**(-1)*prop1257**(-1) * ( za(p2,p4)*
     &    za(p3,p6)*zb(p1,p4)*zb(p2,b7)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)
     &    *izb(p6,b6)*izb(p7,b7)*zab3(p5,p1,p2,p7,b5) - za(p2,p4)*za(p3
     &    ,p6)*zb(p1,p4)*zb(p2,b7)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*zab3(p5,p3,p4,p6,b5) + 2*za(p2,p5)*za(p3,p6
     &    )*zb(p1,b5)*zb(p1,b7)*zb(p4,p6)*izb(p1,p7)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*zab3(p6,p1,p2,p7,b6) - za(p2,p5)*za(p3,p6)*zb(
     &    p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*s346 + za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(
     &    p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p1,p3)
     &     + za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p1,p4) - za(p2,p5)*za(
     &    p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*s(p1,p6) + za(p2,p5)*za(p3,p6)*zb(p1,b5
     &    )*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7
     &    ,b7)*s(p2,p3) )
      b = b + prop34**(-1)*prop127**(-1)*prop1257**(-1) * ( za(p2,p5)*
     &    za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)
     &    *izb(p6,b6)*izb(p7,b7)*s(p2,p4) - za(p2,p5)*za(p3,p6)*zb(p1,
     &    b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*s(p2,p6) + 2*za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*
     &    zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p3,p4
     &    ) + za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p3,p6) + za(p2,p5)*za(
     &    p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*s(p3,p7) + za(p2,p5)*za(p3,p6)*zb(p1,b5
     &    )*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7
     &    ,b7)*s(p4,p6) + za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4
     &    ,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p4,p7) -
     &    za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p6,p7) - 2*za(p2,p5)*zb(p1
     &    ,b5)*zb(p1,b7)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p3,p5,p6,p4)*zab2(p6,p3,p4,b6) )
      b = b + prop34**(-1)*prop127**(-1)*prop1257**(-1) * (  - za(p2,p6
     &    )*za(p3,p6)*zb(p1,b6)*zb(p1,b7)*zb(p4,p6)*izb(p1,p7)*izb(p5,
     &    b5)*izb(p6,b6)*izb(p7,b7)*zab3(p5,p1,p2,p7,b5) + za(p2,p6)*
     &    za(p3,p6)*zb(p1,b6)*zb(p1,b7)*zb(p4,p6)*izb(p1,p7)*izb(p5,b5)
     &    *izb(p6,b6)*izb(p7,b7)*zab3(p5,p3,p4,p6,b5) + 2*za(p2,p7)*za(
     &    p3,p5)*zb(p1,b7)*zb(p4,b5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p6,p3,p4,b6) + 2*za(p2,p7)*za(p3,p6)*za(p5,p6)*zb(p1,b7)
     &    *zb(p4,p6)*zb(b5,b6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7) - 2*za(
     &    p2,p7)*za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p5,p3,p4,b5) + 1.D0/2.D0*za(p2,p7)*za(p3,p6)*
     &    zb(p1,b7)*zb(p4,b6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab3(p5,
     &    p1,p2,p7,b5) + 1.D0/2.D0*za(p2,p7)*za(p3,p6)*zb(p1,b7)*zb(p4,
     &    b6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab3(p5,p3,p4,p6,b5) - 2
     &    *za(p3,p5)*zb(p1,p5)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*zab2(p6,p3,p4,b6) + 2*za(p3,
     &    p5)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*zab2(p2,p5,p7,p1)*zab2(p6,p3,p4,b6) )
      b = b + prop34**(-1)*prop127**(-1)*prop1257**(-1) * (  - za(p3,p6
     &    )*za(p3,p7)*zb(p1,p3)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,
     &    b6)*zab3(p5,p1,p2,p7,b5) + za(p3,p6)*za(p3,p7)*zb(p1,p3)*zb(
     &    p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*zab3(p5,p3,p4,p6,b5)
     &     - za(p3,p6)*za(p4,p7)*zb(p1,p4)*zb(p4,b6)*izb(p2,p7)*izb(p5,
     &    b5)*izb(p6,b6)*zab3(p5,p1,p2,p7,b5) + za(p3,p6)*za(p4,p7)*zb(
     &    p1,p4)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*zab3(p5,p3,
     &    p4,p6,b5) - 2*za(p3,p6)*za(p5,p6)*zb(p1,p5)*zb(p4,p6)*zb(b5,
     &    b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7
     &    ,b7) + 2*za(p3,p6)*za(p5,p6)*zb(p1,b7)*zb(p4,p6)*zb(b5,b6)*
     &    izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p5,p7,p1)
     &     + 2*za(p3,p6)*zb(p1,p5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*zab2(p5,p3,p4,b5) - 1.D0/
     &    2.D0*za(p3,p6)*zb(p1,p5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*zab3(p5,p1,p2,p7,b5) - 1.D
     &    0/2.D0*za(p3,p6)*zb(p1,p5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*zab3(p5,p3,p4,p6,b5)
     &     )
      b = b + prop34**(-1)*prop127**(-1)*prop1257**(-1) * ( 2*za(p3,p6)
     &    *zb(p1,b5)*zb(p4,p6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*zab2(p5,p2,p7,b7)*zab3(p6,p1,p2,p7,b6) - za(p3,p6)*zb(p1,
     &    b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p5,p2,p7,b7)*s346 + za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2
     &    ,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p1,
     &    p3) + za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p1,p4) - za(p3,p6)*
     &    zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*zab2(p5,p2,p7,b7)*s(p1,p6) + za(p3,p6)*zb(p1,b5)*zb(p4,b6
     &    )*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,
     &    b7)*s(p2,p3) + za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(
     &    p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p2,p4) - za(
     &    p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p2,p6) + 2*za(p3,p6)*zb(p1,b5)
     &    *zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(
     &    p5,p2,p7,b7)*s(p3,p4) )
      b = b + prop34**(-1)*prop127**(-1)*prop1257**(-1) * ( za(p3,p6)*
     &    zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*zab2(p5,p2,p7,b7)*s(p3,p6) + za(p3,p6)*zb(p1,b5)*zb(p4,b6
     &    )*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,
     &    b7)*s(p3,p7) + za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(
     &    p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p4,p6) + za(
     &    p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p4,p7) - za(p3,p6)*zb(p1,b5)*
     &    zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5
     &    ,p2,p7,b7)*s(p6,p7) - za(p3,p6)*zb(p1,b6)*zb(p4,p6)*izb(p2,p7
     &    )*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*zab3(p5,
     &    p1,p2,p7,b5) + za(p3,p6)*zb(p1,b6)*zb(p4,p6)*izb(p2,p7)*izb(
     &    p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*zab3(p5,p3,p4,
     &    p6,b5) - 2*za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5
     &    )*izb(p6,b6)*izb(p7,b7)*zab2(p2,p5,p7,p1)*zab2(p5,p3,p4,b5)
     &     + 1.D0/2.D0*za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,
     &    b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p5,p7,p1)*zab3(p5,p1,p2,p7,
     &    b5) )
      b = b + prop34**(-1)*prop127**(-1)*prop1257**(-1) * ( 1.D0/2.D0*
     &    za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6
     &    )*izb(p7,b7)*zab2(p2,p5,p7,p1)*zab3(p5,p3,p4,p6,b5) - zb(p1,
     &    p4)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p2,p7
     &    ,b7)*zab2(p6,p3,p4,b6)*zab3(p5,p1,p2,p7,b5) + zb(p1,p4)*izb(
     &    p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p2,p7,b7)*
     &    zab2(p6,p3,p4,b6)*zab3(p5,p3,p4,p6,b5) - 2*zb(p1,b5)*izb(p2,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p5,p6,p4)*zab2(
     &    p5,p2,p7,b7)*zab2(p6,p3,p4,b6) )
      b = b + prop34**(-1)*prop127**(-1) * ( 2*za(p2,p3)*za(p5,p6)*zb(
     &    p1,p4)*zb(p1,b7)*zb(b5,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7) - za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6
     &    )*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7) - za(p2,p6)*za(
     &    p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7) - za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7) - za(
     &    p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p5,p2,p7,b7) + 2*za(p5,p6)*zb(p1,p4)*zb(b5,b6
     &    )*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p2,p7,
     &    b7) )
      b = b + prop34**(-1)*prop1267**(-1)*Mwsq**(-1) * (  - 1.D0/2.D0*
     &    za(p1,p6)*za(p2,p7)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s167**(-1)*s345 + za(p1,p6)*
     &    za(p2,p7)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*s(p3,p4)*s167**(-1) + 1.D0/2.D0*za(p1,
     &    p6)*za(p2,p7)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p5,
     &    b5)*izb(p6,b6)*izb(p7,b7)*s(p3,p5)*s167**(-1) + 1.D0/2.D0*za(
     &    p1,p6)*za(p2,p7)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(
     &    p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p4,p5)*s167**(-1) + 1.D0/2.D0*
     &    za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s345 - za(p2,p6)*za(p3,p5)*
     &    zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6
     &    )*izb(p7,b7)*s(p3,p4) - 1.D0/2.D0*za(p2,p6)*za(p3,p5)*zb(p1,
     &    b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*s(p3,p5) - 1.D0/2.D0*za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(
     &    p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    s(p4,p5) )
      b = b + prop34**(-1)*prop1267**(-1)*Mwsq**(-1) * ( 1.D0/2.D0*za(
     &    p2,p7)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p6)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p1,p6,p7)*s167**(-1)
     &    *s345 - za(p2,p7)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*
     &    izb(p1,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p1,p6,p7)
     &    *s(p3,p4)*s167**(-1) - 1.D0/2.D0*za(p2,p7)*za(p3,p5)*zb(p1,b6
     &    )*zb(p1,b7)*zb(p4,b5)*izb(p1,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7
     &    ,b7)*zab2(p7,p1,p6,p7)*s(p3,p5)*s167**(-1) - 1.D0/2.D0*za(p2,
     &    p7)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p6)*izb(p5
     &    ,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p1,p6,p7)*s(p4,p5)*
     &    s167**(-1) + 1.D0/2.D0*za(p3,p5)*zb(p1,b6)*zb(p2,b7)*zb(p4,b5
     &    )*izb(p1,p6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p2,p6,p7,p1)*s345 - za(p3,p5)*zb(p1,b6)*zb(p2,b7)*zb(p4,
     &    b5)*izb(p1,p6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p2,p6,p7,p1)*s(p3,p4) - 1.D0/2.D0*za(p3,p5)*zb(p1,b6)*
     &    zb(p2,b7)*zb(p4,b5)*izb(p1,p6)*izb(p2,p7)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*zab2(p2,p6,p7,p1)*s(p3,p5) )
      b = b + prop34**(-1)*prop1267**(-1)*Mwsq**(-1) * (  - 1.D0/2.D0*
     &    za(p3,p5)*zb(p1,b6)*zb(p2,b7)*zb(p4,b5)*izb(p1,p6)*izb(p2,p7)
     &    *izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p6,p7,p1)*s(p4,p5)
     &     + 1.D0/2.D0*za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p1,p6)*izb(p2,
     &    p7)*izb(p5,b5)*izb(p6,b6)*zab2(p7,p2,p6,p1)*s345 - za(p3,p5)*
     &    zb(p1,b6)*zb(p4,b5)*izb(p1,p6)*izb(p2,p7)*izb(p5,b5)*izb(p6,
     &    b6)*zab2(p7,p2,p6,p1)*s(p3,p4) - 1.D0/2.D0*za(p3,p5)*zb(p1,b6
     &    )*zb(p4,b5)*izb(p1,p6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*zab2(
     &    p7,p2,p6,p1)*s(p3,p5) - 1.D0/2.D0*za(p3,p5)*zb(p1,b6)*zb(p4,
     &    b5)*izb(p1,p6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*zab2(p7,p2,p6
     &    ,p1)*s(p4,p5) + 1.D0/2.D0*za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(
     &    p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*
     &    s345 - za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p3,p4) - 1.D0/2.D0*
     &    za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6
     &    )*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p3,p5) )
      b = b + prop34**(-1)*prop1267**(-1)*Mwsq**(-1) * (  - 1.D0/2.D0*
     &    za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6
     &    )*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p4,p5) - 1.D0/2.D0*za(p3,p5)
     &    *zb(p1,b7)*zb(p2,b6)*zb(p4,b5)*izb(p1,p7)*izb(p2,p6)*izb(p5,
     &    b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p6,p7,p1)*s345 + za(p3,p5)*
     &    zb(p1,b7)*zb(p2,b6)*zb(p4,b5)*izb(p1,p7)*izb(p2,p6)*izb(p5,b5
     &    )*izb(p6,b6)*izb(p7,b7)*zab2(p2,p6,p7,p1)*s(p3,p4) + 1.D0/2.D0
     &    *za(p3,p5)*zb(p1,b7)*zb(p2,b6)*zb(p4,b5)*izb(p1,p7)*izb(p2,p6
     &    )*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p6,p7,p1)*s(p3,p5)
     &     + 1.D0/2.D0*za(p3,p5)*zb(p1,b7)*zb(p2,b6)*zb(p4,b5)*izb(p1,
     &    p7)*izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p6,p7
     &    ,p1)*s(p4,p5) - 1.D0/2.D0*za(p3,p5)*zb(p1,b7)*zb(p4,b5)*izb(
     &    p1,p7)*izb(p2,p6)*izb(p5,b5)*izb(p7,b7)*zab2(p6,p2,p7,p1)*
     &    s345 + za(p3,p5)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p2,p6)*
     &    izb(p5,b5)*izb(p7,b7)*zab2(p6,p2,p7,p1)*s(p3,p4) + 1.D0/2.D0*
     &    za(p3,p5)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p2,p6)*izb(p5,b5
     &    )*izb(p7,b7)*zab2(p6,p2,p7,p1)*s(p3,p5) )
      b = b + prop34**(-1)*prop1267**(-1)*Mwsq**(-1) * ( 1.D0/2.D0*za(
     &    p3,p5)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p2,p6)*izb(p5,b5)*
     &    izb(p7,b7)*zab2(p6,p2,p7,p1)*s(p4,p5) + 1.D0/2.D0*za(p3,p5)*
     &    zb(p1,b7)*zb(p4,b5)*izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*zab2(p7,p2,p6,b6)*s345 - za(p3,p5)*zb(p1,b7)*zb(p4,b5)*
     &    izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p6,b6)
     &    *s(p3,p4) - 1.D0/2.D0*za(p3,p5)*zb(p1,b7)*zb(p4,b5)*izb(p2,p6
     &    )*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p6,b6)*s(p3,p5)
     &     - 1.D0/2.D0*za(p3,p5)*zb(p1,b7)*zb(p4,b5)*izb(p2,p6)*izb(p5,
     &    b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p6,b6)*s(p4,p5) )
      b = b + prop34**(-1)*prop1267**(-1)*Mwsq**(-1)*Qsum * (  - 1.D0/2.
     &    D0*za(p1,p6)*za(p2,p7)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5
     &    )*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s167**(-1)*s345 + za(p1,p6
     &    )*za(p2,p7)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p5,b5
     &    )*izb(p6,b6)*izb(p7,b7)*s(p3,p4)*s167**(-1) + 1.D0/2.D0*za(p1
     &    ,p6)*za(p2,p7)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p5
     &    ,b5)*izb(p6,b6)*izb(p7,b7)*s(p3,p5)*s167**(-1) + 1.D0/2.D0*
     &    za(p1,p6)*za(p2,p7)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p4,p5)*s167**(-1) + 1.D0/2.
     &    D0*za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s345 - za(p2,p6)*za(p3,
     &    p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*s(p3,p4) - 1.D0/2.D0*za(p2,p6)*za(p3,p5)*
     &    zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6
     &    )*izb(p7,b7)*s(p3,p5) - 1.D0/2.D0*za(p2,p6)*za(p3,p5)*zb(p1,
     &    b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*s(p4,p5) )
      b = b + prop34**(-1)*prop1267**(-1)*Mwsq**(-1)*Qsum * ( 1.D0/2.D0
     &    *za(p2,p7)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p6)
     &    *izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p1,p6,p7)*
     &    s167**(-1)*s345 - za(p2,p7)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(
     &    p4,b5)*izb(p1,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p1
     &    ,p6,p7)*s(p3,p4)*s167**(-1) - 1.D0/2.D0*za(p2,p7)*za(p3,p5)*
     &    zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p6)*izb(p5,b5)*izb(p6,b6
     &    )*izb(p7,b7)*zab2(p7,p1,p6,p7)*s(p3,p5)*s167**(-1) - 1.D0/2.D0
     &    *za(p2,p7)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p6)
     &    *izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p1,p6,p7)*s(p4,p5)*
     &    s167**(-1) + 1.D0/2.D0*za(p3,p5)*zb(p1,b6)*zb(p2,b7)*zb(p4,b5
     &    )*izb(p1,p6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p2,p6,p7,p1)*s345 - za(p3,p5)*zb(p1,b6)*zb(p2,b7)*zb(p4,
     &    b5)*izb(p1,p6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p2,p6,p7,p1)*s(p3,p4) - 1.D0/2.D0*za(p3,p5)*zb(p1,b6)*
     &    zb(p2,b7)*zb(p4,b5)*izb(p1,p6)*izb(p2,p7)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*zab2(p2,p6,p7,p1)*s(p3,p5) )
      b = b + prop34**(-1)*prop1267**(-1)*Mwsq**(-1)*Qsum * (  - 1.D0/2.
     &    D0*za(p3,p5)*zb(p1,b6)*zb(p2,b7)*zb(p4,b5)*izb(p1,p6)*izb(p2,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p6,p7,p1)*s(p4,
     &    p5) + 1.D0/2.D0*za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p1,p6)*izb(
     &    p2,p7)*izb(p5,b5)*izb(p6,b6)*zab2(p7,p2,p6,p1)*s345 - za(p3,
     &    p5)*zb(p1,b6)*zb(p4,b5)*izb(p1,p6)*izb(p2,p7)*izb(p5,b5)*izb(
     &    p6,b6)*zab2(p7,p2,p6,p1)*s(p3,p4) - 1.D0/2.D0*za(p3,p5)*zb(p1
     &    ,b6)*zb(p4,b5)*izb(p1,p6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*
     &    zab2(p7,p2,p6,p1)*s(p3,p5) - 1.D0/2.D0*za(p3,p5)*zb(p1,b6)*
     &    zb(p4,b5)*izb(p1,p6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*zab2(p7
     &    ,p2,p6,p1)*s(p4,p5) - 1.D0/2.D0*za(p3,p5)*zb(p1,b6)*zb(p4,b5)
     &    *izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7
     &    )*s345 + za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p3,p4) + 1.D0/2.D0*
     &    za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6
     &    )*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p3,p5) )
      b = b + prop34**(-1)*prop1267**(-1)*Mwsq**(-1)*Qsum * ( 1.D0/2.D0
     &    *za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p4,p5) + 1.D0/2.D0*za(p3,
     &    p5)*zb(p1,b7)*zb(p2,b6)*zb(p4,b5)*izb(p1,p7)*izb(p2,p6)*izb(
     &    p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p6,p7,p1)*s345 - za(p3,
     &    p5)*zb(p1,b7)*zb(p2,b6)*zb(p4,b5)*izb(p1,p7)*izb(p2,p6)*izb(
     &    p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p6,p7,p1)*s(p3,p4) - 1.D0
     &    /2.D0*za(p3,p5)*zb(p1,b7)*zb(p2,b6)*zb(p4,b5)*izb(p1,p7)*izb(
     &    p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p6,p7,p1)*s(
     &    p3,p5) - 1.D0/2.D0*za(p3,p5)*zb(p1,b7)*zb(p2,b6)*zb(p4,b5)*
     &    izb(p1,p7)*izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(
     &    p2,p6,p7,p1)*s(p4,p5) + 1.D0/2.D0*za(p3,p5)*zb(p1,b7)*zb(p4,
     &    b5)*izb(p1,p7)*izb(p2,p6)*izb(p5,b5)*izb(p7,b7)*zab2(p6,p2,p7
     &    ,p1)*s345 - za(p3,p5)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p2,
     &    p6)*izb(p5,b5)*izb(p7,b7)*zab2(p6,p2,p7,p1)*s(p3,p4) - 1.D0/2.
     &    D0*za(p3,p5)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p2,p6)*izb(p5
     &    ,b5)*izb(p7,b7)*zab2(p6,p2,p7,p1)*s(p3,p5) )
      b = b + prop34**(-1)*prop1267**(-1)*Mwsq**(-1)*Qsum * (  - 1.D0/2.
     &    D0*za(p3,p5)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p2,p6)*izb(p5
     &    ,b5)*izb(p7,b7)*zab2(p6,p2,p7,p1)*s(p4,p5) - 1.D0/2.D0*za(p3,
     &    p5)*zb(p1,b7)*zb(p4,b5)*izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*zab2(p7,p2,p6,b6)*s345 + za(p3,p5)*zb(p1,b7)*zb(p4,b5)
     &    *izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p6,b6
     &    )*s(p3,p4) + 1.D0/2.D0*za(p3,p5)*zb(p1,b7)*zb(p4,b5)*izb(p2,
     &    p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p6,b6)*s(p3,
     &    p5) + 1.D0/2.D0*za(p3,p5)*zb(p1,b7)*zb(p4,b5)*izb(p2,p6)*izb(
     &    p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p6,b6)*s(p4,p5) )
      b = b + prop34**(-1)*prop1267**(-1) * ( 1.D0/2.D0*za(p1,p6)*za(p2
     &    ,p7)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*s167**(-1) + za(p1,p7)*za(p2,p3)*za(p3,p5)*
     &    zb(p1,p3)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p6)*izb(p5,b5)
     &    *izb(p6,b6)*izb(p7,b7)*s167**(-1) + za(p1,p7)*za(p2,p4)*za(p3
     &    ,p5)*zb(p1,p4)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p6)*izb(
     &    p5,b5)*izb(p6,b6)*izb(p7,b7)*s167**(-1) - za(p2,p3)*za(p3,p5)
     &    *za(p6,p7)*zb(p1,b6)*zb(p1,b7)*zb(p3,p6)*zb(p4,b5)*izb(p1,p6)
     &    *izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s167**(-1) + za(p2,p3)*za(
     &    p3,p5)*zb(p1,p3)*zb(p1,b6)*zb(p2,b7)*zb(p4,b5)*izb(p1,p6)*
     &    izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7) - za(p2,p3)*za(p3
     &    ,p5)*zb(p1,p3)*zb(p1,b7)*zb(p2,b6)*zb(p4,b5)*izb(p1,p7)*izb(
     &    p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7) - za(p2,p3)*za(p3,p5)
     &    *zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*zab2(p6,p1,p7,p3)*s167**(-1) - za(p2,p3)*zb(p1
     &    ,b6)*zb(p1,b7)*izb(p1,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p5,p3,p4,b5)*zab2(p7,p1,p6,p4)*s167**(-1) )
      b = b + prop34**(-1)*prop1267**(-1) * (  - za(p2,p3)*zb(p1,b6)*
     &    zb(p1,b7)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5
     &    ,p3,p4,b5)*zab2(p6,p1,p7,p4)*s167**(-1) - za(p2,p4)*za(p3,p5)
     &    *za(p6,p7)*zb(p1,b6)*zb(p1,b7)*zb(p4,p6)*zb(p4,b5)*izb(p1,p6)
     &    *izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s167**(-1) + za(p2,p4)*za(
     &    p3,p5)*zb(p1,p4)*zb(p1,b6)*zb(p2,b7)*zb(p4,b5)*izb(p1,p6)*
     &    izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7) - za(p2,p4)*za(p3
     &    ,p5)*zb(p1,p4)*zb(p1,b7)*zb(p2,b6)*zb(p4,b5)*izb(p1,p7)*izb(
     &    p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7) - za(p2,p4)*za(p3,p5)
     &    *zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*zab2(p6,p1,p7,p4)*s167**(-1) - za(p2,p5)*za(p3
     &    ,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,p5)*izb(p1,p6)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*zab2(p7,p1,p6,b5)*s167**(-1) - za(p2,p5)*
     &    za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,p5)*izb(p1,p7)*izb(p5,b5)
     &    *izb(p6,b6)*izb(p7,b7)*zab2(p6,p1,p7,b5)*s167**(-1) - 1.D0/2.D
     &    0*za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7
     &    )*izb(p5,b5)*izb(p6,b6)*izb(p7,b7) )
      b = b + prop34**(-1)*prop1267**(-1) * (  - 1.D0/2.D0*za(p2,p7)*
     &    za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p6)*izb(p5,b5)
     &    *izb(p6,b6)*izb(p7,b7)*zab2(p7,p1,p6,p7)*s167**(-1) + za(p3,
     &    p5)*za(p3,p6)*zb(p1,p3)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p2
     &    ,p6)*izb(p5,b5)*izb(p7,b7) - za(p3,p5)*za(p3,p7)*zb(p1,p3)*
     &    zb(p1,b6)*zb(p4,b5)*izb(p1,p6)*izb(p2,p7)*izb(p5,b5)*izb(p6,
     &    b6) + za(p3,p5)*za(p4,p6)*zb(p1,p4)*zb(p1,b7)*zb(p4,b5)*izb(
     &    p1,p7)*izb(p2,p6)*izb(p5,b5)*izb(p7,b7) - za(p3,p5)*za(p4,p7)
     &    *zb(p1,p4)*zb(p1,b6)*zb(p4,b5)*izb(p1,p6)*izb(p2,p7)*izb(p5,
     &    b5)*izb(p6,b6) - za(p3,p5)*zb(p1,p3)*zb(p4,b5)*izb(p2,p6)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p6,b6)*zab3(p3,p2
     &    ,p6,p7,b7)*s267**(-1) - za(p3,p5)*zb(p1,p3)*zb(p4,b5)*izb(p2,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*zab3(
     &    p3,p2,p6,p7,b6)*s267**(-1) - za(p3,p5)*zb(p1,p4)*zb(p4,b5)*
     &    izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p6,b6)
     &    *zab3(p4,p2,p6,p7,b7)*s267**(-1) )
      b = b + prop34**(-1)*prop1267**(-1) * (  - za(p3,p5)*zb(p1,p4)*
     &    zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6
     &    ,p2,p7,b7)*zab3(p4,p2,p6,p7,b6)*s267**(-1) - za(p3,p5)*zb(p1,
     &    b5)*zb(p1,b6)*zb(p4,p5)*izb(p1,p6)*izb(p2,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7) + za(p3,p5)*zb(p1,b5)*zb(
     &    p1,b7)*zb(p4,p5)*izb(p1,p7)*izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p5,p2,p6,b6) - za(p3,p5)*zb(p1,b5)*zb(p4,p5)*
     &    izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p6,b6)
     &    *zab3(p5,p2,p6,p7,b7)*s267**(-1) - za(p3,p5)*zb(p1,b5)*zb(p4,
     &    p5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7
     &    ,b7)*zab3(p5,p2,p6,p7,b6)*s267**(-1) - 1.D0/2.D0*za(p3,p5)*
     &    zb(p1,b6)*zb(p2,b7)*zb(p4,b5)*izb(p1,p6)*izb(p2,p7)*izb(p5,b5
     &    )*izb(p6,b6)*izb(p7,b7)*zab2(p2,p6,p7,p1) - 1.D0/2.D0*za(p3,
     &    p5)*zb(p1,b6)*zb(p4,b5)*izb(p1,p6)*izb(p2,p7)*izb(p5,b5)*izb(
     &    p6,b6)*zab2(p7,p2,p6,p1) - 1.D0/2.D0*za(p3,p5)*zb(p1,b6)*zb(
     &    p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2
     &    ,p7,b7) )
      b = b + prop34**(-1)*prop1267**(-1) * ( 1.D0/2.D0*za(p3,p5)*zb(p1
     &    ,b7)*zb(p2,b6)*zb(p4,b5)*izb(p1,p7)*izb(p2,p6)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab2(p2,p6,p7,p1) + 1.D0/2.D0*za(p3,p5)
     &    *zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p2,p6)*izb(p5,b5)*izb(p7,
     &    b7)*zab2(p6,p2,p7,p1) - 1.D0/2.D0*za(p3,p5)*zb(p1,b7)*zb(p4,
     &    b5)*izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p6
     &    ,b6) - zb(p1,p4)*zb(p1,b6)*izb(p1,p6)*izb(p2,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab2(p3,p2,p7,b7)*zab2(p5,p3,p4,b5) +
     &    zb(p1,p4)*zb(p1,b7)*izb(p1,p7)*izb(p2,p6)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*zab2(p3,p2,p6,b6)*zab2(p5,p3,p4,b5) - zb(p1,p4
     &    )*izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p3,p4,
     &    b5)*zab2(p7,p2,p6,b6)*zab3(p3,p2,p6,p7,b7)*s267**(-1) - zb(p1
     &    ,p4)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p3,
     &    p4,b5)*zab2(p6,p2,p7,b7)*zab3(p3,p2,p6,p7,b6)*s267**(-1) )
      b = b + prop34**(-1)*prop1267**(-1)*Qsum * ( 1.D0/2.D0*za(p1,p6)*
     &    za(p2,p7)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*s167**(-1) + za(p1,p7)*za(p2,p3)*za(p3,
     &    p5)*zb(p1,p3)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p6)*izb(p5
     &    ,b5)*izb(p6,b6)*izb(p7,b7)*s167**(-1) + za(p1,p7)*za(p2,p4)*
     &    za(p3,p5)*zb(p1,p4)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p6)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s167**(-1) - za(p2,p3)*za(p3
     &    ,p5)*za(p6,p7)*zb(p1,b6)*zb(p1,b7)*zb(p3,p6)*zb(p4,b5)*izb(p1
     &    ,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s167**(-1) + za(p2,p3)*
     &    za(p3,p5)*zb(p1,p3)*zb(p1,b6)*zb(p2,b7)*zb(p4,b5)*izb(p1,p6)*
     &    izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7) + za(p2,p3)*za(p3
     &    ,p5)*zb(p1,p3)*zb(p1,b7)*zb(p2,b6)*zb(p4,b5)*izb(p1,p7)*izb(
     &    p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7) - za(p2,p3)*za(p3,p5)
     &    *zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*zab2(p6,p1,p7,p3)*s167**(-1) - za(p2,p3)*zb(p1
     &    ,b6)*zb(p1,b7)*izb(p1,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p5,p3,p4,b5)*zab2(p7,p1,p6,p4)*s167**(-1) )
      b = b + prop34**(-1)*prop1267**(-1)*Qsum * (  - za(p2,p3)*zb(p1,
     &    b6)*zb(p1,b7)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p5,p3,p4,b5)*zab2(p6,p1,p7,p4)*s167**(-1) - za(p2,p4)*
     &    za(p3,p5)*za(p6,p7)*zb(p1,b6)*zb(p1,b7)*zb(p4,p6)*zb(p4,b5)*
     &    izb(p1,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s167**(-1) + za(
     &    p2,p4)*za(p3,p5)*zb(p1,p4)*zb(p1,b6)*zb(p2,b7)*zb(p4,b5)*izb(
     &    p1,p6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7) + za(p2,p4
     &    )*za(p3,p5)*zb(p1,p4)*zb(p1,b7)*zb(p2,b6)*zb(p4,b5)*izb(p1,p7
     &    )*izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7) - za(p2,p4)*za(
     &    p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab2(p6,p1,p7,p4)*s167**(-1) - za(p2,p5
     &    )*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,p5)*izb(p1,p6)*izb(p5,
     &    b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p1,p6,b5)*s167**(-1) - za(
     &    p2,p5)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,p5)*izb(p1,p7)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p1,p7,b5)*s167**(-1)
     &     - 1.D0/2.D0*za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5
     &    )*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7) )
      b = b + prop34**(-1)*prop1267**(-1)*Qsum * (  - 1.D0/2.D0*za(p2,
     &    p7)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p6)*izb(p5
     &    ,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p1,p6,p7)*s167**(-1) - za(
     &    p3,p5)*za(p3,p6)*zb(p1,p3)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*
     &    izb(p2,p6)*izb(p5,b5)*izb(p7,b7) - za(p3,p5)*za(p3,p7)*zb(p1,
     &    p3)*zb(p1,b6)*zb(p4,b5)*izb(p1,p6)*izb(p2,p7)*izb(p5,b5)*izb(
     &    p6,b6) - za(p3,p5)*za(p4,p6)*zb(p1,p4)*zb(p1,b7)*zb(p4,b5)*
     &    izb(p1,p7)*izb(p2,p6)*izb(p5,b5)*izb(p7,b7) - za(p3,p5)*za(p4
     &    ,p7)*zb(p1,p4)*zb(p1,b6)*zb(p4,b5)*izb(p1,p6)*izb(p2,p7)*izb(
     &    p5,b5)*izb(p6,b6) + za(p3,p5)*zb(p1,p3)*zb(p4,b5)*izb(p2,p6)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p6,b6)*zab3(p3,p2
     &    ,p6,p7,b7)*s267**(-1) + za(p3,p5)*zb(p1,p3)*zb(p4,b5)*izb(p2,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*zab3(
     &    p3,p2,p6,p7,b6)*s267**(-1) + za(p3,p5)*zb(p1,p4)*zb(p4,b5)*
     &    izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p6,b6)
     &    *zab3(p4,p2,p6,p7,b7)*s267**(-1) )
      b = b + prop34**(-1)*prop1267**(-1)*Qsum * ( za(p3,p5)*zb(p1,p4)*
     &    zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6
     &    ,p2,p7,b7)*zab3(p4,p2,p6,p7,b6)*s267**(-1) - za(p3,p5)*zb(p1,
     &    b5)*zb(p1,b6)*zb(p4,p5)*izb(p1,p6)*izb(p2,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7) - za(p3,p5)*zb(p1,b5)*zb(
     &    p1,b7)*zb(p4,p5)*izb(p1,p7)*izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p5,p2,p6,b6) + za(p3,p5)*zb(p1,b5)*zb(p4,p5)*
     &    izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p6,b6)
     &    *zab3(p5,p2,p6,p7,b7)*s267**(-1) + za(p3,p5)*zb(p1,b5)*zb(p4,
     &    p5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7
     &    ,b7)*zab3(p5,p2,p6,p7,b6)*s267**(-1) - 1.D0/2.D0*za(p3,p5)*
     &    zb(p1,b6)*zb(p2,b7)*zb(p4,b5)*izb(p1,p6)*izb(p2,p7)*izb(p5,b5
     &    )*izb(p6,b6)*izb(p7,b7)*zab2(p2,p6,p7,p1) - 1.D0/2.D0*za(p3,
     &    p5)*zb(p1,b6)*zb(p4,b5)*izb(p1,p6)*izb(p2,p7)*izb(p5,b5)*izb(
     &    p6,b6)*zab2(p7,p2,p6,p1) + 1.D0/2.D0*za(p3,p5)*zb(p1,b6)*zb(
     &    p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2
     &    ,p7,b7) )
      b = b + prop34**(-1)*prop1267**(-1)*Qsum * (  - 1.D0/2.D0*za(p3,
     &    p5)*zb(p1,b7)*zb(p2,b6)*zb(p4,b5)*izb(p1,p7)*izb(p2,p6)*izb(
     &    p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p6,p7,p1) - 1.D0/2.D0*
     &    za(p3,p5)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p2,p6)*izb(p5,b5
     &    )*izb(p7,b7)*zab2(p6,p2,p7,p1) + 1.D0/2.D0*za(p3,p5)*zb(p1,b7
     &    )*zb(p4,b5)*izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(
     &    p7,p2,p6,b6) - zb(p1,p4)*zb(p1,b6)*izb(p1,p6)*izb(p2,p7)*izb(
     &    p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p2,p7,b7)*zab2(p5,p3,p4,
     &    b5) - zb(p1,p4)*zb(p1,b7)*izb(p1,p7)*izb(p2,p6)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab2(p3,p2,p6,b6)*zab2(p5,p3,p4,b5) +
     &    zb(p1,p4)*izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5
     &    ,p3,p4,b5)*zab2(p7,p2,p6,b6)*zab3(p3,p2,p6,p7,b7)*s267**(-1)
     &     + zb(p1,p4)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p5,p3,p4,b5)*zab2(p6,p2,p7,b7)*zab3(p3,p2,p6,p7,b6)*
     &    s267**(-1) )
      b = b + prop34**(-1)*prop1257**(-1)*Mwsq**(-1) * (  - 1.D0/2.D0*
     &    za(p1,p5)*za(p2,p7)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s157**(-1)*s346 + za(p1,p5)*
     &    za(p2,p7)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*s(p3,p4)*s157**(-1) + 1.D0/2.D0*za(p1,
     &    p5)*za(p2,p7)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p5,
     &    b5)*izb(p6,b6)*izb(p7,b7)*s(p3,p6)*s157**(-1) + 1.D0/2.D0*za(
     &    p1,p5)*za(p2,p7)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(
     &    p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p4,p6)*s157**(-1) + 1.D0/2.D0*
     &    za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s346 - za(p2,p5)*za(p3,p6)*
     &    zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6
     &    )*izb(p7,b7)*s(p3,p4) - 1.D0/2.D0*za(p2,p5)*za(p3,p6)*zb(p1,
     &    b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*s(p3,p6) - 1.D0/2.D0*za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(
     &    p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    s(p4,p6) )
      b = b + prop34**(-1)*prop1257**(-1)*Mwsq**(-1) * ( 1.D0/2.D0*za(
     &    p2,p7)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p5)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p1,p5,p7)*s157**(-1)
     &    *s346 - za(p2,p7)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*
     &    izb(p1,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p1,p5,p7)
     &    *s(p3,p4)*s157**(-1) - 1.D0/2.D0*za(p2,p7)*za(p3,p6)*zb(p1,b5
     &    )*zb(p1,b7)*zb(p4,b6)*izb(p1,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7
     &    ,b7)*zab2(p7,p1,p5,p7)*s(p3,p6)*s157**(-1) - 1.D0/2.D0*za(p2,
     &    p7)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p5)*izb(p5
     &    ,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p1,p5,p7)*s(p4,p6)*
     &    s157**(-1) + 1.D0/2.D0*za(p3,p6)*zb(p1,b5)*zb(p2,b7)*zb(p4,b6
     &    )*izb(p1,p5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p2,p5,p7,p1)*s346 - za(p3,p6)*zb(p1,b5)*zb(p2,b7)*zb(p4,
     &    b6)*izb(p1,p5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p2,p5,p7,p1)*s(p3,p4) - 1.D0/2.D0*za(p3,p6)*zb(p1,b5)*
     &    zb(p2,b7)*zb(p4,b6)*izb(p1,p5)*izb(p2,p7)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*zab2(p2,p5,p7,p1)*s(p3,p6) )
      b = b + prop34**(-1)*prop1257**(-1)*Mwsq**(-1) * (  - 1.D0/2.D0*
     &    za(p3,p6)*zb(p1,b5)*zb(p2,b7)*zb(p4,b6)*izb(p1,p5)*izb(p2,p7)
     &    *izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p5,p7,p1)*s(p4,p6)
     &     + 1.D0/2.D0*za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p1,p5)*izb(p2,
     &    p7)*izb(p5,b5)*izb(p6,b6)*zab2(p7,p2,p5,p1)*s346 - za(p3,p6)*
     &    zb(p1,b5)*zb(p4,b6)*izb(p1,p5)*izb(p2,p7)*izb(p5,b5)*izb(p6,
     &    b6)*zab2(p7,p2,p5,p1)*s(p3,p4) - 1.D0/2.D0*za(p3,p6)*zb(p1,b5
     &    )*zb(p4,b6)*izb(p1,p5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*zab2(
     &    p7,p2,p5,p1)*s(p3,p6) - 1.D0/2.D0*za(p3,p6)*zb(p1,b5)*zb(p4,
     &    b6)*izb(p1,p5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*zab2(p7,p2,p5
     &    ,p1)*s(p4,p6) + 1.D0/2.D0*za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(
     &    p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*
     &    s346 - za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p3,p4) - 1.D0/2.D0*
     &    za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6
     &    )*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p3,p6) )
      b = b + prop34**(-1)*prop1257**(-1)*Mwsq**(-1) * (  - 1.D0/2.D0*
     &    za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6
     &    )*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p4,p6) - 1.D0/2.D0*za(p3,p6)
     &    *zb(p1,b7)*zb(p2,b5)*zb(p4,b6)*izb(p1,p7)*izb(p2,p5)*izb(p5,
     &    b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p5,p7,p1)*s346 + za(p3,p6)*
     &    zb(p1,b7)*zb(p2,b5)*zb(p4,b6)*izb(p1,p7)*izb(p2,p5)*izb(p5,b5
     &    )*izb(p6,b6)*izb(p7,b7)*zab2(p2,p5,p7,p1)*s(p3,p4) + 1.D0/2.D0
     &    *za(p3,p6)*zb(p1,b7)*zb(p2,b5)*zb(p4,b6)*izb(p1,p7)*izb(p2,p5
     &    )*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p5,p7,p1)*s(p3,p6)
     &     + 1.D0/2.D0*za(p3,p6)*zb(p1,b7)*zb(p2,b5)*zb(p4,b6)*izb(p1,
     &    p7)*izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p5,p7
     &    ,p1)*s(p4,p6) - 1.D0/2.D0*za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(
     &    p1,p7)*izb(p2,p5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,p1)*
     &    s346 + za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p2,p5)*
     &    izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,p1)*s(p3,p4) + 1.D0/2.D0*
     &    za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p2,p5)*izb(p6,b6
     &    )*izb(p7,b7)*zab2(p5,p2,p7,p1)*s(p3,p6) )
      b = b + prop34**(-1)*prop1257**(-1)*Mwsq**(-1) * ( 1.D0/2.D0*za(
     &    p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p2,p5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p5,p2,p7,p1)*s(p4,p6) + 1.D0/2.D0*za(p3,p6)*
     &    zb(p1,b7)*zb(p4,b6)*izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*zab2(p7,p2,p5,b5)*s346 - za(p3,p6)*zb(p1,b7)*zb(p4,b6)*
     &    izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p5,b5)
     &    *s(p3,p4) - 1.D0/2.D0*za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p2,p5
     &    )*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p5,b5)*s(p3,p6)
     &     - 1.D0/2.D0*za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p2,p5)*izb(p5,
     &    b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p5,b5)*s(p4,p6) )
      b = b + prop34**(-1)*prop1257**(-1)*Mwsq**(-1)*Qsum * (  - 1.D0/2.
     &    D0*za(p1,p5)*za(p2,p7)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6
     &    )*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s157**(-1)*s346 + za(p1,p5
     &    )*za(p2,p7)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p5,b5
     &    )*izb(p6,b6)*izb(p7,b7)*s(p3,p4)*s157**(-1) + 1.D0/2.D0*za(p1
     &    ,p5)*za(p2,p7)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p5
     &    ,b5)*izb(p6,b6)*izb(p7,b7)*s(p3,p6)*s157**(-1) + 1.D0/2.D0*
     &    za(p1,p5)*za(p2,p7)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p4,p6)*s157**(-1) + 1.D0/2.
     &    D0*za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s346 - za(p2,p5)*za(p3,
     &    p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*s(p3,p4) - 1.D0/2.D0*za(p2,p5)*za(p3,p6)*
     &    zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6
     &    )*izb(p7,b7)*s(p3,p6) - 1.D0/2.D0*za(p2,p5)*za(p3,p6)*zb(p1,
     &    b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*s(p4,p6) )
      b = b + prop34**(-1)*prop1257**(-1)*Mwsq**(-1)*Qsum * ( 1.D0/2.D0
     &    *za(p2,p7)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p5)
     &    *izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p1,p5,p7)*
     &    s157**(-1)*s346 - za(p2,p7)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(
     &    p4,b6)*izb(p1,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p1
     &    ,p5,p7)*s(p3,p4)*s157**(-1) - 1.D0/2.D0*za(p2,p7)*za(p3,p6)*
     &    zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p5)*izb(p5,b5)*izb(p6,b6
     &    )*izb(p7,b7)*zab2(p7,p1,p5,p7)*s(p3,p6)*s157**(-1) - 1.D0/2.D0
     &    *za(p2,p7)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p5)
     &    *izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p1,p5,p7)*s(p4,p6)*
     &    s157**(-1) + 1.D0/2.D0*za(p3,p6)*zb(p1,b5)*zb(p2,b7)*zb(p4,b6
     &    )*izb(p1,p5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p2,p5,p7,p1)*s346 - za(p3,p6)*zb(p1,b5)*zb(p2,b7)*zb(p4,
     &    b6)*izb(p1,p5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p2,p5,p7,p1)*s(p3,p4) - 1.D0/2.D0*za(p3,p6)*zb(p1,b5)*
     &    zb(p2,b7)*zb(p4,b6)*izb(p1,p5)*izb(p2,p7)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*zab2(p2,p5,p7,p1)*s(p3,p6) )
      b = b + prop34**(-1)*prop1257**(-1)*Mwsq**(-1)*Qsum * (  - 1.D0/2.
     &    D0*za(p3,p6)*zb(p1,b5)*zb(p2,b7)*zb(p4,b6)*izb(p1,p5)*izb(p2,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p5,p7,p1)*s(p4,
     &    p6) + 1.D0/2.D0*za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p1,p5)*izb(
     &    p2,p7)*izb(p5,b5)*izb(p6,b6)*zab2(p7,p2,p5,p1)*s346 - za(p3,
     &    p6)*zb(p1,b5)*zb(p4,b6)*izb(p1,p5)*izb(p2,p7)*izb(p5,b5)*izb(
     &    p6,b6)*zab2(p7,p2,p5,p1)*s(p3,p4) - 1.D0/2.D0*za(p3,p6)*zb(p1
     &    ,b5)*zb(p4,b6)*izb(p1,p5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*
     &    zab2(p7,p2,p5,p1)*s(p3,p6) - 1.D0/2.D0*za(p3,p6)*zb(p1,b5)*
     &    zb(p4,b6)*izb(p1,p5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*zab2(p7
     &    ,p2,p5,p1)*s(p4,p6) - 1.D0/2.D0*za(p3,p6)*zb(p1,b5)*zb(p4,b6)
     &    *izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7
     &    )*s346 + za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p3,p4) + 1.D0/2.D0*
     &    za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6
     &    )*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p3,p6) )
      b = b + prop34**(-1)*prop1257**(-1)*Mwsq**(-1)*Qsum * ( 1.D0/2.D0
     &    *za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p4,p6) + 1.D0/2.D0*za(p3,
     &    p6)*zb(p1,b7)*zb(p2,b5)*zb(p4,b6)*izb(p1,p7)*izb(p2,p5)*izb(
     &    p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p5,p7,p1)*s346 - za(p3,
     &    p6)*zb(p1,b7)*zb(p2,b5)*zb(p4,b6)*izb(p1,p7)*izb(p2,p5)*izb(
     &    p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p5,p7,p1)*s(p3,p4) - 1.D0
     &    /2.D0*za(p3,p6)*zb(p1,b7)*zb(p2,b5)*zb(p4,b6)*izb(p1,p7)*izb(
     &    p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p5,p7,p1)*s(
     &    p3,p6) - 1.D0/2.D0*za(p3,p6)*zb(p1,b7)*zb(p2,b5)*zb(p4,b6)*
     &    izb(p1,p7)*izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(
     &    p2,p5,p7,p1)*s(p4,p6) + 1.D0/2.D0*za(p3,p6)*zb(p1,b7)*zb(p4,
     &    b6)*izb(p1,p7)*izb(p2,p5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7
     &    ,p1)*s346 - za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p2,
     &    p5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,p1)*s(p3,p4) - 1.D0/2.
     &    D0*za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p2,p5)*izb(p6
     &    ,b6)*izb(p7,b7)*zab2(p5,p2,p7,p1)*s(p3,p6) )
      b = b + prop34**(-1)*prop1257**(-1)*Mwsq**(-1)*Qsum * (  - 1.D0/2.
     &    D0*za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p2,p5)*izb(p6
     &    ,b6)*izb(p7,b7)*zab2(p5,p2,p7,p1)*s(p4,p6) - 1.D0/2.D0*za(p3,
     &    p6)*zb(p1,b7)*zb(p4,b6)*izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*zab2(p7,p2,p5,b5)*s346 + za(p3,p6)*zb(p1,b7)*zb(p4,b6)
     &    *izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p5,b5
     &    )*s(p3,p4) + 1.D0/2.D0*za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p2,
     &    p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p5,b5)*s(p3,
     &    p6) + 1.D0/2.D0*za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p2,p5)*izb(
     &    p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p5,b5)*s(p4,p6) )
      b = b + prop34**(-1)*prop1257**(-1) * ( 1.D0/2.D0*za(p1,p5)*za(p2
     &    ,p7)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*s157**(-1) + za(p1,p7)*za(p2,p3)*za(p3,p6)*
     &    zb(p1,p3)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p5)*izb(p5,b5)
     &    *izb(p6,b6)*izb(p7,b7)*s157**(-1) + za(p1,p7)*za(p2,p4)*za(p3
     &    ,p6)*zb(p1,p4)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p5)*izb(
     &    p5,b5)*izb(p6,b6)*izb(p7,b7)*s157**(-1) - za(p2,p3)*za(p3,p6)
     &    *za(p5,p7)*zb(p1,b5)*zb(p1,b7)*zb(p3,p5)*zb(p4,b6)*izb(p1,p5)
     &    *izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s157**(-1) + za(p2,p3)*za(
     &    p3,p6)*zb(p1,p3)*zb(p1,b5)*zb(p2,b7)*zb(p4,b6)*izb(p1,p5)*
     &    izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7) - za(p2,p3)*za(p3
     &    ,p6)*zb(p1,p3)*zb(p1,b7)*zb(p2,b5)*zb(p4,b6)*izb(p1,p7)*izb(
     &    p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7) - za(p2,p3)*za(p3,p6)
     &    *zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*zab2(p5,p1,p7,p3)*s157**(-1) - za(p2,p3)*zb(p1
     &    ,b5)*zb(p1,b7)*izb(p1,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p6,p3,p4,b6)*zab2(p7,p1,p5,p4)*s157**(-1) )
      b = b + prop34**(-1)*prop1257**(-1) * (  - za(p2,p3)*zb(p1,b5)*
     &    zb(p1,b7)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5
     &    ,p1,p7,p4)*zab2(p6,p3,p4,b6)*s157**(-1) - za(p2,p4)*za(p3,p6)
     &    *za(p5,p7)*zb(p1,b5)*zb(p1,b7)*zb(p4,p5)*zb(p4,b6)*izb(p1,p5)
     &    *izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s157**(-1) + za(p2,p4)*za(
     &    p3,p6)*zb(p1,p4)*zb(p1,b5)*zb(p2,b7)*zb(p4,b6)*izb(p1,p5)*
     &    izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7) - za(p2,p4)*za(p3
     &    ,p6)*zb(p1,p4)*zb(p1,b7)*zb(p2,b5)*zb(p4,b6)*izb(p1,p7)*izb(
     &    p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7) - za(p2,p4)*za(p3,p6)
     &    *zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*zab2(p5,p1,p7,p4)*s157**(-1) - 1.D0/2.D0*za(p2
     &    ,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(
     &    p5,b5)*izb(p6,b6)*izb(p7,b7) - za(p2,p6)*za(p3,p6)*zb(p1,b5)*
     &    zb(p1,b7)*zb(p4,p6)*izb(p1,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*zab2(p7,p1,p5,b6)*s157**(-1) - za(p2,p6)*za(p3,p6)*zb(p1,
     &    b5)*zb(p1,b7)*zb(p4,p6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*zab2(p5,p1,p7,b6)*s157**(-1) )
      b = b + prop34**(-1)*prop1257**(-1) * (  - 1.D0/2.D0*za(p2,p7)*
     &    za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p5)*izb(p5,b5)
     &    *izb(p6,b6)*izb(p7,b7)*zab2(p7,p1,p5,p7)*s157**(-1) + za(p3,
     &    p5)*za(p3,p6)*zb(p1,p3)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p2
     &    ,p5)*izb(p6,b6)*izb(p7,b7) - za(p3,p6)*za(p3,p7)*zb(p1,p3)*
     &    zb(p1,b5)*zb(p4,b6)*izb(p1,p5)*izb(p2,p7)*izb(p5,b5)*izb(p6,
     &    b6) + za(p3,p6)*za(p4,p5)*zb(p1,p4)*zb(p1,b7)*zb(p4,b6)*izb(
     &    p1,p7)*izb(p2,p5)*izb(p6,b6)*izb(p7,b7) - za(p3,p6)*za(p4,p7)
     &    *zb(p1,p4)*zb(p1,b5)*zb(p4,b6)*izb(p1,p5)*izb(p2,p7)*izb(p5,
     &    b5)*izb(p6,b6) - za(p3,p6)*zb(p1,p3)*zb(p4,b6)*izb(p2,p5)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p5,b5)*zab3(p3,p2
     &    ,p5,p7,b7)*s257**(-1) - za(p3,p6)*zb(p1,p3)*zb(p4,b6)*izb(p2,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*zab3(
     &    p3,p2,p5,p7,b5)*s257**(-1) - za(p3,p6)*zb(p1,p4)*zb(p4,b6)*
     &    izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p5,b5)
     &    *zab3(p4,p2,p5,p7,b7)*s257**(-1) )
      b = b + prop34**(-1)*prop1257**(-1) * (  - za(p3,p6)*zb(p1,p4)*
     &    zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5
     &    ,p2,p7,b7)*zab3(p4,p2,p5,p7,b5)*s257**(-1) - za(p3,p6)*zb(p1,
     &    b5)*zb(p1,b6)*zb(p4,p6)*izb(p1,p5)*izb(p2,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7) - 1.D0/2.D0*za(p3,p6)*zb(
     &    p1,b5)*zb(p2,b7)*zb(p4,b6)*izb(p1,p5)*izb(p2,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab2(p2,p5,p7,p1) - 1.D0/2.D0*za(p3,p6)
     &    *zb(p1,b5)*zb(p4,b6)*izb(p1,p5)*izb(p2,p7)*izb(p5,b5)*izb(p6,
     &    b6)*zab2(p7,p2,p5,p1) - 1.D0/2.D0*za(p3,p6)*zb(p1,b5)*zb(p4,
     &    b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7
     &    ,b7) + za(p3,p6)*zb(p1,b6)*zb(p1,b7)*zb(p4,p6)*izb(p1,p7)*
     &    izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p5,b5)
     &     - za(p3,p6)*zb(p1,b6)*zb(p4,p6)*izb(p2,p5)*izb(p5,b5)*izb(p6
     &    ,b6)*izb(p7,b7)*zab2(p7,p2,p5,b5)*zab3(p6,p2,p5,p7,b7)*
     &    s257**(-1) - za(p3,p6)*zb(p1,b6)*zb(p4,p6)*izb(p2,p7)*izb(p5,
     &    b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*zab3(p6,p2,p5,p7,
     &    b5)*s257**(-1) )
      b = b + prop34**(-1)*prop1257**(-1) * ( 1.D0/2.D0*za(p3,p6)*zb(p1
     &    ,b7)*zb(p2,b5)*zb(p4,b6)*izb(p1,p7)*izb(p2,p5)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab2(p2,p5,p7,p1) + 1.D0/2.D0*za(p3,p6)
     &    *zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p2,p5)*izb(p6,b6)*izb(p7,
     &    b7)*zab2(p5,p2,p7,p1) - 1.D0/2.D0*za(p3,p6)*zb(p1,b7)*zb(p4,
     &    b6)*izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p5
     &    ,b5) - zb(p1,p4)*zb(p1,b5)*izb(p1,p5)*izb(p2,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab2(p3,p2,p7,b7)*zab2(p6,p3,p4,b6) +
     &    zb(p1,p4)*zb(p1,b7)*izb(p1,p7)*izb(p2,p5)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*zab2(p3,p2,p5,b5)*zab2(p6,p3,p4,b6) - zb(p1,p4
     &    )*izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p3,p4,
     &    b6)*zab2(p7,p2,p5,b5)*zab3(p3,p2,p5,p7,b7)*s257**(-1) - zb(p1
     &    ,p4)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,
     &    p7,b7)*zab2(p6,p3,p4,b6)*zab3(p3,p2,p5,p7,b5)*s257**(-1) )
      b = b + prop34**(-1)*prop1257**(-1)*Qsum * ( 1.D0/2.D0*za(p1,p5)*
     &    za(p2,p7)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*s157**(-1) + za(p1,p7)*za(p2,p3)*za(p3,
     &    p6)*zb(p1,p3)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p5)*izb(p5
     &    ,b5)*izb(p6,b6)*izb(p7,b7)*s157**(-1) + za(p1,p7)*za(p2,p4)*
     &    za(p3,p6)*zb(p1,p4)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p5)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s157**(-1) - za(p2,p3)*za(p3
     &    ,p6)*za(p5,p7)*zb(p1,b5)*zb(p1,b7)*zb(p3,p5)*zb(p4,b6)*izb(p1
     &    ,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s157**(-1) + za(p2,p3)*
     &    za(p3,p6)*zb(p1,p3)*zb(p1,b5)*zb(p2,b7)*zb(p4,b6)*izb(p1,p5)*
     &    izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7) + za(p2,p3)*za(p3
     &    ,p6)*zb(p1,p3)*zb(p1,b7)*zb(p2,b5)*zb(p4,b6)*izb(p1,p7)*izb(
     &    p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7) - za(p2,p3)*za(p3,p6)
     &    *zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*zab2(p5,p1,p7,p3)*s157**(-1) - za(p2,p3)*zb(p1
     &    ,b5)*zb(p1,b7)*izb(p1,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p6,p3,p4,b6)*zab2(p7,p1,p5,p4)*s157**(-1) )
      b = b + prop34**(-1)*prop1257**(-1)*Qsum * (  - za(p2,p3)*zb(p1,
     &    b5)*zb(p1,b7)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p5,p1,p7,p4)*zab2(p6,p3,p4,b6)*s157**(-1) - za(p2,p4)*
     &    za(p3,p6)*za(p5,p7)*zb(p1,b5)*zb(p1,b7)*zb(p4,p5)*zb(p4,b6)*
     &    izb(p1,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s157**(-1) + za(
     &    p2,p4)*za(p3,p6)*zb(p1,p4)*zb(p1,b5)*zb(p2,b7)*zb(p4,b6)*izb(
     &    p1,p5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7) + za(p2,p4
     &    )*za(p3,p6)*zb(p1,p4)*zb(p1,b7)*zb(p2,b5)*zb(p4,b6)*izb(p1,p7
     &    )*izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7) - za(p2,p4)*za(
     &    p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab2(p5,p1,p7,p4)*s157**(-1) - 1.D0/2.D0
     &    *za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)
     &    *izb(p5,b5)*izb(p6,b6)*izb(p7,b7) - za(p2,p6)*za(p3,p6)*zb(p1
     &    ,b5)*zb(p1,b7)*zb(p4,p6)*izb(p1,p5)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p7,p1,p5,b6)*s157**(-1) - za(p2,p6)*za(p3,p6)
     &    *zb(p1,b5)*zb(p1,b7)*zb(p4,p6)*izb(p1,p7)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*zab2(p5,p1,p7,b6)*s157**(-1) )
      b = b + prop34**(-1)*prop1257**(-1)*Qsum * (  - 1.D0/2.D0*za(p2,
     &    p7)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p5)*izb(p5
     &    ,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p1,p5,p7)*s157**(-1) - za(
     &    p3,p5)*za(p3,p6)*zb(p1,p3)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*
     &    izb(p2,p5)*izb(p6,b6)*izb(p7,b7) - za(p3,p6)*za(p3,p7)*zb(p1,
     &    p3)*zb(p1,b5)*zb(p4,b6)*izb(p1,p5)*izb(p2,p7)*izb(p5,b5)*izb(
     &    p6,b6) - za(p3,p6)*za(p4,p5)*zb(p1,p4)*zb(p1,b7)*zb(p4,b6)*
     &    izb(p1,p7)*izb(p2,p5)*izb(p6,b6)*izb(p7,b7) - za(p3,p6)*za(p4
     &    ,p7)*zb(p1,p4)*zb(p1,b5)*zb(p4,b6)*izb(p1,p5)*izb(p2,p7)*izb(
     &    p5,b5)*izb(p6,b6) + za(p3,p6)*zb(p1,p3)*zb(p4,b6)*izb(p2,p5)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p5,b5)*zab3(p3,p2
     &    ,p5,p7,b7)*s257**(-1) + za(p3,p6)*zb(p1,p3)*zb(p4,b6)*izb(p2,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*zab3(
     &    p3,p2,p5,p7,b5)*s257**(-1) + za(p3,p6)*zb(p1,p4)*zb(p4,b6)*
     &    izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p5,b5)
     &    *zab3(p4,p2,p5,p7,b7)*s257**(-1) )
      b = b + prop34**(-1)*prop1257**(-1)*Qsum * ( za(p3,p6)*zb(p1,p4)*
     &    zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5
     &    ,p2,p7,b7)*zab3(p4,p2,p5,p7,b5)*s257**(-1) - za(p3,p6)*zb(p1,
     &    b5)*zb(p1,b6)*zb(p4,p6)*izb(p1,p5)*izb(p2,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7) - 1.D0/2.D0*za(p3,p6)*zb(
     &    p1,b5)*zb(p2,b7)*zb(p4,b6)*izb(p1,p5)*izb(p2,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab2(p2,p5,p7,p1) - 1.D0/2.D0*za(p3,p6)
     &    *zb(p1,b5)*zb(p4,b6)*izb(p1,p5)*izb(p2,p7)*izb(p5,b5)*izb(p6,
     &    b6)*zab2(p7,p2,p5,p1) + 1.D0/2.D0*za(p3,p6)*zb(p1,b5)*zb(p4,
     &    b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7
     &    ,b7) - za(p3,p6)*zb(p1,b6)*zb(p1,b7)*zb(p4,p6)*izb(p1,p7)*
     &    izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p5,b5)
     &     + za(p3,p6)*zb(p1,b6)*zb(p4,p6)*izb(p2,p5)*izb(p5,b5)*izb(p6
     &    ,b6)*izb(p7,b7)*zab2(p7,p2,p5,b5)*zab3(p6,p2,p5,p7,b7)*
     &    s257**(-1) + za(p3,p6)*zb(p1,b6)*zb(p4,p6)*izb(p2,p7)*izb(p5,
     &    b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*zab3(p6,p2,p5,p7,
     &    b5)*s257**(-1) )
      b = b + prop34**(-1)*prop1257**(-1)*Qsum * (  - 1.D0/2.D0*za(p3,
     &    p6)*zb(p1,b7)*zb(p2,b5)*zb(p4,b6)*izb(p1,p7)*izb(p2,p5)*izb(
     &    p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p5,p7,p1) - 1.D0/2.D0*
     &    za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p2,p5)*izb(p6,b6
     &    )*izb(p7,b7)*zab2(p5,p2,p7,p1) + 1.D0/2.D0*za(p3,p6)*zb(p1,b7
     &    )*zb(p4,b6)*izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(
     &    p7,p2,p5,b5) - zb(p1,p4)*zb(p1,b5)*izb(p1,p5)*izb(p2,p7)*izb(
     &    p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p2,p7,b7)*zab2(p6,p3,p4,
     &    b6) - zb(p1,p4)*zb(p1,b7)*izb(p1,p7)*izb(p2,p5)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab2(p3,p2,p5,b5)*zab2(p6,p3,p4,b6) +
     &    zb(p1,p4)*izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6
     &    ,p3,p4,b6)*zab2(p7,p2,p5,b5)*zab3(p3,p2,p5,p7,b7)*s257**(-1)
     &     + zb(p1,p4)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p5,p2,p7,b7)*zab2(p6,p3,p4,b6)*zab3(p3,p2,p5,p7,b5)*
     &    s257**(-1) )
      b = b + prop34**(-1) * (  - 1.D0/2.D0*za(p2,p3)*zb(p1,b5)*zb(p1,
     &    b6)*izb(p1,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p1,p5
     &    ,b7)*zab2(p7,p2,p3,p4)*s156**(-1)*s234**(-1) - 1.D0/2.D0*za(
     &    p2,p3)*zb(p1,b5)*zb(p1,b6)*izb(p1,p6)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p5,p1,p6,b7)*zab2(p7,p2,p3,p4)*s156**(-1)*
     &    s234**(-1) - 1.D0/2.D0*za(p2,p3)*zb(p1,b5)*zb(p1,b7)*izb(p1,
     &    p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p3,p4)*zab2(
     &    p7,p1,p5,b6)*s157**(-1)*s234**(-1) - 1.D0/2.D0*za(p2,p3)*zb(
     &    p1,b5)*zb(p1,b7)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p5,p1,p7,b6)*zab2(p6,p2,p3,p4)*s157**(-1)*s234**(-1) - 1.
     &    D0/2.D0*za(p2,p3)*zb(p1,b6)*zb(p1,b7)*izb(p1,p6)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p3,p4)*zab2(p7,p1,p6,b5)*
     &    s167**(-1)*s234**(-1) - 1.D0/2.D0*za(p2,p3)*zb(p1,b6)*zb(p1,
     &    b7)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p3
     &    ,p4)*zab2(p6,p1,p7,b5)*s167**(-1)*s234**(-1) + 1.D0/2.D0*zb(
     &    p1,p4)*zb(p1,b5)*izb(p1,p5)*izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p7,p2,p6,b6)*zab3(p3,p2,p6,p7,b7)*s267**(-1)
     &     )
      b = b + prop34**(-1) * ( 1.D0/2.D0*zb(p1,p4)*zb(p1,b5)*izb(p1,p5)
     &    *izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7
     &    )*zab3(p3,p2,p6,p7,b6)*s267**(-1) + 1.D0/2.D0*zb(p1,p4)*zb(p1
     &    ,b6)*izb(p1,p6)*izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p7,p2,p5,b5)*zab3(p3,p2,p5,p7,b7)*s257**(-1) + 1.D0/2.D0
     &    *zb(p1,p4)*zb(p1,b6)*izb(p1,p6)*izb(p2,p7)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*zab3(p3,p2,p5,p7,b5)*
     &    s257**(-1) - 1.D0/2.D0*zb(p1,p4)*zb(p1,b7)*izb(p1,p7)*izb(p2,
     &    p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p5,b5)*zab3(
     &    p3,p2,p5,p6,b6)*s256**(-1) - 1.D0/2.D0*zb(p1,p4)*zb(p1,b7)*
     &    izb(p1,p7)*izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(
     &    p5,p2,p6,b6)*zab3(p3,p2,p5,p6,b5)*s256**(-1) - 1.D0/2.D0*zb(
     &    p1,p4)*izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p1
     &    ,p4,b6)*zab2(p7,p2,p5,b5)*zab3(p6,p2,p5,p7,b7)*s134**(-1)*
     &    s257**(-1) - 1.D0/2.D0*zb(p1,p4)*izb(p2,p5)*izb(p5,b5)*izb(p6
     &    ,b6)*izb(p7,b7)*zab2(p3,p1,p4,b7)*zab2(p6,p2,p5,b5)*zab3(p7,
     &    p2,p5,p6,b6)*s134**(-1)*s256**(-1) )
      b = b + prop34**(-1) * (  - 1.D0/2.D0*zb(p1,p4)*izb(p2,p6)*izb(p5
     &    ,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p1,p4,b5)*zab2(p7,p2,p6,b6
     &    )*zab3(p5,p2,p6,p7,b7)*s134**(-1)*s267**(-1) - 1.D0/2.D0*zb(
     &    p1,p4)*izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p1
     &    ,p4,b7)*zab2(p5,p2,p6,b6)*zab3(p7,p2,p5,p6,b5)*s134**(-1)*
     &    s256**(-1) - 1.D0/2.D0*zb(p1,p4)*izb(p2,p7)*izb(p5,b5)*izb(p6
     &    ,b6)*izb(p7,b7)*zab2(p3,p1,p4,b5)*zab2(p6,p2,p7,b7)*zab3(p5,
     &    p2,p6,p7,b6)*s134**(-1)*s267**(-1) - 1.D0/2.D0*zb(p1,p4)*izb(
     &    p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p1,p4,b6)*
     &    zab2(p5,p2,p7,b7)*zab3(p6,p2,p5,p7,b5)*s134**(-1)*s257**(-1)
     &     + 1.D0/2.D0*zb(p1,b5)*zb(p1,b6)*izb(p1,p5)*izb(p2,p7)*izb(p5
     &    ,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p2,p7,b7)*zab2(p6,p1,p5,p4
     &    )*s156**(-1) + 1.D0/2.D0*zb(p1,b5)*zb(p1,b6)*izb(p1,p6)*izb(
     &    p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p2,p7,b7)*
     &    zab2(p5,p1,p6,p4)*s156**(-1) - 1.D0/2.D0*zb(p1,b5)*zb(p1,b7)*
     &    izb(p1,p5)*izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(
     &    p3,p2,p6,b6)*zab2(p7,p1,p5,p4)*s157**(-1) )
      b = b + prop34**(-1) * (  - 1.D0/2.D0*zb(p1,b5)*zb(p1,b7)*izb(p1,
     &    p7)*izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p2,p6
     &    ,b6)*zab2(p5,p1,p7,p4)*s157**(-1) - 1.D0/2.D0*zb(p1,b6)*zb(p1
     &    ,b7)*izb(p1,p6)*izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p3,p2,p5,b5)*zab2(p7,p1,p6,p4)*s167**(-1) - 1.D0/2.D0*
     &    zb(p1,b6)*zb(p1,b7)*izb(p1,p7)*izb(p2,p5)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*zab2(p3,p2,p5,b5)*zab2(p6,p1,p7,p4)*s167**(-1)
     &     )
      b = b + prop34**(-1)*Qsum * (  - za(p2,p3)*zb(p1,b5)*zb(p1,b6)*
     &    izb(p1,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p1,p5,b7)
     &    *zab2(p7,p2,p3,p4)*s156**(-1)*s234**(-1) - za(p2,p3)*zb(p1,b5
     &    )*zb(p1,b6)*izb(p1,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(
     &    p5,p1,p6,b7)*zab2(p7,p2,p3,p4)*s156**(-1)*s234**(-1) - za(p2,
     &    p3)*zb(p1,b5)*zb(p1,b7)*izb(p1,p5)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*zab2(p6,p2,p3,p4)*zab2(p7,p1,p5,b6)*s157**(-1)*
     &    s234**(-1) - za(p2,p3)*zb(p1,b5)*zb(p1,b7)*izb(p1,p7)*izb(p5,
     &    b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p1,p7,b6)*zab2(p6,p2,p3,p4)
     &    *s157**(-1)*s234**(-1) - za(p2,p3)*zb(p1,b6)*zb(p1,b7)*izb(p1
     &    ,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p3,p4)*zab2(
     &    p7,p1,p6,b5)*s167**(-1)*s234**(-1) - za(p2,p3)*zb(p1,b6)*zb(
     &    p1,b7)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2
     &    ,p3,p4)*zab2(p6,p1,p7,b5)*s167**(-1)*s234**(-1) + zb(p1,p4)*
     &    zb(p1,b7)*izb(p1,p7)*izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*zab2(p6,p2,p5,b5)*zab3(p3,p2,p5,p6,b6)*s256**(-1) )
      b = b + prop34**(-1)*Qsum * ( zb(p1,p4)*zb(p1,b7)*izb(p1,p7)*izb(
     &    p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p6,b6)*
     &    zab3(p3,p2,p5,p6,b5)*s256**(-1) + zb(p1,p4)*izb(p2,p5)*izb(p5
     &    ,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p1,p4,b6)*zab2(p7,p2,p5,b5
     &    )*zab3(p6,p2,p5,p7,b7)*s134**(-1)*s257**(-1) + zb(p1,p4)*izb(
     &    p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p1,p4,b7)*
     &    zab2(p6,p2,p5,b5)*zab3(p7,p2,p5,p6,b6)*s134**(-1)*s256**(-1)
     &     + zb(p1,p4)*izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p3,p1,p4,b5)*zab2(p7,p2,p6,b6)*zab3(p5,p2,p6,p7,b7)*
     &    s134**(-1)*s267**(-1) + zb(p1,p4)*izb(p2,p6)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*zab2(p3,p1,p4,b7)*zab2(p5,p2,p6,b6)*zab3(p7
     &    ,p2,p5,p6,b5)*s134**(-1)*s256**(-1) + zb(p1,p4)*izb(p2,p7)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p1,p4,b5)*zab2(p6,p2
     &    ,p7,b7)*zab3(p5,p2,p6,p7,b6)*s134**(-1)*s267**(-1) + zb(p1,p4
     &    )*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p1,p4,
     &    b6)*zab2(p5,p2,p7,b7)*zab3(p6,p2,p5,p7,b5)*s134**(-1)*
     &    s257**(-1) )
      b = b + prop34**(-1)*Qsum * ( zb(p1,b5)*zb(p1,b6)*izb(p1,p5)*izb(
     &    p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p2,p7,b7)*
     &    zab2(p6,p1,p5,p4)*s156**(-1) + zb(p1,b5)*zb(p1,b6)*izb(p1,p6)
     &    *izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p2,p7,b7
     &    )*zab2(p5,p1,p6,p4)*s156**(-1) )
      b = b + prop34**(-1)*Qsum**2 * (  - 1.D0/2.D0*za(p2,p3)*zb(p1,b5)
     &    *zb(p1,b6)*izb(p1,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(
     &    p6,p1,p5,b7)*zab2(p7,p2,p3,p4)*s156**(-1)*s234**(-1) - 1.D0/2.
     &    D0*za(p2,p3)*zb(p1,b5)*zb(p1,b6)*izb(p1,p6)*izb(p5,b5)*izb(p6
     &    ,b6)*izb(p7,b7)*zab2(p5,p1,p6,b7)*zab2(p7,p2,p3,p4)*
     &    s156**(-1)*s234**(-1) - 1.D0/2.D0*za(p2,p3)*zb(p1,b5)*zb(p1,
     &    b7)*izb(p1,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p3
     &    ,p4)*zab2(p7,p1,p5,b6)*s157**(-1)*s234**(-1) - 1.D0/2.D0*za(
     &    p2,p3)*zb(p1,b5)*zb(p1,b7)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p5,p1,p7,b6)*zab2(p6,p2,p3,p4)*s157**(-1)*
     &    s234**(-1) - 1.D0/2.D0*za(p2,p3)*zb(p1,b6)*zb(p1,b7)*izb(p1,
     &    p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p3,p4)*zab2(
     &    p7,p1,p6,b5)*s167**(-1)*s234**(-1) - 1.D0/2.D0*za(p2,p3)*zb(
     &    p1,b6)*zb(p1,b7)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p5,p2,p3,p4)*zab2(p6,p1,p7,b5)*s167**(-1)*s234**(-1) - 1.
     &    D0/2.D0*zb(p1,p4)*zb(p1,b5)*izb(p1,p5)*izb(p2,p6)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p6,b6)*zab3(p3,p2,p6,p7,b7)*
     &    s267**(-1) )
      b = b + prop34**(-1)*Qsum**2 * (  - 1.D0/2.D0*zb(p1,p4)*zb(p1,b5)
     &    *izb(p1,p5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(
     &    p6,p2,p7,b7)*zab3(p3,p2,p6,p7,b6)*s267**(-1) - 1.D0/2.D0*zb(
     &    p1,p4)*zb(p1,b6)*izb(p1,p6)*izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p7,p2,p5,b5)*zab3(p3,p2,p5,p7,b7)*s257**(-1)
     &     - 1.D0/2.D0*zb(p1,p4)*zb(p1,b6)*izb(p1,p6)*izb(p2,p7)*izb(p5
     &    ,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*zab3(p3,p2,p5,p7
     &    ,b5)*s257**(-1) - 1.D0/2.D0*zb(p1,p4)*zb(p1,b7)*izb(p1,p7)*
     &    izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p5,b5)
     &    *zab3(p3,p2,p5,p6,b6)*s256**(-1) - 1.D0/2.D0*zb(p1,p4)*zb(p1,
     &    b7)*izb(p1,p7)*izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p5,p2,p6,b6)*zab3(p3,p2,p5,p6,b5)*s256**(-1) - 1.D0/2.D0
     &    *zb(p1,p4)*izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(
     &    p3,p1,p4,b6)*zab2(p7,p2,p5,b5)*zab3(p6,p2,p5,p7,b7)*
     &    s134**(-1)*s257**(-1) - 1.D0/2.D0*zb(p1,p4)*izb(p2,p5)*izb(p5
     &    ,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p1,p4,b7)*zab2(p6,p2,p5,b5
     &    )*zab3(p7,p2,p5,p6,b6)*s134**(-1)*s256**(-1) )
      b = b + prop34**(-1)*Qsum**2 * (  - 1.D0/2.D0*zb(p1,p4)*izb(p2,p6
     &    )*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p1,p4,b5)*zab2(p7,
     &    p2,p6,b6)*zab3(p5,p2,p6,p7,b7)*s134**(-1)*s267**(-1) - 1.D0/2.
     &    D0*zb(p1,p4)*izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p3,p1,p4,b7)*zab2(p5,p2,p6,b6)*zab3(p7,p2,p5,p6,b5)*
     &    s134**(-1)*s256**(-1) - 1.D0/2.D0*zb(p1,p4)*izb(p2,p7)*izb(p5
     &    ,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p1,p4,b5)*zab2(p6,p2,p7,b7
     &    )*zab3(p5,p2,p6,p7,b6)*s134**(-1)*s267**(-1) - 1.D0/2.D0*zb(
     &    p1,p4)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p1
     &    ,p4,b6)*zab2(p5,p2,p7,b7)*zab3(p6,p2,p5,p7,b5)*s134**(-1)*
     &    s257**(-1) + 1.D0/2.D0*zb(p1,b5)*zb(p1,b6)*izb(p1,p5)*izb(p2,
     &    p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p2,p7,b7)*zab2(
     &    p6,p1,p5,p4)*s156**(-1) + 1.D0/2.D0*zb(p1,b5)*zb(p1,b6)*izb(
     &    p1,p6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p2
     &    ,p7,b7)*zab2(p5,p1,p6,p4)*s156**(-1) + 1.D0/2.D0*zb(p1,b5)*
     &    zb(p1,b7)*izb(p1,p5)*izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*zab2(p3,p2,p6,b6)*zab2(p7,p1,p5,p4)*s157**(-1) )
      b = b + prop34**(-1)*Qsum**2 * ( 1.D0/2.D0*zb(p1,b5)*zb(p1,b7)*
     &    izb(p1,p7)*izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(
     &    p3,p2,p6,b6)*zab2(p5,p1,p7,p4)*s157**(-1) + 1.D0/2.D0*zb(p1,
     &    b6)*zb(p1,b7)*izb(p1,p6)*izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p3,p2,p5,b5)*zab2(p7,p1,p6,p4)*s167**(-1) + 1.
     &    D0/2.D0*zb(p1,b6)*zb(p1,b7)*izb(p1,p7)*izb(p2,p5)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab2(p3,p2,p5,b5)*zab2(p6,p1,p7,p4)*
     &    s167**(-1) )
      b = b + prop127**(-1)*prop1267**(-1)*Mwsq**(-1) * ( za(p2,p6)*za(
     &    p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*s345 + za(p2,p6)*za(p3,p5)*zb(p1,b6)*
     &    zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*s(p1,p3) + za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,
     &    b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p1,p4) +
     &    za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p1,p5) + za(p2,p6)*za(p3,
     &    p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*s(p2,p3) + za(p2,p6)*za(p3,p5)*zb(p1,b6)*
     &    zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*s(p2,p4) + za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,
     &    b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p2,p5) +
     &    za(p2,p6)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p3,p7) + za(p2,p6)*za(p3,
     &    p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*s(p4,p7) )
      b = b + prop127**(-1)*prop1267**(-1)*Mwsq**(-1) * ( za(p2,p6)*za(
     &    p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*s(p5,p7) + za(p3,p5)*zb(p1,b6)*zb(p4,b5
     &    )*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,
     &    b7)*s345 + za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5
     &    )*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p1,p3) + za(p3,p5
     &    )*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7
     &    ,b7)*zab2(p6,p2,p7,b7)*s(p1,p4) + za(p3,p5)*zb(p1,b6)*zb(p4,
     &    b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7
     &    ,b7)*s(p1,p5) + za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(
     &    p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p2,p3) + za(
     &    p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p2,p4) + za(p3,p5)*zb(p1,b6)*
     &    zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6
     &    ,p2,p7,b7)*s(p2,p5) + za(p3,p5)*zb(p1,b6)*zb(p4,b5)*izb(p2,p7
     &    )*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)*s(p3,p7)
     &     )
      b = b + prop127**(-1)*prop1267**(-1)*Mwsq**(-1) * ( za(p3,p5)*zb(
     &    p1,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p6,p2,p7,b7)*s(p4,p7) + za(p3,p5)*zb(p1,b6)*zb(p4,b5)*
     &    izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)
     &    *s(p5,p7) )
      b = b + prop127**(-1)*prop1267**(-1) * (  - 2*za(p2,p3)*zb(p1,p4)
     &    *zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p4,p5)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*zab3(p6,p1,p2,p7,b6) - za(p2,p6)*za(p3,p5)*zb(
     &    p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7) - 2*za(p2,p6)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1
     &    ,p7)*izb(p4,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p5,
     &    p6,p4) + 2*za(p2,p7)*za(p3,p6)*zb(p1,b7)*zb(p4,b5)*zb(p4,b6)*
     &    izb(p4,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7) - za(p3,p5)*zb(p1
     &    ,b6)*zb(p4,b5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p6,p2,p7,b7) - 2*za(p3,p6)*zb(p1,p6)*zb(p4,b5)*zb(p4,b6)
     &    *izb(p2,p7)*izb(p4,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(
     &    p6,p2,p7,b7) + 2*za(p3,p6)*zb(p1,b7)*zb(p4,b5)*zb(p4,b6)*izb(
     &    p1,p7)*izb(p4,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p6
     &    ,p7,p1) - 2*zb(p1,p4)*zb(p4,b5)*izb(p2,p7)*izb(p4,p5)*izb(p5,
     &    b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p2,p7,b7)*zab3(p6,p1,p2,p7,
     &    b6) )
      b = b + prop127**(-1)*prop1267**(-1) * (  - 2*zb(p1,b6)*zb(p4,b5)
     &    *izb(p2,p7)*izb(p4,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(
     &    p3,p5,p6,p4)*zab2(p6,p2,p7,b7) )
      b = b + prop127**(-1)*prop1257**(-1)*Mwsq**(-1) * ( za(p2,p5)*za(
     &    p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*s346 + za(p2,p5)*za(p3,p6)*zb(p1,b5)*
     &    zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*s(p1,p3) + za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,
     &    b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p1,p4) +
     &    za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p1,p6) + za(p2,p5)*za(p3,
     &    p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*s(p2,p3) + za(p2,p5)*za(p3,p6)*zb(p1,b5)*
     &    zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*s(p2,p4) + za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,
     &    b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p2,p6) +
     &    za(p2,p5)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s(p3,p7) + za(p2,p5)*za(p3,
     &    p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*s(p4,p7) )
      b = b + prop127**(-1)*prop1257**(-1)*Mwsq**(-1) * ( za(p2,p5)*za(
     &    p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*s(p6,p7) - 1.D0/2.D0*za(p2,p7)*za(p3,p6
     &    )*zb(p1,b7)*zb(p4,b6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab3(
     &    p5,p1,p2,p7,b5) - 1.D0/2.D0*za(p2,p7)*za(p3,p6)*zb(p1,b7)*zb(
     &    p4,b6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab3(p5,p3,p4,p6,b5)
     &     + 1.D0/2.D0*za(p3,p6)*zb(p1,p5)*zb(p4,b6)*izb(p2,p7)*izb(p5,
     &    b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*zab3(p5,p1,p2,p7,
     &    b5) + 1.D0/2.D0*za(p3,p6)*zb(p1,p5)*zb(p4,b6)*izb(p2,p7)*izb(
     &    p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*zab3(p5,p3,p4,
     &    p6,b5) + za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s346 + za(p3,p6)*zb(
     &    p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p5,p2,p7,b7)*s(p1,p3) + za(p3,p6)*zb(p1,b5)*zb(p4,b6)*
     &    izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)
     &    *s(p1,p4) )
      b = b + prop127**(-1)*prop1257**(-1)*Mwsq**(-1) * ( za(p3,p6)*zb(
     &    p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p5,p2,p7,b7)*s(p1,p6) + za(p3,p6)*zb(p1,b5)*zb(p4,b6)*
     &    izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)
     &    *s(p2,p3) + za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,
     &    b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p2,p4) + za(p3,
     &    p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*zab2(p5,p2,p7,b7)*s(p2,p6) + za(p3,p6)*zb(p1,b5)*zb(p4
     &    ,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,
     &    p7,b7)*s(p3,p7) + za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p4,p7)
     &     + za(p3,p6)*zb(p1,b5)*zb(p4,b6)*izb(p2,p7)*izb(p5,b5)*izb(p6
     &    ,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)*s(p6,p7) - 1.D0/2.D0*za(p3,
     &    p6)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7)*zab2(p2,p5,p7,p1)*zab3(p5,p1,p2,p7,b5) - 1.D0/2.D0*za(
     &    p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p2,p5,p7,p1)*zab3(p5,p3,p4,p6,b5) )
      b = b + prop127**(-1)*prop1257**(-1) * (  - za(p2,p3)*zb(p1,p4)*
     &    zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p4,p6)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*zab3(p5,p1,p2,p7,b5) + za(p2,p3)*zb(p1,p4)*zb(
     &    p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p4,p6)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab3(p5,p3,p4,p6,b5) - za(p2,p5)*za(p3,p6)*zb(p1,
     &    b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(p6,b6)*izb(
     &    p7,b7) - 2*za(p2,p5)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)
     &    *izb(p4,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p5,p6,p4
     &    ) + 2*za(p2,p7)*za(p3,p5)*zb(p1,b7)*zb(p4,b5)*zb(p4,b6)*izb(
     &    p4,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7) - 2*za(p3,p5)*zb(p1,
     &    p5)*zb(p4,b5)*zb(p4,b6)*izb(p2,p7)*izb(p4,p6)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7) + 2*za(p3,p5)*zb(p1,b7)*
     &    zb(p4,b5)*zb(p4,b6)*izb(p1,p7)*izb(p4,p6)*izb(p5,b5)*izb(p6,
     &    b6)*izb(p7,b7)*zab2(p2,p5,p7,p1) - za(p3,p6)*zb(p1,b5)*zb(p4,
     &    b6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7
     &    ,b7) )
      b = b + prop127**(-1)*prop1257**(-1) * (  - zb(p1,p4)*zb(p4,b6)*
     &    izb(p2,p7)*izb(p4,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(
     &    p3,p2,p7,b7)*zab3(p5,p1,p2,p7,b5) + zb(p1,p4)*zb(p4,b6)*izb(
     &    p2,p7)*izb(p4,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p2
     &    ,p7,b7)*zab3(p5,p3,p4,p6,b5) - 2*zb(p1,b5)*zb(p4,b6)*izb(p2,
     &    p7)*izb(p4,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p3,p5,p6
     &    ,p4)*zab2(p5,p2,p7,b7) )
      b = b + prop127**(-1) * (  - 2*za(p2,p3)*za(p4,p5)*zb(p1,p4)*zb(
     &    p1,b7)*zb(p4,b5)*zb(p4,b6)*izb(p1,p7)*izb(p4,p6)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*s456**(-1) + 2*za(p2,p3)*za(p4,p5)*zb(
     &    p1,p4)*zb(p2,b7)*zb(p4,b5)*zb(p4,b6)*izb(p2,p7)*izb(p4,p6)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s456**(-1) - 2*za(p2,p3)*za(
     &    p4,p6)*zb(p1,p4)*zb(p1,b7)*zb(p4,b5)*zb(p4,b6)*izb(p1,p7)*
     &    izb(p4,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s456**(-1) + 2*
     &    za(p2,p3)*za(p4,p6)*zb(p1,p4)*zb(p2,b7)*zb(p4,b5)*zb(p4,b6)*
     &    izb(p2,p7)*izb(p4,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    s456**(-1) - 2*za(p2,p3)*za(p5,p6)*zb(p1,p5)*zb(p1,b7)*zb(p4,
     &    b5)*zb(p4,b6)*izb(p1,p7)*izb(p4,p5)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*s456**(-1) + 2*za(p2,p3)*za(p5,p6)*zb(p1,p5)*zb(p2
     &    ,b7)*zb(p4,b5)*zb(p4,b6)*izb(p2,p7)*izb(p4,p5)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*s456**(-1) + 2*za(p2,p3)*za(p5,p6)*zb(
     &    p1,p6)*zb(p1,b7)*zb(p4,b5)*zb(p4,b6)*izb(p1,p7)*izb(p4,p6)*
     &    izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*s456**(-1) )
      b = b + prop127**(-1) * (  - 2*za(p2,p3)*za(p5,p6)*zb(p1,p6)*zb(
     &    p2,b7)*zb(p4,b5)*zb(p4,b6)*izb(p2,p7)*izb(p4,p6)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*s456**(-1) - 2*za(p3,p7)*za(p4,p5)*zb(
     &    p1,p4)*zb(p4,b5)*zb(p4,b6)*izb(p2,p7)*izb(p4,p6)*izb(p5,b5)*
     &    izb(p6,b6)*s456**(-1) - 2*za(p3,p7)*za(p4,p6)*zb(p1,p4)*zb(p4
     &    ,b5)*zb(p4,b6)*izb(p2,p7)*izb(p4,p5)*izb(p5,b5)*izb(p6,b6)*
     &    s456**(-1) - 2*za(p3,p7)*za(p5,p6)*zb(p1,p5)*zb(p4,b5)*zb(p4,
     &    b6)*izb(p2,p7)*izb(p4,p5)*izb(p5,b5)*izb(p6,b6)*s456**(-1) +
     &    2*za(p3,p7)*za(p5,p6)*zb(p1,p6)*zb(p4,b5)*zb(p4,b6)*izb(p2,p7
     &    )*izb(p4,p6)*izb(p5,b5)*izb(p6,b6)*s456**(-1) )
      b = b + prop1267**(-1)*Mwsq**(-1) * (  - 1.D0/2.D0*za(p1,p6)*za(
     &    p2,p7)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*s167**(-1) + 1.D0/2.D0*za(p2,p6)*za(p3,
     &    p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7) + 1.D0/2.D0*za(p2,p7)*za(p3,p5)*zb(p1,b6)*
     &    zb(p1,b7)*zb(p4,b5)*izb(p1,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*zab2(p7,p1,p6,p7)*s167**(-1) + 1.D0/2.D0*za(p3,p5)*zb(p1,
     &    b6)*zb(p2,b7)*zb(p4,b5)*izb(p1,p6)*izb(p2,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*zab2(p2,p6,p7,p1) + 1.D0/2.D0*za(p3,p5)*zb(
     &    p1,b6)*zb(p4,b5)*izb(p1,p6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*
     &    zab2(p7,p2,p6,p1) + 1.D0/2.D0*za(p3,p5)*zb(p1,b6)*zb(p4,b5)*
     &    izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)
     &     - 1.D0/2.D0*za(p3,p5)*zb(p1,b7)*zb(p2,b6)*zb(p4,b5)*izb(p1,
     &    p7)*izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p6,p7
     &    ,p1) - 1.D0/2.D0*za(p3,p5)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*
     &    izb(p2,p6)*izb(p5,b5)*izb(p7,b7)*zab2(p6,p2,p7,p1) )
      b = b + prop1267**(-1)*Mwsq**(-1) * ( 1.D0/2.D0*za(p3,p5)*zb(p1,
     &    b7)*zb(p4,b5)*izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p7,p2,p6,b6) )
      b = b + prop1267**(-1)*Mwsq**(-1)*Qsum * (  - 1.D0/2.D0*za(p1,p6)
     &    *za(p2,p7)*za(p3,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p5,b5)
     &    *izb(p6,b6)*izb(p7,b7)*s167**(-1) + 1.D0/2.D0*za(p2,p6)*za(p3
     &    ,p5)*zb(p1,b6)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7) + 1.D0/2.D0*za(p2,p7)*za(p3,p5)*zb(p1,b6)*
     &    zb(p1,b7)*zb(p4,b5)*izb(p1,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*zab2(p7,p1,p6,p7)*s167**(-1) + 1.D0/2.D0*za(p3,p5)*zb(p1,
     &    b6)*zb(p2,b7)*zb(p4,b5)*izb(p1,p6)*izb(p2,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*zab2(p2,p6,p7,p1) + 1.D0/2.D0*za(p3,p5)*zb(
     &    p1,b6)*zb(p4,b5)*izb(p1,p6)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*
     &    zab2(p7,p2,p6,p1) - 1.D0/2.D0*za(p3,p5)*zb(p1,b6)*zb(p4,b5)*
     &    izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7,b7)
     &     + 1.D0/2.D0*za(p3,p5)*zb(p1,b7)*zb(p2,b6)*zb(p4,b5)*izb(p1,
     &    p7)*izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p6,p7
     &    ,p1) + 1.D0/2.D0*za(p3,p5)*zb(p1,b7)*zb(p4,b5)*izb(p1,p7)*
     &    izb(p2,p6)*izb(p5,b5)*izb(p7,b7)*zab2(p6,p2,p7,p1) )
      b = b + prop1267**(-1)*Mwsq**(-1)*Qsum * (  - 1.D0/2.D0*za(p3,p5)
     &    *zb(p1,b7)*zb(p4,b5)*izb(p2,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*zab2(p7,p2,p6,b6) )
      b = b + prop1267**(-1) * (  - za(p2,p3)*zb(p1,b6)*zb(p1,b7)*zb(p4
     &    ,b5)*izb(p1,p6)*izb(p4,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p7,p1,p6,p4)*s167**(-1) - za(p2,p3)*zb(p1,b6)*zb(p1,b7)*
     &    zb(p4,b5)*izb(p1,p7)*izb(p4,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*zab2(p6,p1,p7,p4)*s167**(-1) - zb(p1,p4)*zb(p1,b6)*zb(p4,
     &    b5)*izb(p1,p6)*izb(p2,p7)*izb(p4,p5)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p3,p2,p7,b7) + zb(p1,p4)*zb(p1,b7)*zb(p4,b5)*
     &    izb(p1,p7)*izb(p2,p6)*izb(p4,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7
     &    ,b7)*zab2(p3,p2,p6,b6) - zb(p1,p4)*zb(p4,b5)*izb(p2,p6)*izb(
     &    p4,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p6,b6)*
     &    zab3(p3,p2,p6,p7,b7)*s267**(-1) - zb(p1,p4)*zb(p4,b5)*izb(p2,
     &    p7)*izb(p4,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2,p7
     &    ,b7)*zab3(p3,p2,p6,p7,b6)*s267**(-1) )
      b = b + prop1267**(-1)*Qsum * (  - za(p2,p3)*zb(p1,b6)*zb(p1,b7)*
     &    zb(p4,b5)*izb(p1,p6)*izb(p4,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*zab2(p7,p1,p6,p4)*s167**(-1) - za(p2,p3)*zb(p1,b6)*zb(p1,
     &    b7)*zb(p4,b5)*izb(p1,p7)*izb(p4,p5)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p6,p1,p7,p4)*s167**(-1) - zb(p1,p4)*zb(p1,b6)
     &    *zb(p4,b5)*izb(p1,p6)*izb(p2,p7)*izb(p4,p5)*izb(p5,b5)*izb(p6
     &    ,b6)*izb(p7,b7)*zab2(p3,p2,p7,b7) - zb(p1,p4)*zb(p1,b7)*zb(p4
     &    ,b5)*izb(p1,p7)*izb(p2,p6)*izb(p4,p5)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p3,p2,p6,b6) + zb(p1,p4)*zb(p4,b5)*izb(p2,p6)
     &    *izb(p4,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p6,b6
     &    )*zab3(p3,p2,p6,p7,b7)*s267**(-1) + zb(p1,p4)*zb(p4,b5)*izb(
     &    p2,p7)*izb(p4,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p6,p2
     &    ,p7,b7)*zab3(p3,p2,p6,p7,b6)*s267**(-1) )
      b = b + prop1257**(-1)*Mwsq**(-1) * (  - 1.D0/2.D0*za(p1,p5)*za(
     &    p2,p7)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p5,b5)*
     &    izb(p6,b6)*izb(p7,b7)*s157**(-1) + 1.D0/2.D0*za(p2,p5)*za(p3,
     &    p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7) + 1.D0/2.D0*za(p2,p7)*za(p3,p6)*zb(p1,b5)*
     &    zb(p1,b7)*zb(p4,b6)*izb(p1,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*zab2(p7,p1,p5,p7)*s157**(-1) + 1.D0/2.D0*za(p3,p6)*zb(p1,
     &    b5)*zb(p2,b7)*zb(p4,b6)*izb(p1,p5)*izb(p2,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*zab2(p2,p5,p7,p1) + 1.D0/2.D0*za(p3,p6)*zb(
     &    p1,b5)*zb(p4,b6)*izb(p1,p5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*
     &    zab2(p7,p2,p5,p1) + 1.D0/2.D0*za(p3,p6)*zb(p1,b5)*zb(p4,b6)*
     &    izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)
     &     - 1.D0/2.D0*za(p3,p6)*zb(p1,b7)*zb(p2,b5)*zb(p4,b6)*izb(p1,
     &    p7)*izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p5,p7
     &    ,p1) - 1.D0/2.D0*za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*
     &    izb(p2,p5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,p1) )
      b = b + prop1257**(-1)*Mwsq**(-1) * ( 1.D0/2.D0*za(p3,p6)*zb(p1,
     &    b7)*zb(p4,b6)*izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p7,p2,p5,b5) )
      b = b + prop1257**(-1)*Mwsq**(-1)*Qsum * (  - 1.D0/2.D0*za(p1,p5)
     &    *za(p2,p7)*za(p3,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p5,b5)
     &    *izb(p6,b6)*izb(p7,b7)*s157**(-1) + 1.D0/2.D0*za(p2,p5)*za(p3
     &    ,p6)*zb(p1,b5)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7) + 1.D0/2.D0*za(p2,p7)*za(p3,p6)*zb(p1,b5)*
     &    zb(p1,b7)*zb(p4,b6)*izb(p1,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*zab2(p7,p1,p5,p7)*s157**(-1) + 1.D0/2.D0*za(p3,p6)*zb(p1,
     &    b5)*zb(p2,b7)*zb(p4,b6)*izb(p1,p5)*izb(p2,p7)*izb(p5,b5)*izb(
     &    p6,b6)*izb(p7,b7)*zab2(p2,p5,p7,p1) + 1.D0/2.D0*za(p3,p6)*zb(
     &    p1,b5)*zb(p4,b6)*izb(p1,p5)*izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*
     &    zab2(p7,p2,p5,p1) - 1.D0/2.D0*za(p3,p6)*zb(p1,b5)*zb(p4,b6)*
     &    izb(p2,p7)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,b7)
     &     + 1.D0/2.D0*za(p3,p6)*zb(p1,b7)*zb(p2,b5)*zb(p4,b6)*izb(p1,
     &    p7)*izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p2,p5,p7
     &    ,p1) + 1.D0/2.D0*za(p3,p6)*zb(p1,b7)*zb(p4,b6)*izb(p1,p7)*
     &    izb(p2,p5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7,p1) )
      b = b + prop1257**(-1)*Mwsq**(-1)*Qsum * (  - 1.D0/2.D0*za(p3,p6)
     &    *zb(p1,b7)*zb(p4,b6)*izb(p2,p5)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*zab2(p7,p2,p5,b5) )
      b = b + prop1257**(-1) * (  - za(p2,p3)*zb(p1,b5)*zb(p1,b7)*zb(p4
     &    ,b6)*izb(p1,p5)*izb(p4,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*
     &    zab2(p7,p1,p5,p4)*s157**(-1) - za(p2,p3)*zb(p1,b5)*zb(p1,b7)*
     &    zb(p4,b6)*izb(p1,p7)*izb(p4,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*zab2(p5,p1,p7,p4)*s157**(-1) - zb(p1,p4)*zb(p1,b5)*zb(p4,
     &    b6)*izb(p1,p5)*izb(p2,p7)*izb(p4,p6)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p3,p2,p7,b7) + zb(p1,p4)*zb(p1,b7)*zb(p4,b6)*
     &    izb(p1,p7)*izb(p2,p5)*izb(p4,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7
     &    ,b7)*zab2(p3,p2,p5,b5) - zb(p1,p4)*zb(p4,b6)*izb(p2,p5)*izb(
     &    p4,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p5,b5)*
     &    zab3(p3,p2,p5,p7,b7)*s257**(-1) - zb(p1,p4)*zb(p4,b6)*izb(p2,
     &    p7)*izb(p4,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2,p7
     &    ,b7)*zab3(p3,p2,p5,p7,b5)*s257**(-1) )
      b = b + prop1257**(-1)*Qsum * (  - za(p2,p3)*zb(p1,b5)*zb(p1,b7)*
     &    zb(p4,b6)*izb(p1,p5)*izb(p4,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,
     &    b7)*zab2(p7,p1,p5,p4)*s157**(-1) - za(p2,p3)*zb(p1,b5)*zb(p1,
     &    b7)*zb(p4,b6)*izb(p1,p7)*izb(p4,p6)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p5,p1,p7,p4)*s157**(-1) - zb(p1,p4)*zb(p1,b5)
     &    *zb(p4,b6)*izb(p1,p5)*izb(p2,p7)*izb(p4,p6)*izb(p5,b5)*izb(p6
     &    ,b6)*izb(p7,b7)*zab2(p3,p2,p7,b7) - zb(p1,p4)*zb(p1,b7)*zb(p4
     &    ,b6)*izb(p1,p7)*izb(p2,p5)*izb(p4,p6)*izb(p5,b5)*izb(p6,b6)*
     &    izb(p7,b7)*zab2(p3,p2,p5,b5) + zb(p1,p4)*zb(p4,b6)*izb(p2,p5)
     &    *izb(p4,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p7,p2,p5,b5
     &    )*zab3(p3,p2,p5,p7,b7)*s257**(-1) + zb(p1,p4)*zb(p4,b6)*izb(
     &    p2,p7)*izb(p4,p6)*izb(p5,b5)*izb(p6,b6)*izb(p7,b7)*zab2(p5,p2
     &    ,p7,b7)*zab3(p3,p2,p5,p7,b5)*s257**(-1) )
      return
      end
