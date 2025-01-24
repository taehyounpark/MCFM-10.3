!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine A6Wpgamgam(p1,p2,p3,p4,p5,p6,za,zb,b)
      implicit none
c---- Matrix element for Wgamgamma radiation from line 12
c  u(-p1) +dbar(-p2) --> ve(p3)+e^+(p4)+gamgam(p5)+g(p6)
c  including radiation from decay electron

c Overall factor of  i_*e^2*gw^2
c     including QED type propagators
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'cplx.h'
      integer::p1,p2,p3,p4,p5,p6
      real(dp):: s3,s12,s34,s125,s126,s156,s256,s456,Qsum
      complex(dp):: iza,izb,zab2,zab3,b(2,2),Mwsq,
     & prop34,prop12,prop125,prop126
c      amplitude b(h5,h6)

      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p1,p3)
      zab2(p1,p2,p3,p4)=+za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zab3(p1,p2,p3,p4,p5)=+za(p1,p2)*zb(p2,p5)+za(p1,p3)*zb(p3,p5)+za(p1,p4)*zb(p4,p5)
      iza(p1,p2)=cone/za(p1,p2)
      izb(p1,p2)=cone/zb(p1,p2)

      Qsum=1._dp/3._dp
      s156=s3(p1,p5,p6)
      s126=s3(p1,p2,p6)
      s125=s3(p1,p2,p5)
      s256=s3(p2,p5,p6)
      s456=s3(p4,p5,p6)
      s12=s(p1,p2)
      s34=s(p3,p4)
      Mwsq=cplx2(wmass**2,-wmass*wwidth)
      prop34=cplx1(s34)-Mwsq
      prop12=cplx1(s12)-Mwsq
      prop125=cplx1(s125)-Mwsq
      prop126=cplx1(s126)-Mwsq

c Commented-out gauge check

c      do b5m=p1,p2,p2-p1
c      do b6m=p1,p2,p2-p1
c      b(1,1)= - 2*za(p2,p3)*za(p5,p6)*zb(p1,p4)*zb(b5m,b6m)*izb(p5,b5m)
c     &    *izb(p6,b6m)*prop34**(-1)*prop12**(-1) + za(p2,p3)*zb(p1,p4)*
c     &    zb(p1,b5m)*zb(p4,b6m)*izb(p1,p5)*izb(p4,p6)*izb(p5,b5m)*izb(
c     &    p6,b6m)*prop125**(-1) + za(p2,p3)*zb(p1,p4)*zb(p1,b5m)*zb(p4,
c     &    b6m)*izb(p1,p5)*izb(p4,p6)*izb(p5,b5m)*izb(p6,b6m)*
c     &    prop125**(-1)*Qsum + za(p2,p3)*zb(p1,p4)*zb(p1,b5m)*izb(p1,p5
c     &    )*izb(p5,b5m)*izb(p6,b6m)*zab2(p6,p3,p4,b6m)*prop34**(-1)*
c     &    prop125**(-1) + za(p2,p3)*zb(p1,p4)*zb(p1,b5m)*izb(p1,p5)*
c     &    izb(p5,b5m)*izb(p6,b6m)*zab2(p6,p3,p4,b6m)*prop34**(-1)*
c     &    prop125**(-1)*Qsum + za(p2,p3)*zb(p1,p4)*zb(p1,b6m)*zb(p4,b5m
c     &    )*izb(p1,p6)*izb(p4,p5)*izb(p5,b5m)*izb(p6,b6m)*prop126**(-1)
c     &     + za(p2,p3)*zb(p1,p4)*zb(p1,b6m)*zb(p4,b5m)*izb(p1,p6)*izb(
c     &    p4,p5)*izb(p5,b5m)*izb(p6,b6m)*prop126**(-1)*Qsum + za(p2,p3)
c     &    *zb(p1,p4)*zb(p1,b6m)*izb(p1,p6)*izb(p5,b5m)*izb(p6,b6m)*
c     &    zab2(p5,p3,p4,b5m)*prop34**(-1)*prop126**(-1) + za(p2,p3)*zb(
c     &    p1,p4)*zb(p1,b6m)*izb(p1,p6)*izb(p5,b5m)*izb(p6,b6m)*zab2(p5,
c     &    p3,p4,b5m)*prop34**(-1)*prop126**(-1)*Qsum
c      b(1,1) = b(1,1) + 2*za(p2,p3)*zb(p1,p4)*zb(p4,b5m)*izb(p4,p5)*
c     & izb(p5,b5m)*izb(p6,b6m)*zab2(p6,p1,p2,b6m)*prop12**(-1)*
c     & prop126**(-1) + 2*za(p2,p3)*zb(p1,p4)*zb(p4,b6m)*izb(p4,p6)*izb(
c     &    p5,b5m)*izb(p6,b6m)*zab2(p5,p1,p2,b5m)*prop12**(-1)*
c     &    prop125**(-1) + 2*za(p2,p3)*zb(p1,p4)*izb(p5,b5m)*izb(p6,b6m)
c     &    *zab2(p5,p1,p2,b5m)*zab2(p6,p3,p4,b6m)*prop34**(-1)*
c     &    prop12**(-1)*prop125**(-1) + 2*za(p2,p3)*zb(p1,p4)*izb(p5,b5m
c     &    )*izb(p6,b6m)*zab2(p5,p3,p4,b5m)*zab2(p6,p1,p2,b6m)*
c     &    prop34**(-1)*prop12**(-1)*prop126**(-1) - 1.D0/2.D0*za(p2,p3)
c     &    *zb(p1,b5m)*zb(p1,b6m)*izb(p1,p5)*izb(p5,b5m)*izb(p6,b6m)*
c     &    zab2(p6,p1,p5,p4)*prop34**(-1)*s156**(-1) - za(p2,p3)*zb(p1,
c     &    b5m)*zb(p1,b6m)*izb(p1,p5)*izb(p5,b5m)*izb(p6,b6m)*zab2(p6,p1
c     &    ,p5,p4)*prop34**(-1)*Qsum*s156**(-1) - 1.D0/2.D0*za(p2,p3)*
c     &    zb(p1,b5m)*zb(p1,b6m)*izb(p1,p5)*izb(p5,b5m)*izb(p6,b6m)*
c     &    zab2(p6,p1,p5,p4)*prop34**(-1)*Qsum**2*s156**(-1) - 1.D0/2.D0
c     &    *za(p2,p3)*zb(p1,b5m)*zb(p1,b6m)*izb(p1,p6)*izb(p5,b5m)*izb(
c     &    p6,b6m)*zab2(p5,p1,p6,p4)*prop34**(-1)*s156**(-1)
c      b(1,1) = b(1,1) - za(p2,p3)*zb(p1,b5m)*zb(p1,b6m)*izb(p1,p6)*izb(
c     & p5,b5m)*izb(p6,b6m)*zab2(p5,p1,p6,p4)*prop34**(-1)*Qsum*
c     & s156**(-1) - 1.D0/2.D0*za(p2,p3)*zb(p1,b5m)*zb(p1,b6m)*izb(p1,p6
c     &    )*izb(p5,b5m)*izb(p6,b6m)*zab2(p5,p1,p6,p4)*prop34**(-1)*
c     &    Qsum**2*s156**(-1) + 2*za(p2,p3)*zb(p4,b5m)*zb(p4,b6m)*izb(p4
c     &    ,p5)*izb(p5,b5m)*izb(p6,b6m)*zab2(p6,p4,p5,p1)*prop12**(-1)*
c     &    s456**(-1) + 2*za(p2,p3)*zb(p4,b5m)*zb(p4,b6m)*izb(p4,p6)*
c     &    izb(p5,b5m)*izb(p6,b6m)*zab2(p5,p4,p6,p1)*prop12**(-1)*
c     &    s456**(-1) + 2*za(p2,p5)*za(p3,p5)*zb(p1,p5)*zb(p4,b5m)*zb(p4
c     &    ,b6m)*izb(p4,p6)*izb(p5,b5m)*izb(p6,b6m)*prop12**(-1)*
c     &    prop125**(-1) + 2*za(p2,p5)*za(p3,p5)*zb(p1,p5)*zb(p4,b5m)*
c     &    izb(p5,b5m)*izb(p6,b6m)*zab2(p6,p3,p4,b6m)*prop34**(-1)*
c     &    prop12**(-1)*prop125**(-1) + za(p2,p5)*za(p3,p5)*zb(p1,b5m)*
c     &    zb(p1,b6m)*zb(p4,p5)*izb(p1,p6)*izb(p5,b5m)*izb(p6,b6m)*
c     &    prop34**(-1)*prop126**(-1) + za(p2,p5)*za(p3,p5)*zb(p1,b5m)*
c     &    zb(p1,b6m)*zb(p4,p5)*izb(p1,p6)*izb(p5,b5m)*izb(p6,b6m)*
c     &    prop34**(-1)*prop126**(-1)*Qsum
c      b(1,1) = b(1,1) + 2*za(p2,p5)*za(p3,p5)*zb(p1,b5m)*zb(p4,p5)*izb(
c     & p5,b5m)*izb(p6,b6m)*zab2(p6,p1,p2,b6m)*prop34**(-1)*prop12**(-1)
c     & *prop126**(-1) + 2*za(p2,p5)*za(p3,p6)*za(p5,p6)*zb(p1,p5)*zb(p4
c     &    ,p6)*zb(b5m,b6m)*izb(p5,b5m)*izb(p6,b6m)*prop34**(-1)*
c     &    prop12**(-1)*prop125**(-1) - 2*za(p2,p5)*za(p3,p6)*zb(p1,p5)*
c     &    zb(p4,b6m)*izb(p5,b5m)*izb(p6,b6m)*zab2(p5,p3,p4,b5m)*
c     &    prop34**(-1)*prop12**(-1)*prop125**(-1) - 2*za(p2,p5)*za(p3,
c     &    p6)*zb(p1,b5m)*zb(p4,p6)*izb(p5,b5m)*izb(p6,b6m)*zab2(p6,p1,
c     &    p2,b6m)*prop34**(-1)*prop12**(-1)*prop125**(-1) - za(p2,p5)*
c     &    za(p3,p6)*zb(p1,b5m)*zb(p4,b6m)*izb(p5,b5m)*izb(p6,b6m)*
c     &    prop34**(-1)*prop12**(-1)*prop125**(-1)*Mwsq**(-1)*s34*s12 -
c     &    2*za(p2,p5)*za(p3,p6)*zb(p1,b5m)*zb(p4,b6m)*izb(p5,b5m)*izb(
c     &    p6,b6m)*prop34**(-1)*prop12**(-1)*prop125**(-1)*s56 + za(p2,
c     &    p5)*za(p3,p6)*zb(p1,b5m)*zb(p4,b6m)*izb(p5,b5m)*izb(p6,b6m)*
c     &    prop34**(-1)*prop12**(-1)*prop125**(-1)*s12 + za(p2,p5)*za(p3
c     &    ,p6)*zb(p1,b5m)*zb(p4,b6m)*izb(p5,b5m)*izb(p6,b6m)*
c     &    prop34**(-1)*prop12**(-1)*prop125**(-1)*s34
c      b(1,1) = b(1,1) - za(p2,p5)*za(p3,p6)*zb(p1,b5m)*zb(p4,b6m)*izb(
c     & p5,b5m)*izb(p6,b6m)*prop34**(-1)*prop12**(-1)*prop125**(-1)*s125
c     &     + za(p2,p5)*za(p3,p6)*zb(p1,b5m)*zb(p4,b6m)*izb(p5,b5m)*izb(
c     &    p6,b6m)*prop34**(-1)*prop12**(-1) + za(p2,p5)*za(p3,p6)*zb(p1
c     &    ,b5m)*zb(p4,b6m)*izb(p5,b5m)*izb(p6,b6m)*prop34**(-1)*
c     &    prop125**(-1)*Mwsq**(-1)*s34 + za(p2,p5)*za(p3,p6)*zb(p1,b5m)
c     &    *zb(p4,b6m)*izb(p5,b5m)*izb(p6,b6m)*prop34**(-1)*
c     &    prop125**(-1)*Qsum + za(p2,p5)*za(p3,p6)*zb(p1,b5m)*zb(p4,b6m
c     &    )*izb(p5,b5m)*izb(p6,b6m)*prop12**(-1)*prop125**(-1)*
c     &    Mwsq**(-1)*s12 + za(p2,p5)*za(p3,p6)*zb(p1,b5m)*zb(p4,b6m)*
c     &    izb(p5,b5m)*izb(p6,b6m)*prop12**(-1)*prop125**(-1) - za(p2,p5
c     &    )*za(p3,p6)*zb(p1,b5m)*zb(p4,b6m)*izb(p5,b5m)*izb(p6,b6m)*
c     &    prop125**(-1)*Mwsq**(-1) - 2*za(p2,p5)*zb(p1,b5m)*zb(p4,b6m)*
c     &    izb(p4,p6)*izb(p5,b5m)*izb(p6,b6m)*zab2(p3,p1,p2,p4)*
c     &    prop12**(-1)*prop125**(-1) - 2*za(p2,p5)*zb(p1,b5m)*izb(p5,
c     &    b5m)*izb(p6,b6m)*zab2(p3,p1,p2,p4)*zab2(p6,p3,p4,b6m)*
c     &    prop34**(-1)*prop12**(-1)*prop125**(-1)
c      b(1,1) = b(1,1) + 2*za(p2,p6)*za(p3,p5)*za(p5,p6)*zb(p1,p6)*zb(p4
c     & ,p5)*zb(b5m,b6m)*izb(p5,b5m)*izb(p6,b6m)*prop34**(-1)*
c     & prop12**(-1)*prop126**(-1) - 2*za(p2,p6)*za(p3,p5)*zb(p1,p6)*zb(
c     &    p4,b5m)*izb(p5,b5m)*izb(p6,b6m)*zab2(p6,p3,p4,b6m)*
c     &    prop34**(-1)*prop12**(-1)*prop126**(-1) - 2*za(p2,p6)*za(p3,
c     &    p5)*zb(p1,b6m)*zb(p4,p5)*izb(p5,b5m)*izb(p6,b6m)*zab2(p5,p1,
c     &    p2,b5m)*prop34**(-1)*prop12**(-1)*prop126**(-1) - za(p2,p6)*
c     &    za(p3,p5)*zb(p1,b6m)*zb(p4,b5m)*izb(p5,b5m)*izb(p6,b6m)*
c     &    prop34**(-1)*prop12**(-1)*prop126**(-1)*Mwsq**(-1)*s34*s12 -
c     &    2*za(p2,p6)*za(p3,p5)*zb(p1,b6m)*zb(p4,b5m)*izb(p5,b5m)*izb(
c     &    p6,b6m)*prop34**(-1)*prop12**(-1)*prop126**(-1)*s56 + za(p2,
c     &    p6)*za(p3,p5)*zb(p1,b6m)*zb(p4,b5m)*izb(p5,b5m)*izb(p6,b6m)*
c     &    prop34**(-1)*prop12**(-1)*prop126**(-1)*s12 + za(p2,p6)*za(p3
c     &    ,p5)*zb(p1,b6m)*zb(p4,b5m)*izb(p5,b5m)*izb(p6,b6m)*
c     &    prop34**(-1)*prop12**(-1)*prop126**(-1)*s34 - za(p2,p6)*za(p3
c     &    ,p5)*zb(p1,b6m)*zb(p4,b5m)*izb(p5,b5m)*izb(p6,b6m)*
c     &    prop34**(-1)*prop12**(-1)*prop126**(-1)*s126
c      b(1,1) = b(1,1) + za(p2,p6)*za(p3,p5)*zb(p1,b6m)*zb(p4,b5m)*izb(
c     & p5,b5m)*izb(p6,b6m)*prop34**(-1)*prop12**(-1) + za(p2,p6)*za(p3,
c     &    p5)*zb(p1,b6m)*zb(p4,b5m)*izb(p5,b5m)*izb(p6,b6m)*
c     &    prop34**(-1)*prop126**(-1)*Mwsq**(-1)*s34 + za(p2,p6)*za(p3,
c     &    p5)*zb(p1,b6m)*zb(p4,b5m)*izb(p5,b5m)*izb(p6,b6m)*
c     &    prop34**(-1)*prop126**(-1)*Qsum + za(p2,p6)*za(p3,p5)*zb(p1,
c     &    b6m)*zb(p4,b5m)*izb(p5,b5m)*izb(p6,b6m)*prop12**(-1)*
c     &    prop126**(-1)*Mwsq**(-1)*s12 + za(p2,p6)*za(p3,p5)*zb(p1,b6m)
c     &    *zb(p4,b5m)*izb(p5,b5m)*izb(p6,b6m)*prop12**(-1)*
c     &    prop126**(-1) - za(p2,p6)*za(p3,p5)*zb(p1,b6m)*zb(p4,b5m)*
c     &    izb(p5,b5m)*izb(p6,b6m)*prop126**(-1)*Mwsq**(-1) + 2*za(p2,p6
c     &    )*za(p3,p6)*zb(p1,p6)*zb(p4,b5m)*zb(p4,b6m)*izb(p4,p5)*izb(p5
c     &    ,b5m)*izb(p6,b6m)*prop12**(-1)*prop126**(-1) + 2*za(p2,p6)*
c     &    za(p3,p6)*zb(p1,p6)*zb(p4,b6m)*izb(p5,b5m)*izb(p6,b6m)*zab2(
c     &    p5,p3,p4,b5m)*prop34**(-1)*prop12**(-1)*prop126**(-1) + za(p2
c     &    ,p6)*za(p3,p6)*zb(p1,b5m)*zb(p1,b6m)*zb(p4,p6)*izb(p1,p5)*
c     &    izb(p5,b5m)*izb(p6,b6m)*prop34**(-1)*prop125**(-1)
c      b(1,1) = b(1,1) + za(p2,p6)*za(p3,p6)*zb(p1,b5m)*zb(p1,b6m)*zb(p4
c     & ,p6)*izb(p1,p5)*izb(p5,b5m)*izb(p6,b6m)*prop34**(-1)*
c     & prop125**(-1)*Qsum + 2*za(p2,p6)*za(p3,p6)*zb(p1,b6m)*zb(p4,p6)*
c     &    izb(p5,b5m)*izb(p6,b6m)*zab2(p5,p1,p2,b5m)*prop34**(-1)*
c     &    prop12**(-1)*prop125**(-1) - 2*za(p2,p6)*zb(p1,b6m)*zb(p4,b5m
c     &    )*izb(p4,p5)*izb(p5,b5m)*izb(p6,b6m)*zab2(p3,p1,p2,p4)*
c     &    prop12**(-1)*prop126**(-1) - 2*za(p2,p6)*zb(p1,b6m)*izb(p5,
c     &    b5m)*izb(p6,b6m)*zab2(p3,p1,p2,p4)*zab2(p5,p3,p4,b5m)*
c     &    prop34**(-1)*prop12**(-1)*prop126**(-1) + za(p3,p5)*zb(p1,p5)
c     &    *zb(p4,b5m)*izb(p2,p6)*izb(p5,b5m)*izb(p6,b6m)*zab2(p5,p2,p6,
c     &    b6m)*prop34**(-1)*prop126**(-1) - za(p3,p5)*zb(p1,p5)*zb(p4,
c     &    b5m)*izb(p2,p6)*izb(p5,b5m)*izb(p6,b6m)*zab2(p5,p2,p6,b6m)*
c     &    prop34**(-1)*prop126**(-1)*Qsum - za(p3,p5)*zb(p1,b5m)*zb(p4,
c     &    p5)*izb(p2,p6)*izb(p5,b5m)*izb(p6,b6m)*zab2(p5,p2,p6,b6m)*
c     &    prop34**(-1)*prop126**(-1) + za(p3,p5)*zb(p1,b5m)*zb(p4,p5)*
c     &    izb(p2,p6)*izb(p5,b5m)*izb(p6,b6m)*zab2(p5,p2,p6,b6m)*
c     &    prop34**(-1)*prop126**(-1)*Qsum
c      b(1,1) = b(1,1) - za(p3,p5)*zb(p1,b6m)*zb(p4,b5m)*izb(p1,p6)*izb(
c     & p5,b5m)*izb(p6,b6m)*zab2(p2,p3,p4,p1)*prop34**(-1)*prop126**(-1)
c     &     - za(p3,p5)*zb(p1,b6m)*zb(p4,b5m)*izb(p1,p6)*izb(p5,b5m)*
c     &    izb(p6,b6m)*zab2(p2,p3,p4,p1)*prop34**(-1)*prop126**(-1)*Qsum
c     &     - 2*za(p3,p5)*zb(p4,b5m)*izb(p5,b5m)*izb(p6,b6m)*zab2(p2,p3,
c     &    p4,p1)*zab2(p6,p1,p2,b6m)*prop34**(-1)*prop12**(-1)*
c     &    prop126**(-1) + za(p3,p6)*zb(p1,p6)*zb(p4,b6m)*izb(p2,p5)*
c     &    izb(p5,b5m)*izb(p6,b6m)*zab2(p6,p2,p5,b5m)*prop34**(-1)*
c     &    prop125**(-1) - za(p3,p6)*zb(p1,p6)*zb(p4,b6m)*izb(p2,p5)*
c     &    izb(p5,b5m)*izb(p6,b6m)*zab2(p6,p2,p5,b5m)*prop34**(-1)*
c     &    prop125**(-1)*Qsum - za(p3,p6)*zb(p1,b5m)*zb(p4,b6m)*izb(p1,
c     &    p5)*izb(p5,b5m)*izb(p6,b6m)*zab2(p2,p3,p4,p1)*prop34**(-1)*
c     &    prop125**(-1) - za(p3,p6)*zb(p1,b5m)*zb(p4,b6m)*izb(p1,p5)*
c     &    izb(p5,b5m)*izb(p6,b6m)*zab2(p2,p3,p4,p1)*prop34**(-1)*
c     &    prop125**(-1)*Qsum - za(p3,p6)*zb(p1,b6m)*zb(p4,p6)*izb(p2,p5
c     &    )*izb(p5,b5m)*izb(p6,b6m)*zab2(p6,p2,p5,b5m)*prop34**(-1)*
c     &    prop125**(-1)
c      b(1,1) = b(1,1) + za(p3,p6)*zb(p1,b6m)*zb(p4,p6)*izb(p2,p5)*izb(
c     & p5,b5m)*izb(p6,b6m)*zab2(p6,p2,p5,b5m)*prop34**(-1)*
c     & prop125**(-1)*Qsum - 2*za(p3,p6)*zb(p4,b6m)*izb(p5,b5m)*izb(p6,
c     &    b6m)*zab2(p2,p3,p4,p1)*zab2(p5,p1,p2,b5m)*prop34**(-1)*
c     &    prop12**(-1)*prop125**(-1) + 1.D0/2.D0*zb(p1,p4)*zb(p1,b5m)*
c     &    izb(p1,p5)*izb(p2,p6)*izb(p5,b5m)*izb(p6,b6m)*zab2(p3,p2,p6,
c     &    b6m)*prop34**(-1) - 1.D0/2.D0*zb(p1,p4)*zb(p1,b5m)*izb(p1,p5)
c     &    *izb(p2,p6)*izb(p5,b5m)*izb(p6,b6m)*zab2(p3,p2,p6,b6m)*
c     &    prop34**(-1)*Qsum**2 + 1.D0/2.D0*zb(p1,p4)*zb(p1,b6m)*izb(p1,
c     &    p6)*izb(p2,p5)*izb(p5,b5m)*izb(p6,b6m)*zab2(p3,p2,p5,b5m)*
c     &    prop34**(-1) - 1.D0/2.D0*zb(p1,p4)*zb(p1,b6m)*izb(p1,p6)*izb(
c     &    p2,p5)*izb(p5,b5m)*izb(p6,b6m)*zab2(p3,p2,p5,b5m)*
c     &    prop34**(-1)*Qsum**2 - zb(p1,p4)*zb(p4,b5m)*izb(p2,p6)*izb(p4
c     &    ,p5)*izb(p5,b5m)*izb(p6,b6m)*zab2(p3,p2,p6,b6m)*prop126**(-1)
c     &     + zb(p1,p4)*zb(p4,b5m)*izb(p2,p6)*izb(p4,p5)*izb(p5,b5m)*
c     &    izb(p6,b6m)*zab2(p3,p2,p6,b6m)*prop126**(-1)*Qsum
c      b(1,1) = b(1,1) - zb(p1,p4)*zb(p4,b6m)*izb(p2,p5)*izb(p4,p6)*izb(
c     & p5,b5m)*izb(p6,b6m)*zab2(p3,p2,p5,b5m)*prop125**(-1) + zb(p1,p4)
c     &    *zb(p4,b6m)*izb(p2,p5)*izb(p4,p6)*izb(p5,b5m)*izb(p6,b6m)*
c     &    zab2(p3,p2,p5,b5m)*prop125**(-1)*Qsum - zb(p1,p4)*izb(p2,p5)*
c     &    izb(p5,b5m)*izb(p6,b6m)*zab2(p3,p2,p5,b5m)*zab2(p6,p3,p4,b6m)
c     &    *prop34**(-1)*prop125**(-1) + zb(p1,p4)*izb(p2,p5)*izb(p5,b5m
c     &    )*izb(p6,b6m)*zab2(p3,p2,p5,b5m)*zab2(p6,p3,p4,b6m)*
c     &    prop34**(-1)*prop125**(-1)*Qsum + 1.D0/2.D0*zb(p1,p4)*izb(p2,
c     &    p5)*izb(p5,b5m)*izb(p6,b6m)*zab2(p6,p2,p5,b5m)*zab3(p3,p2,p5,
c     &    p6,b6m)*prop34**(-1)*s256**(-1) - zb(p1,p4)*izb(p2,p5)*izb(p5
c     &    ,b5m)*izb(p6,b6m)*zab2(p6,p2,p5,b5m)*zab3(p3,p2,p5,p6,b6m)*
c     &    prop34**(-1)*Qsum*s256**(-1) + 1.D0/2.D0*zb(p1,p4)*izb(p2,p5)
c     &    *izb(p5,b5m)*izb(p6,b6m)*zab2(p6,p2,p5,b5m)*zab3(p3,p2,p5,p6,
c     &    b6m)*prop34**(-1)*Qsum**2*s256**(-1) - zb(p1,p4)*izb(p2,p6)*
c     &    izb(p5,b5m)*izb(p6,b6m)*zab2(p3,p2,p6,b6m)*zab2(p5,p3,p4,b5m)
c     &    *prop34**(-1)*prop126**(-1)
c      b(1,1) = b(1,1) + zb(p1,p4)*izb(p2,p6)*izb(p5,b5m)*izb(p6,b6m)*
c     & zab2(p3,p2,p6,b6m)*zab2(p5,p3,p4,b5m)*prop34**(-1)*prop126**(-1)
c     & *Qsum + 1.D0/2.D0*zb(p1,p4)*izb(p2,p6)*izb(p5,b5m)*izb(p6,b6m)*
c     &    zab2(p5,p2,p6,b6m)*zab3(p3,p2,p5,p6,b5m)*prop34**(-1)*
c     &    s256**(-1) - zb(p1,p4)*izb(p2,p6)*izb(p5,b5m)*izb(p6,b6m)*
c     &    zab2(p5,p2,p6,b6m)*zab3(p3,p2,p5,p6,b5m)*prop34**(-1)*Qsum*
c     &    s256**(-1) + 1.D0/2.D0*zb(p1,p4)*izb(p2,p6)*izb(p5,b5m)*izb(
c     &    p6,b6m)*zab2(p5,p2,p6,b6m)*zab3(p3,p2,p5,p6,b5m)*prop34**(-1)
c     &    *Qsum**2*s256**(-1)
c      write(6,*) 'b5m,b6m,b(1,1)',b5m,b6m,b(1,1)
c      enddo
c      enddo
c!
c      write(6,*)
c      do b5m=p1,p2,p2-p1
c      do b6p=p1,p2,p2-p1
c      b(1,2)= + 1.D0/2.D0*za(p2,p3)*za(p2,b6p)*zb(p1,p4)*zb(p1,b5m)*
c     &    iza(p2,p6)*iza(p6,b6p)*izb(p1,p5)*izb(p5,b5m)*prop34**(-1) -
c     &    1.D0/2.D0*za(p2,p3)*za(p2,b6p)*zb(p1,p4)*zb(p1,b5m)*iza(p2,p6
c     &    )*iza(p6,b6p)*izb(p1,p5)*izb(p5,b5m)*prop34**(-1)*Qsum**2 -
c     &    za(p2,p3)*za(p2,b6p)*zb(p1,p4)*zb(p4,b5m)*iza(p2,p6)*iza(p6,
c     &    b6p)*izb(p4,p5)*izb(p5,b5m)*prop126**(-1) + za(p2,p3)*za(p2,
c     &    b6p)*zb(p1,p4)*zb(p4,b5m)*iza(p2,p6)*iza(p6,b6p)*izb(p4,p5)*
c     &    izb(p5,b5m)*prop126**(-1)*Qsum - za(p2,p3)*za(p2,b6p)*zb(p1,
c     &    p4)*iza(p2,p6)*iza(p6,b6p)*izb(p5,b5m)*zab2(p5,p3,p4,b5m)*
c     &    prop34**(-1)*prop126**(-1) + za(p2,p3)*za(p2,b6p)*zb(p1,p4)*
c     &    iza(p2,p6)*iza(p6,b6p)*izb(p5,b5m)*zab2(p5,p3,p4,b5m)*
c     &    prop34**(-1)*prop126**(-1)*Qsum - 2*za(p2,p3)*za(p5,b6p)*zb(
c     &    p1,p4)*zb(p6,b5m)*iza(p6,b6p)*izb(p5,b5m)*prop34**(-1)*
c     &    prop12**(-1) - za(p2,p3)*zb(p1,p4)*zb(p1,b5m)*iza(p6,b6p)*
c     &    izb(p1,p5)*izb(p5,b5m)*zab2(b6p,p3,p4,p6)*prop34**(-1)*
c     &    prop125**(-1)
c      b(1,2) = b(1,2) - za(p2,p3)*zb(p1,p4)*zb(p1,b5m)*iza(p6,b6p)*izb(
c     & p1,p5)*izb(p5,b5m)*zab2(b6p,p3,p4,p6)*prop34**(-1)*prop125**(-1)
c     & *Qsum - 2*za(p2,p3)*zb(p1,p4)*zb(p4,b5m)*iza(p6,b6p)*izb(p4,p5)*
c     &    izb(p5,b5m)*zab2(b6p,p1,p2,p6)*prop12**(-1)*prop126**(-1) - 2
c     &    *za(p2,p3)*zb(p1,p4)*iza(p6,b6p)*izb(p5,b5m)*zab2(p5,p1,p2,
c     &    b5m)*zab2(b6p,p3,p4,p6)*prop34**(-1)*prop12**(-1)*
c     &    prop125**(-1) - 2*za(p2,p3)*zb(p1,p4)*iza(p6,b6p)*izb(p5,b5m)
c     &    *zab2(p5,p3,p4,b5m)*zab2(b6p,p1,p2,p6)*prop34**(-1)*
c     &    prop12**(-1)*prop126**(-1) + 1.D0/2.D0*za(p2,p3)*zb(p1,p6)*
c     &    zb(p1,b5m)*iza(p6,b6p)*izb(p1,p5)*izb(p5,b5m)*zab3(b6p,p1,p5,
c     &    p6,p4)*prop34**(-1)*s156**(-1) + za(p2,p3)*zb(p1,p6)*zb(p1,
c     &    b5m)*iza(p6,b6p)*izb(p1,p5)*izb(p5,b5m)*zab3(b6p,p1,p5,p6,p4)
c     &    *prop34**(-1)*Qsum*s156**(-1) + 1.D0/2.D0*za(p2,p3)*zb(p1,p6)
c     &    *zb(p1,b5m)*iza(p6,b6p)*izb(p1,p5)*izb(p5,b5m)*zab3(b6p,p1,p5
c     &    ,p6,p4)*prop34**(-1)*Qsum**2*s156**(-1) - za(p2,p3)*zb(p1,b5m
c     &    )*iza(p4,p6)*iza(p6,b6p)*izb(p1,p5)*izb(p5,b5m)*zab2(b6p,p4,
c     &    p6,p1)*prop125**(-1)
c      b(1,2) = b(1,2) - za(p2,p3)*zb(p1,b5m)*iza(p4,p6)*iza(p6,b6p)*
c     & izb(p1,p5)*izb(p5,b5m)*zab2(b6p,p4,p6,p1)*prop125**(-1)*Qsum -
c     &    za(p2,p3)*zb(p2,b5m)*iza(p4,p6)*iza(p6,b6p)*izb(p2,p5)*izb(p5
c     &    ,b5m)*zab2(b6p,p4,p6,p1)*prop125**(-1) + za(p2,p3)*zb(p2,b5m)
c     &    *iza(p4,p6)*iza(p6,b6p)*izb(p2,p5)*izb(p5,b5m)*zab2(b6p,p4,p6
c     &    ,p1)*prop125**(-1)*Qsum + 2*za(p2,p3)*zb(p4,p6)*zb(p4,b5m)*
c     &    iza(p6,b6p)*izb(p4,p5)*izb(p5,b5m)*zab2(b6p,p2,p3,p1)*
c     &    prop12**(-1)*s456**(-1) + za(p2,p3)*zb(p4,b5m)*iza(p1,p6)*
c     &    iza(p6,b6p)*izb(p4,p5)*izb(p5,b5m)*zab2(b6p,p1,p6,p4)*
c     &    prop126**(-1) + za(p2,p3)*zb(p4,b5m)*iza(p1,p6)*iza(p6,b6p)*
c     &    izb(p4,p5)*izb(p5,b5m)*zab2(b6p,p1,p6,p4)*prop126**(-1)*Qsum
c     &     - 1.D0/2.D0*za(p2,p3)*iza(p1,p6)*iza(p6,b6p)*izb(p5,b5m)*
c     &    zab2(p5,p1,p6,p4)*zab2(b6p,p1,p6,b5m)*prop34**(-1)*s156**(-1)
c     &     - za(p2,p3)*iza(p1,p6)*iza(p6,b6p)*izb(p5,b5m)*zab2(p5,p1,p6
c     &    ,p4)*zab2(b6p,p1,p6,b5m)*prop34**(-1)*Qsum*s156**(-1) - 1.D0/
c     &    2.D0*za(p2,p3)*iza(p1,p6)*iza(p6,b6p)*izb(p5,b5m)*zab2(p5,p1,
c     &    p6,p4)*zab2(b6p,p1,p6,b5m)*prop34**(-1)*Qsum**2*s156**(-1)
c      b(1,2) = b(1,2) + za(p2,p3)*iza(p1,p6)*iza(p6,b6p)*izb(p5,b5m)*
c     & zab2(p5,p3,p4,b5m)*zab2(b6p,p1,p6,p4)*prop34**(-1)*prop126**(-1)
c     &     + za(p2,p3)*iza(p1,p6)*iza(p6,b6p)*izb(p5,b5m)*zab2(p5,p3,p4
c     &    ,b5m)*zab2(b6p,p1,p6,p4)*prop34**(-1)*prop126**(-1)*Qsum - 2*
c     &    za(p2,p3)*iza(p4,p6)*iza(p6,b6p)*izb(p5,b5m)*zab2(p5,p1,p2,
c     &    b5m)*zab2(b6p,p4,p6,p1)*prop12**(-1)*prop125**(-1) + 2*za(p2,
c     &    p3)*iza(p4,p6)*iza(p6,b6p)*izb(p5,b5m)*zab2(p5,p4,p6,p1)*
c     &    zab2(b6p,p4,p6,b5m)*prop12**(-1)*s456**(-1) + za(p2,p5)*za(p2
c     &    ,b6p)*za(p3,p5)*zb(p1,p5)*zb(p4,b5m)*iza(p2,p6)*iza(p6,b6p)*
c     &    izb(p5,b5m)*prop34**(-1)*prop126**(-1) - za(p2,p5)*za(p2,b6p)
c     &    *za(p3,p5)*zb(p1,p5)*zb(p4,b5m)*iza(p2,p6)*iza(p6,b6p)*izb(p5
c     &    ,b5m)*prop34**(-1)*prop126**(-1)*Qsum - za(p2,p5)*za(p2,b6p)*
c     &    za(p3,p5)*zb(p1,b5m)*zb(p4,p5)*iza(p2,p6)*iza(p6,b6p)*izb(p5,
c     &    b5m)*prop34**(-1)*prop126**(-1) + za(p2,p5)*za(p2,b6p)*za(p3,
c     &    p5)*zb(p1,b5m)*zb(p4,p5)*iza(p2,p6)*iza(p6,b6p)*izb(p5,b5m)*
c     &    prop34**(-1)*prop126**(-1)*Qsum
c      b(1,2) = b(1,2) + 1.D0/2.D0*za(p2,p5)*za(p2,b6p)*zb(p1,p4)*iza(p2
c     & ,p6)*iza(p6,b6p)*izb(p5,b5m)*zab3(p3,p2,p5,p6,b5m)*prop34**(-1)*
c     & s256**(-1) - za(p2,p5)*za(p2,b6p)*zb(p1,p4)*iza(p2,p6)*iza(p6,
c     &    b6p)*izb(p5,b5m)*zab3(p3,p2,p5,p6,b5m)*prop34**(-1)*Qsum*
c     &    s256**(-1) + 1.D0/2.D0*za(p2,p5)*za(p2,b6p)*zb(p1,p4)*iza(p2,
c     &    p6)*iza(p6,b6p)*izb(p5,b5m)*zab3(p3,p2,p5,p6,b5m)*
c     &    prop34**(-1)*Qsum**2*s256**(-1) - 2*za(p2,p5)*za(p3,p5)*zb(p1
c     &    ,p5)*zb(p4,b5m)*iza(p6,b6p)*izb(p5,b5m)*zab2(b6p,p3,p4,p6)*
c     &    prop34**(-1)*prop12**(-1)*prop125**(-1) + 2*za(p2,p5)*za(p3,
c     &    p5)*zb(p1,p5)*iza(p4,p6)*iza(p6,b6p)*izb(p5,b5m)*zab2(b6p,p4,
c     &    p6,b5m)*prop12**(-1)*prop125**(-1) - 2*za(p2,p5)*za(p3,p5)*
c     &    zb(p1,b5m)*zb(p4,p5)*iza(p6,b6p)*izb(p5,b5m)*zab2(b6p,p1,p2,
c     &    p6)*prop34**(-1)*prop12**(-1)*prop126**(-1) - 2*za(p2,p5)*za(
c     &    p3,p5)*zb(p1,b5m)*iza(p4,p6)*iza(p6,b6p)*izb(p5,b5m)*zab2(b6p
c     &    ,p4,p6,p5)*prop12**(-1)*prop125**(-1) + za(p2,p5)*za(p3,p5)*
c     &    zb(p4,p5)*iza(p1,p6)*iza(p6,b6p)*izb(p5,b5m)*zab2(b6p,p1,p6,
c     &    b5m)*prop34**(-1)*prop126**(-1)
c      b(1,2) = b(1,2) + za(p2,p5)*za(p3,p5)*zb(p4,p5)*iza(p1,p6)*iza(p6
c     & ,b6p)*izb(p5,b5m)*zab2(b6p,p1,p6,b5m)*prop34**(-1)*prop126**(-1)
c     & *Qsum - za(p2,p5)*za(p3,p5)*zb(p4,b5m)*iza(p1,p6)*iza(p6,b6p)*
c     &    izb(p5,b5m)*zab2(b6p,p1,p6,p5)*prop34**(-1)*prop126**(-1) -
c     &    za(p2,p5)*za(p3,p5)*zb(p4,b5m)*iza(p1,p6)*iza(p6,b6p)*izb(p5,
c     &    b5m)*zab2(b6p,p1,p6,p5)*prop34**(-1)*prop126**(-1)*Qsum + 2*
c     &    za(p2,p5)*za(p3,p6)*za(p5,b6p)*zb(p1,p5)*zb(p4,p6)*zb(p6,b5m)
c     &    *iza(p6,b6p)*izb(p5,b5m)*prop34**(-1)*prop12**(-1)*
c     &    prop125**(-1) + 2*za(p2,p5)*za(p3,p6)*zb(p1,b5m)*zb(p4,p6)*
c     &    iza(p6,b6p)*izb(p5,b5m)*zab2(b6p,p1,p2,p6)*prop34**(-1)*
c     &    prop12**(-1)*prop125**(-1) + 2*za(p2,p5)*za(p3,b6p)*zb(p1,p5)
c     &    *zb(p4,p6)*iza(p6,b6p)*izb(p5,b5m)*zab2(p5,p3,p4,b5m)*
c     &    prop34**(-1)*prop12**(-1)*prop125**(-1) + za(p2,p5)*za(p3,b6p
c     &    )*zb(p1,b5m)*zb(p4,p6)*iza(p6,b6p)*izb(p5,b5m)*prop34**(-1)*
c     &    prop12**(-1)*prop125**(-1)*Mwsq**(-1)*s34*s12 + 2*za(p2,p5)*
c     &    za(p3,b6p)*zb(p1,b5m)*zb(p4,p6)*iza(p6,b6p)*izb(p5,b5m)*
c     &    prop34**(-1)*prop12**(-1)*prop125**(-1)*s56
c      b(1,2) = b(1,2) - za(p2,p5)*za(p3,b6p)*zb(p1,b5m)*zb(p4,p6)*iza(
c     & p6,b6p)*izb(p5,b5m)*prop34**(-1)*prop12**(-1)*prop125**(-1)*s12
c     &     - za(p2,p5)*za(p3,b6p)*zb(p1,b5m)*zb(p4,p6)*iza(p6,b6p)*izb(
c     &    p5,b5m)*prop34**(-1)*prop12**(-1)*prop125**(-1)*s34 + za(p2,
c     &    p5)*za(p3,b6p)*zb(p1,b5m)*zb(p4,p6)*iza(p6,b6p)*izb(p5,b5m)*
c     &    prop34**(-1)*prop12**(-1)*prop125**(-1)*s125 - za(p2,p5)*za(
c     &    p3,b6p)*zb(p1,b5m)*zb(p4,p6)*iza(p6,b6p)*izb(p5,b5m)*
c     &    prop34**(-1)*prop12**(-1) - za(p2,p5)*za(p3,b6p)*zb(p1,b5m)*
c     &    zb(p4,p6)*iza(p6,b6p)*izb(p5,b5m)*prop34**(-1)*prop125**(-1)*
c     &    Mwsq**(-1)*s34 - za(p2,p5)*za(p3,b6p)*zb(p1,b5m)*zb(p4,p6)*
c     &    iza(p6,b6p)*izb(p5,b5m)*prop34**(-1)*prop125**(-1)*Qsum - za(
c     &    p2,p5)*za(p3,b6p)*zb(p1,b5m)*zb(p4,p6)*iza(p6,b6p)*izb(p5,b5m
c     &    )*prop12**(-1)*prop125**(-1)*Mwsq**(-1)*s12 + za(p2,p5)*za(p3
c     &    ,b6p)*zb(p1,b5m)*zb(p4,p6)*iza(p6,b6p)*izb(p5,b5m)*
c     &    prop12**(-1)*prop125**(-1) + za(p2,p5)*za(p3,b6p)*zb(p1,b5m)*
c     &    zb(p4,p6)*iza(p6,b6p)*izb(p5,b5m)*prop125**(-1)*Mwsq**(-1)
c      b(1,2) = b(1,2) + 2*za(p2,p5)*zb(p1,b5m)*iza(p6,b6p)*izb(p5,b5m)*
c     & zab2(p3,p1,p2,p4)*zab2(b6p,p3,p4,p6)*prop34**(-1)*prop12**(-1)*
c     & prop125**(-1) + 2*za(p2,p6)*za(p3,p5)*za(p5,b6p)*zb(p1,p6)*zb(p4
c     &    ,p5)*zb(p6,b5m)*iza(p6,b6p)*izb(p5,b5m)*prop34**(-1)*
c     &    prop12**(-1)*prop126**(-1) + 2*za(p2,p6)*za(p3,p5)*zb(p1,p6)*
c     &    zb(p4,b5m)*iza(p6,b6p)*izb(p5,b5m)*zab2(b6p,p3,p4,p6)*
c     &    prop34**(-1)*prop12**(-1)*prop126**(-1) - 2*za(p2,p6)*za(p3,
c     &    b6p)*zb(p1,p6)*zb(p4,p6)*zb(p4,b5m)*iza(p6,b6p)*izb(p4,p5)*
c     &    izb(p5,b5m)*prop12**(-1)*prop126**(-1) - 2*za(p2,p6)*za(p3,
c     &    b6p)*zb(p1,p6)*zb(p4,p6)*iza(p6,b6p)*izb(p5,b5m)*zab2(p5,p3,
c     &    p4,b5m)*prop34**(-1)*prop12**(-1)*prop126**(-1) + 2*za(p2,b6p
c     &    )*za(p3,p5)*zb(p1,p6)*zb(p4,p5)*iza(p6,b6p)*izb(p5,b5m)*zab2(
c     &    p5,p1,p2,b5m)*prop34**(-1)*prop12**(-1)*prop126**(-1) + za(p2
c     &    ,b6p)*za(p3,p5)*zb(p1,p6)*zb(p4,b5m)*iza(p6,b6p)*izb(p5,b5m)*
c     &    prop34**(-1)*prop12**(-1)*prop126**(-1)*Mwsq**(-1)*s34*s12 +
c     &    2*za(p2,b6p)*za(p3,p5)*zb(p1,p6)*zb(p4,b5m)*iza(p6,b6p)*izb(
c     &    p5,b5m)*prop34**(-1)*prop12**(-1)*prop126**(-1)*s56
c      b(1,2) = b(1,2) - za(p2,b6p)*za(p3,p5)*zb(p1,p6)*zb(p4,b5m)*iza(
c     & p6,b6p)*izb(p5,b5m)*prop34**(-1)*prop12**(-1)*prop126**(-1)*s12
c     &     - za(p2,b6p)*za(p3,p5)*zb(p1,p6)*zb(p4,b5m)*iza(p6,b6p)*izb(
c     &    p5,b5m)*prop34**(-1)*prop12**(-1)*prop126**(-1)*s34 + za(p2,
c     &    b6p)*za(p3,p5)*zb(p1,p6)*zb(p4,b5m)*iza(p6,b6p)*izb(p5,b5m)*
c     &    prop34**(-1)*prop12**(-1)*prop126**(-1)*s126 - za(p2,b6p)*za(
c     &    p3,p5)*zb(p1,p6)*zb(p4,b5m)*iza(p6,b6p)*izb(p5,b5m)*
c     &    prop34**(-1)*prop12**(-1) - za(p2,b6p)*za(p3,p5)*zb(p1,p6)*
c     &    zb(p4,b5m)*iza(p6,b6p)*izb(p5,b5m)*prop34**(-1)*prop126**(-1)
c     &    *Mwsq**(-1)*s34 + za(p2,b6p)*za(p3,p5)*zb(p1,p6)*zb(p4,b5m)*
c     &    iza(p6,b6p)*izb(p5,b5m)*prop34**(-1)*prop126**(-1) - za(p2,
c     &    b6p)*za(p3,p5)*zb(p1,p6)*zb(p4,b5m)*iza(p6,b6p)*izb(p5,b5m)*
c     &    prop12**(-1)*prop126**(-1)*Mwsq**(-1)*s12 - za(p2,b6p)*za(p3,
c     &    p5)*zb(p1,p6)*zb(p4,b5m)*iza(p6,b6p)*izb(p5,b5m)*prop12**(-1)
c     &    *prop126**(-1) + za(p2,b6p)*za(p3,p5)*zb(p1,p6)*zb(p4,b5m)*
c     &    iza(p6,b6p)*izb(p5,b5m)*prop126**(-1)*Mwsq**(-1)
c      b(1,2) = b(1,2) - za(p2,b6p)*za(p3,p6)*zb(p1,p6)*zb(p1,b5m)*zb(p4
c     & ,p6)*iza(p6,b6p)*izb(p1,p5)*izb(p5,b5m)*prop34**(-1)*
c     & prop125**(-1) - za(p2,b6p)*za(p3,p6)*zb(p1,p6)*zb(p1,b5m)*zb(p4,
c     &    p6)*iza(p6,b6p)*izb(p1,p5)*izb(p5,b5m)*prop34**(-1)*
c     &    prop125**(-1)*Qsum - 2*za(p2,b6p)*za(p3,p6)*zb(p1,p6)*zb(p4,
c     &    p6)*iza(p6,b6p)*izb(p5,b5m)*zab2(p5,p1,p2,b5m)*prop34**(-1)*
c     &    prop12**(-1)*prop125**(-1) + 2*za(p2,b6p)*zb(p1,p6)*zb(p4,b5m
c     &    )*iza(p6,b6p)*izb(p4,p5)*izb(p5,b5m)*zab2(p3,p1,p2,p4)*
c     &    prop12**(-1)*prop126**(-1) + 2*za(p2,b6p)*zb(p1,p6)*iza(p6,
c     &    b6p)*izb(p5,b5m)*zab2(p3,p1,p2,p4)*zab2(p5,p3,p4,b5m)*
c     &    prop34**(-1)*prop12**(-1)*prop126**(-1) + 2*za(p3,p5)*zb(p4,
c     &    b5m)*iza(p6,b6p)*izb(p5,b5m)*zab2(p2,p3,p4,p1)*zab2(b6p,p1,p2
c     &    ,p6)*prop34**(-1)*prop12**(-1)*prop126**(-1) + za(p3,p5)*iza(
c     &    p4,p6)*iza(p6,b6p)*izb(p2,p5)*zab2(b6p,p4,p6,p1)*
c     &    prop125**(-1) - za(p3,p5)*iza(p4,p6)*iza(p6,b6p)*izb(p2,p5)*
c     &    zab2(b6p,p4,p6,p1)*prop125**(-1)*Qsum
c      b(1,2) = b(1,2) + za(p3,p6)*zb(p1,p6)*zb(p4,p6)*iza(p6,b6p)*izb(
c     & p2,p5)*izb(p5,b5m)*zab2(b6p,p2,p5,b5m)*prop34**(-1)*
c     & prop125**(-1) - za(p3,p6)*zb(p1,p6)*zb(p4,p6)*iza(p6,b6p)*izb(p2
c     &    ,p5)*izb(p5,b5m)*zab2(b6p,p2,p5,b5m)*prop34**(-1)*
c     &    prop125**(-1)*Qsum - za(p3,b6p)*zb(p1,p6)*zb(p4,p6)*iza(p6,
c     &    b6p)*izb(p2,p5)*izb(p5,b5m)*zab2(p6,p2,p5,b5m)*prop34**(-1)*
c     &    prop125**(-1) + za(p3,b6p)*zb(p1,p6)*zb(p4,p6)*iza(p6,b6p)*
c     &    izb(p2,p5)*izb(p5,b5m)*zab2(p6,p2,p5,b5m)*prop34**(-1)*
c     &    prop125**(-1)*Qsum + za(p3,b6p)*zb(p1,b5m)*zb(p4,p6)*iza(p6,
c     &    b6p)*izb(p1,p5)*izb(p5,b5m)*zab2(p2,p3,p4,p1)*prop34**(-1)*
c     &    prop125**(-1) + za(p3,b6p)*zb(p1,b5m)*zb(p4,p6)*iza(p6,b6p)*
c     &    izb(p1,p5)*izb(p5,b5m)*zab2(p2,p3,p4,p1)*prop34**(-1)*
c     &    prop125**(-1)*Qsum + 2*za(p3,b6p)*zb(p4,p6)*iza(p6,b6p)*izb(
c     &    p5,b5m)*zab2(p2,p3,p4,p1)*zab2(p5,p1,p2,b5m)*prop34**(-1)*
c     &    prop12**(-1)*prop125**(-1) - 1.D0/2.D0*zb(p1,p4)*iza(p6,b6p)*
c     &    izb(p2,p5)*izb(p5,b5m)*zab2(p3,p2,p5,p6)*zab2(b6p,p2,p5,b5m)*
c     &    prop34**(-1)*s256**(-1)
c      b(1,2) = b(1,2) + zb(p1,p4)*iza(p6,b6p)*izb(p2,p5)*izb(p5,b5m)*
c     & zab2(p3,p2,p5,p6)*zab2(b6p,p2,p5,b5m)*prop34**(-1)*Qsum*
c     & s256**(-1) - 1.D0/2.D0*zb(p1,p4)*iza(p6,b6p)*izb(p2,p5)*izb(p5,
c     &    b5m)*zab2(p3,p2,p5,p6)*zab2(b6p,p2,p5,b5m)*prop34**(-1)*
c     &    Qsum**2*s256**(-1) + zb(p1,p4)*iza(p6,b6p)*izb(p2,p5)*izb(p5,
c     &    b5m)*zab2(p3,p2,p5,b5m)*zab2(b6p,p3,p4,p6)*prop34**(-1)*
c     &    prop125**(-1) - zb(p1,p4)*iza(p6,b6p)*izb(p2,p5)*izb(p5,b5m)*
c     &    zab2(p3,p2,p5,b5m)*zab2(b6p,p3,p4,p6)*prop34**(-1)*
c     &    prop125**(-1)*Qsum + 1.D0/2.D0*iza(p1,p6)*iza(p6,b6p)*izb(p2,
c     &    p5)*izb(p5,b5m)*zab2(p3,p2,p5,b5m)*zab2(b6p,p1,p6,p4)*
c     &    prop34**(-1) - 1.D0/2.D0*iza(p1,p6)*iza(p6,b6p)*izb(p2,p5)*
c     &    izb(p5,b5m)*zab2(p3,p2,p5,b5m)*zab2(b6p,p1,p6,p4)*
c     &    prop34**(-1)*Qsum**2
c      write(6,*) 'b5m,b6p,b(1,2)',b5m,b6p,b(1,2)
c      enddo
c      enddo
c!
c      write(6,*)
c      do b5p=p1,p2,p2-p1
c      do b6m=p1,p2,p2-p1
c      b(2,1)= + 1.D0/2.D0*za(p2,p3)*za(p2,b5p)*zb(p1,p4)*zb(p1,b6m)*
c     &    iza(p2,p5)*iza(p5,b5p)*izb(p1,p6)*izb(p6,b6m)*prop34**(-1) -
c     &    1.D0/2.D0*za(p2,p3)*za(p2,b5p)*zb(p1,p4)*zb(p1,b6m)*iza(p2,p5
c     &    )*iza(p5,b5p)*izb(p1,p6)*izb(p6,b6m)*prop34**(-1)*Qsum**2 -
c     &    za(p2,p3)*za(p2,b5p)*zb(p1,p4)*zb(p4,b6m)*iza(p2,p5)*iza(p5,
c     &    b5p)*izb(p4,p6)*izb(p6,b6m)*prop125**(-1) + za(p2,p3)*za(p2,
c     &    b5p)*zb(p1,p4)*zb(p4,b6m)*iza(p2,p5)*iza(p5,b5p)*izb(p4,p6)*
c     &    izb(p6,b6m)*prop125**(-1)*Qsum - za(p2,p3)*za(p2,b5p)*zb(p1,
c     &    p4)*iza(p2,p5)*iza(p5,b5p)*izb(p6,b6m)*zab2(p6,p3,p4,b6m)*
c     &    prop34**(-1)*prop125**(-1) + za(p2,p3)*za(p2,b5p)*zb(p1,p4)*
c     &    iza(p2,p5)*iza(p5,b5p)*izb(p6,b6m)*zab2(p6,p3,p4,b6m)*
c     &    prop34**(-1)*prop125**(-1)*Qsum - 2*za(p2,p3)*za(p6,b5p)*zb(
c     &    p1,p4)*zb(p5,b6m)*iza(p5,b5p)*izb(p6,b6m)*prop34**(-1)*
c     &    prop12**(-1) - za(p2,p3)*zb(p1,p4)*zb(p1,b6m)*iza(p5,b5p)*
c     &    izb(p1,p6)*izb(p6,b6m)*zab2(b5p,p3,p4,p5)*prop34**(-1)*
c     &    prop126**(-1)
c      b(2,1) = b(2,1) - za(p2,p3)*zb(p1,p4)*zb(p1,b6m)*iza(p5,b5p)*izb(
c     & p1,p6)*izb(p6,b6m)*zab2(b5p,p3,p4,p5)*prop34**(-1)*prop126**(-1)
c     & *Qsum - 2*za(p2,p3)*zb(p1,p4)*zb(p4,b6m)*iza(p5,b5p)*izb(p4,p6)*
c     &    izb(p6,b6m)*zab2(b5p,p1,p2,p5)*prop12**(-1)*prop125**(-1) - 2
c     &    *za(p2,p3)*zb(p1,p4)*iza(p5,b5p)*izb(p6,b6m)*zab2(p6,p1,p2,
c     &    b6m)*zab2(b5p,p3,p4,p5)*prop34**(-1)*prop12**(-1)*
c     &    prop126**(-1) - 2*za(p2,p3)*zb(p1,p4)*iza(p5,b5p)*izb(p6,b6m)
c     &    *zab2(p6,p3,p4,b6m)*zab2(b5p,p1,p2,p5)*prop34**(-1)*
c     &    prop12**(-1)*prop125**(-1) + 1.D0/2.D0*za(p2,p3)*zb(p1,p5)*
c     &    zb(p1,b6m)*iza(p5,b5p)*izb(p1,p6)*izb(p6,b6m)*zab3(b5p,p1,p5,
c     &    p6,p4)*prop34**(-1)*s156**(-1) + za(p2,p3)*zb(p1,p5)*zb(p1,
c     &    b6m)*iza(p5,b5p)*izb(p1,p6)*izb(p6,b6m)*zab3(b5p,p1,p5,p6,p4)
c     &    *prop34**(-1)*Qsum*s156**(-1) + 1.D0/2.D0*za(p2,p3)*zb(p1,p5)
c     &    *zb(p1,b6m)*iza(p5,b5p)*izb(p1,p6)*izb(p6,b6m)*zab3(b5p,p1,p5
c     &    ,p6,p4)*prop34**(-1)*Qsum**2*s156**(-1) - za(p2,p3)*zb(p1,b6m
c     &    )*iza(p4,p5)*iza(p5,b5p)*izb(p1,p6)*izb(p6,b6m)*zab2(b5p,p4,
c     &    p5,p1)*prop126**(-1)
c      b(2,1) = b(2,1) - za(p2,p3)*zb(p1,b6m)*iza(p4,p5)*iza(p5,b5p)*
c     & izb(p1,p6)*izb(p6,b6m)*zab2(b5p,p4,p5,p1)*prop126**(-1)*Qsum -
c     &    za(p2,p3)*zb(p2,b6m)*iza(p4,p5)*iza(p5,b5p)*izb(p2,p6)*izb(p6
c     &    ,b6m)*zab2(b5p,p4,p5,p1)*prop126**(-1) + za(p2,p3)*zb(p2,b6m)
c     &    *iza(p4,p5)*iza(p5,b5p)*izb(p2,p6)*izb(p6,b6m)*zab2(b5p,p4,p5
c     &    ,p1)*prop126**(-1)*Qsum + 2*za(p2,p3)*zb(p4,p5)*zb(p4,b6m)*
c     &    iza(p5,b5p)*izb(p4,p6)*izb(p6,b6m)*zab2(b5p,p2,p3,p1)*
c     &    prop12**(-1)*s456**(-1) + za(p2,p3)*zb(p4,b6m)*iza(p1,p5)*
c     &    iza(p5,b5p)*izb(p4,p6)*izb(p6,b6m)*zab2(b5p,p1,p5,p4)*
c     &    prop125**(-1) + za(p2,p3)*zb(p4,b6m)*iza(p1,p5)*iza(p5,b5p)*
c     &    izb(p4,p6)*izb(p6,b6m)*zab2(b5p,p1,p5,p4)*prop125**(-1)*Qsum
c     &     - 1.D0/2.D0*za(p2,p3)*iza(p1,p5)*iza(p5,b5p)*izb(p6,b6m)*
c     &    zab2(p6,p1,p5,p4)*zab2(b5p,p1,p5,b6m)*prop34**(-1)*s156**(-1)
c     &     - za(p2,p3)*iza(p1,p5)*iza(p5,b5p)*izb(p6,b6m)*zab2(p6,p1,p5
c     &    ,p4)*zab2(b5p,p1,p5,b6m)*prop34**(-1)*Qsum*s156**(-1) - 1.D0/
c     &    2.D0*za(p2,p3)*iza(p1,p5)*iza(p5,b5p)*izb(p6,b6m)*zab2(p6,p1,
c     &    p5,p4)*zab2(b5p,p1,p5,b6m)*prop34**(-1)*Qsum**2*s156**(-1)
c      b(2,1) = b(2,1) + za(p2,p3)*iza(p1,p5)*iza(p5,b5p)*izb(p6,b6m)*
c     & zab2(p6,p3,p4,b6m)*zab2(b5p,p1,p5,p4)*prop34**(-1)*prop125**(-1)
c     &     + za(p2,p3)*iza(p1,p5)*iza(p5,b5p)*izb(p6,b6m)*zab2(p6,p3,p4
c     &    ,b6m)*zab2(b5p,p1,p5,p4)*prop34**(-1)*prop125**(-1)*Qsum - 2*
c     &    za(p2,p3)*iza(p4,p5)*iza(p5,b5p)*izb(p6,b6m)*zab2(p6,p1,p2,
c     &    b6m)*zab2(b5p,p4,p5,p1)*prop12**(-1)*prop126**(-1) + 2*za(p2,
c     &    p3)*iza(p4,p5)*iza(p5,b5p)*izb(p6,b6m)*zab2(p6,p4,p5,p1)*
c     &    zab2(b5p,p4,p5,b6m)*prop12**(-1)*s456**(-1) + 2*za(p2,p5)*za(
c     &    p3,p6)*za(p6,b5p)*zb(p1,p5)*zb(p4,p6)*zb(p5,b6m)*iza(p5,b5p)*
c     &    izb(p6,b6m)*prop34**(-1)*prop12**(-1)*prop125**(-1) + 2*za(p2
c     &    ,p5)*za(p3,p6)*zb(p1,p5)*zb(p4,b6m)*iza(p5,b5p)*izb(p6,b6m)*
c     &    zab2(b5p,p3,p4,p5)*prop34**(-1)*prop12**(-1)*prop125**(-1) -
c     &    2*za(p2,p5)*za(p3,b5p)*zb(p1,p5)*zb(p4,p5)*zb(p4,b6m)*iza(p5,
c     &    b5p)*izb(p4,p6)*izb(p6,b6m)*prop12**(-1)*prop125**(-1) - 2*
c     &    za(p2,p5)*za(p3,b5p)*zb(p1,p5)*zb(p4,p5)*iza(p5,b5p)*izb(p6,
c     &    b6m)*zab2(p6,p3,p4,b6m)*prop34**(-1)*prop12**(-1)*
c     &    prop125**(-1)
c      b(2,1) = b(2,1) + za(p2,p6)*za(p2,b5p)*za(p3,p6)*zb(p1,p6)*zb(p4,
c     & b6m)*iza(p2,p5)*iza(p5,b5p)*izb(p6,b6m)*prop34**(-1)*
c     & prop125**(-1) - za(p2,p6)*za(p2,b5p)*za(p3,p6)*zb(p1,p6)*zb(p4,
c     &    b6m)*iza(p2,p5)*iza(p5,b5p)*izb(p6,b6m)*prop34**(-1)*
c     &    prop125**(-1)*Qsum - za(p2,p6)*za(p2,b5p)*za(p3,p6)*zb(p1,b6m
c     &    )*zb(p4,p6)*iza(p2,p5)*iza(p5,b5p)*izb(p6,b6m)*prop34**(-1)*
c     &    prop125**(-1) + za(p2,p6)*za(p2,b5p)*za(p3,p6)*zb(p1,b6m)*zb(
c     &    p4,p6)*iza(p2,p5)*iza(p5,b5p)*izb(p6,b6m)*prop34**(-1)*
c     &    prop125**(-1)*Qsum + 1.D0/2.D0*za(p2,p6)*za(p2,b5p)*zb(p1,p4)
c     &    *iza(p2,p5)*iza(p5,b5p)*izb(p6,b6m)*zab3(p3,p2,p5,p6,b6m)*
c     &    prop34**(-1)*s256**(-1) - za(p2,p6)*za(p2,b5p)*zb(p1,p4)*iza(
c     &    p2,p5)*iza(p5,b5p)*izb(p6,b6m)*zab3(p3,p2,p5,p6,b6m)*
c     &    prop34**(-1)*Qsum*s256**(-1) + 1.D0/2.D0*za(p2,p6)*za(p2,b5p)
c     &    *zb(p1,p4)*iza(p2,p5)*iza(p5,b5p)*izb(p6,b6m)*zab3(p3,p2,p5,
c     &    p6,b6m)*prop34**(-1)*Qsum**2*s256**(-1) + 2*za(p2,p6)*za(p3,
c     &    p5)*za(p6,b5p)*zb(p1,p6)*zb(p4,p5)*zb(p5,b6m)*iza(p5,b5p)*
c     &    izb(p6,b6m)*prop34**(-1)*prop12**(-1)*prop126**(-1)
c      b(2,1) = b(2,1) + 2*za(p2,p6)*za(p3,p5)*zb(p1,b6m)*zb(p4,p5)*iza(
c     & p5,b5p)*izb(p6,b6m)*zab2(b5p,p1,p2,p5)*prop34**(-1)*prop12**(-1)
c     & *prop126**(-1) - 2*za(p2,p6)*za(p3,p6)*zb(p1,p6)*zb(p4,b6m)*iza(
c     &    p5,b5p)*izb(p6,b6m)*zab2(b5p,p3,p4,p5)*prop34**(-1)*
c     &    prop12**(-1)*prop126**(-1) + 2*za(p2,p6)*za(p3,p6)*zb(p1,p6)*
c     &    iza(p4,p5)*iza(p5,b5p)*izb(p6,b6m)*zab2(b5p,p4,p5,b6m)*
c     &    prop12**(-1)*prop126**(-1) - 2*za(p2,p6)*za(p3,p6)*zb(p1,b6m)
c     &    *zb(p4,p6)*iza(p5,b5p)*izb(p6,b6m)*zab2(b5p,p1,p2,p5)*
c     &    prop34**(-1)*prop12**(-1)*prop125**(-1) - 2*za(p2,p6)*za(p3,
c     &    p6)*zb(p1,b6m)*iza(p4,p5)*iza(p5,b5p)*izb(p6,b6m)*zab2(b5p,p4
c     &    ,p5,p6)*prop12**(-1)*prop126**(-1) + za(p2,p6)*za(p3,p6)*zb(
c     &    p4,p6)*iza(p1,p5)*iza(p5,b5p)*izb(p6,b6m)*zab2(b5p,p1,p5,b6m)
c     &    *prop34**(-1)*prop125**(-1) + za(p2,p6)*za(p3,p6)*zb(p4,p6)*
c     &    iza(p1,p5)*iza(p5,b5p)*izb(p6,b6m)*zab2(b5p,p1,p5,b6m)*
c     &    prop34**(-1)*prop125**(-1)*Qsum - za(p2,p6)*za(p3,p6)*zb(p4,
c     &    b6m)*iza(p1,p5)*iza(p5,b5p)*izb(p6,b6m)*zab2(b5p,p1,p5,p6)*
c     &    prop34**(-1)*prop125**(-1)
c      b(2,1) = b(2,1) - za(p2,p6)*za(p3,p6)*zb(p4,b6m)*iza(p1,p5)*iza(
c     & p5,b5p)*izb(p6,b6m)*zab2(b5p,p1,p5,p6)*prop34**(-1)*
c     & prop125**(-1)*Qsum + 2*za(p2,p6)*za(p3,b5p)*zb(p1,p6)*zb(p4,p5)*
c     &    iza(p5,b5p)*izb(p6,b6m)*zab2(p6,p3,p4,b6m)*prop34**(-1)*
c     &    prop12**(-1)*prop126**(-1) + za(p2,p6)*za(p3,b5p)*zb(p1,b6m)*
c     &    zb(p4,p5)*iza(p5,b5p)*izb(p6,b6m)*prop34**(-1)*prop12**(-1)*
c     &    prop126**(-1)*Mwsq**(-1)*s34*s12 + 2*za(p2,p6)*za(p3,b5p)*zb(
c     &    p1,b6m)*zb(p4,p5)*iza(p5,b5p)*izb(p6,b6m)*prop34**(-1)*
c     &    prop12**(-1)*prop126**(-1)*s56 - za(p2,p6)*za(p3,b5p)*zb(p1,
c     &    b6m)*zb(p4,p5)*iza(p5,b5p)*izb(p6,b6m)*prop34**(-1)*
c     &    prop12**(-1)*prop126**(-1)*s12 - za(p2,p6)*za(p3,b5p)*zb(p1,
c     &    b6m)*zb(p4,p5)*iza(p5,b5p)*izb(p6,b6m)*prop34**(-1)*
c     &    prop12**(-1)*prop126**(-1)*s34 + za(p2,p6)*za(p3,b5p)*zb(p1,
c     &    b6m)*zb(p4,p5)*iza(p5,b5p)*izb(p6,b6m)*prop34**(-1)*
c     &    prop12**(-1)*prop126**(-1)*s126 - za(p2,p6)*za(p3,b5p)*zb(p1,
c     &    b6m)*zb(p4,p5)*iza(p5,b5p)*izb(p6,b6m)*prop34**(-1)*
c     &    prop12**(-1)
c      b(2,1) = b(2,1) - za(p2,p6)*za(p3,b5p)*zb(p1,b6m)*zb(p4,p5)*iza(
c     & p5,b5p)*izb(p6,b6m)*prop34**(-1)*prop126**(-1)*Mwsq**(-1)*s34 -
c     &    za(p2,p6)*za(p3,b5p)*zb(p1,b6m)*zb(p4,p5)*iza(p5,b5p)*izb(p6,
c     &    b6m)*prop34**(-1)*prop126**(-1)*Qsum - za(p2,p6)*za(p3,b5p)*
c     &    zb(p1,b6m)*zb(p4,p5)*iza(p5,b5p)*izb(p6,b6m)*prop12**(-1)*
c     &    prop126**(-1)*Mwsq**(-1)*s12 + za(p2,p6)*za(p3,b5p)*zb(p1,b6m
c     &    )*zb(p4,p5)*iza(p5,b5p)*izb(p6,b6m)*prop12**(-1)*
c     &    prop126**(-1) + za(p2,p6)*za(p3,b5p)*zb(p1,b6m)*zb(p4,p5)*
c     &    iza(p5,b5p)*izb(p6,b6m)*prop126**(-1)*Mwsq**(-1) + 2*za(p2,p6
c     &    )*zb(p1,b6m)*iza(p5,b5p)*izb(p6,b6m)*zab2(p3,p1,p2,p4)*zab2(
c     &    b5p,p3,p4,p5)*prop34**(-1)*prop12**(-1)*prop126**(-1) - za(p2
c     &    ,b5p)*za(p3,p5)*zb(p1,p5)*zb(p1,b6m)*zb(p4,p5)*iza(p5,b5p)*
c     &    izb(p1,p6)*izb(p6,b6m)*prop34**(-1)*prop126**(-1) - za(p2,b5p
c     &    )*za(p3,p5)*zb(p1,p5)*zb(p1,b6m)*zb(p4,p5)*iza(p5,b5p)*izb(p1
c     &    ,p6)*izb(p6,b6m)*prop34**(-1)*prop126**(-1)*Qsum - 2*za(p2,
c     &    b5p)*za(p3,p5)*zb(p1,p5)*zb(p4,p5)*iza(p5,b5p)*izb(p6,b6m)*
c     &    zab2(p6,p1,p2,b6m)*prop34**(-1)*prop12**(-1)*prop126**(-1)
c      b(2,1) = b(2,1) + 2*za(p2,b5p)*za(p3,p6)*zb(p1,p5)*zb(p4,p6)*iza(
c     & p5,b5p)*izb(p6,b6m)*zab2(p6,p1,p2,b6m)*prop34**(-1)*prop12**(-1)
c     & *prop125**(-1) + za(p2,b5p)*za(p3,p6)*zb(p1,p5)*zb(p4,b6m)*iza(
c     &    p5,b5p)*izb(p6,b6m)*prop34**(-1)*prop12**(-1)*prop125**(-1)*
c     &    Mwsq**(-1)*s34*s12 + 2*za(p2,b5p)*za(p3,p6)*zb(p1,p5)*zb(p4,
c     &    b6m)*iza(p5,b5p)*izb(p6,b6m)*prop34**(-1)*prop12**(-1)*
c     &    prop125**(-1)*s56 - za(p2,b5p)*za(p3,p6)*zb(p1,p5)*zb(p4,b6m)
c     &    *iza(p5,b5p)*izb(p6,b6m)*prop34**(-1)*prop12**(-1)*
c     &    prop125**(-1)*s12 - za(p2,b5p)*za(p3,p6)*zb(p1,p5)*zb(p4,b6m)
c     &    *iza(p5,b5p)*izb(p6,b6m)*prop34**(-1)*prop12**(-1)*
c     &    prop125**(-1)*s34 + za(p2,b5p)*za(p3,p6)*zb(p1,p5)*zb(p4,b6m)
c     &    *iza(p5,b5p)*izb(p6,b6m)*prop34**(-1)*prop12**(-1)*
c     &    prop125**(-1)*s125 - za(p2,b5p)*za(p3,p6)*zb(p1,p5)*zb(p4,b6m
c     &    )*iza(p5,b5p)*izb(p6,b6m)*prop34**(-1)*prop12**(-1) - za(p2,
c     &    b5p)*za(p3,p6)*zb(p1,p5)*zb(p4,b6m)*iza(p5,b5p)*izb(p6,b6m)*
c     &    prop34**(-1)*prop125**(-1)*Mwsq**(-1)*s34
c      b(2,1) = b(2,1) + za(p2,b5p)*za(p3,p6)*zb(p1,p5)*zb(p4,b6m)*iza(
c     & p5,b5p)*izb(p6,b6m)*prop34**(-1)*prop125**(-1) - za(p2,b5p)*za(
c     &    p3,p6)*zb(p1,p5)*zb(p4,b6m)*iza(p5,b5p)*izb(p6,b6m)*
c     &    prop12**(-1)*prop125**(-1)*Mwsq**(-1)*s12 - za(p2,b5p)*za(p3,
c     &    p6)*zb(p1,p5)*zb(p4,b6m)*iza(p5,b5p)*izb(p6,b6m)*prop12**(-1)
c     &    *prop125**(-1) + za(p2,b5p)*za(p3,p6)*zb(p1,p5)*zb(p4,b6m)*
c     &    iza(p5,b5p)*izb(p6,b6m)*prop125**(-1)*Mwsq**(-1) + 2*za(p2,
c     &    b5p)*zb(p1,p5)*zb(p4,b6m)*iza(p5,b5p)*izb(p4,p6)*izb(p6,b6m)*
c     &    zab2(p3,p1,p2,p4)*prop12**(-1)*prop125**(-1) + 2*za(p2,b5p)*
c     &    zb(p1,p5)*iza(p5,b5p)*izb(p6,b6m)*zab2(p3,p1,p2,p4)*zab2(p6,
c     &    p3,p4,b6m)*prop34**(-1)*prop12**(-1)*prop125**(-1) + za(p3,p5
c     &    )*zb(p1,p5)*zb(p4,p5)*iza(p5,b5p)*izb(p2,p6)*izb(p6,b6m)*
c     &    zab2(b5p,p2,p6,b6m)*prop34**(-1)*prop126**(-1) - za(p3,p5)*
c     &    zb(p1,p5)*zb(p4,p5)*iza(p5,b5p)*izb(p2,p6)*izb(p6,b6m)*zab2(
c     &    b5p,p2,p6,b6m)*prop34**(-1)*prop126**(-1)*Qsum + 2*za(p3,p6)*
c     &    zb(p4,b6m)*iza(p5,b5p)*izb(p6,b6m)*zab2(p2,p3,p4,p1)*zab2(b5p
c     &    ,p1,p2,p5)*prop34**(-1)*prop12**(-1)*prop125**(-1)
c      b(2,1) = b(2,1) + za(p3,p6)*iza(p4,p5)*iza(p5,b5p)*izb(p2,p6)*
c     & zab2(b5p,p4,p5,p1)*prop126**(-1) - za(p3,p6)*iza(p4,p5)*iza(p5,
c     &    b5p)*izb(p2,p6)*zab2(b5p,p4,p5,p1)*prop126**(-1)*Qsum - za(p3
c     &    ,b5p)*zb(p1,p5)*zb(p4,p5)*iza(p5,b5p)*izb(p2,p6)*izb(p6,b6m)*
c     &    zab2(p5,p2,p6,b6m)*prop34**(-1)*prop126**(-1) + za(p3,b5p)*
c     &    zb(p1,p5)*zb(p4,p5)*iza(p5,b5p)*izb(p2,p6)*izb(p6,b6m)*zab2(
c     &    p5,p2,p6,b6m)*prop34**(-1)*prop126**(-1)*Qsum + za(p3,b5p)*
c     &    zb(p1,b6m)*zb(p4,p5)*iza(p5,b5p)*izb(p1,p6)*izb(p6,b6m)*zab2(
c     &    p2,p3,p4,p1)*prop34**(-1)*prop126**(-1) + za(p3,b5p)*zb(p1,
c     &    b6m)*zb(p4,p5)*iza(p5,b5p)*izb(p1,p6)*izb(p6,b6m)*zab2(p2,p3,
c     &    p4,p1)*prop34**(-1)*prop126**(-1)*Qsum + 2*za(p3,b5p)*zb(p4,
c     &    p5)*iza(p5,b5p)*izb(p6,b6m)*zab2(p2,p3,p4,p1)*zab2(p6,p1,p2,
c     &    b6m)*prop34**(-1)*prop12**(-1)*prop126**(-1) - 1.D0/2.D0*zb(
c     &    p1,p4)*iza(p5,b5p)*izb(p2,p6)*izb(p6,b6m)*zab2(p3,p2,p6,p5)*
c     &    zab2(b5p,p2,p6,b6m)*prop34**(-1)*s256**(-1) + zb(p1,p4)*iza(
c     &    p5,b5p)*izb(p2,p6)*izb(p6,b6m)*zab2(p3,p2,p6,p5)*zab2(b5p,p2,
c     &    p6,b6m)*prop34**(-1)*Qsum*s256**(-1)
c      b(2,1) = b(2,1) - 1.D0/2.D0*zb(p1,p4)*iza(p5,b5p)*izb(p2,p6)*izb(
c     & p6,b6m)*zab2(p3,p2,p6,p5)*zab2(b5p,p2,p6,b6m)*prop34**(-1)*
c     & Qsum**2*s256**(-1) + zb(p1,p4)*iza(p5,b5p)*izb(p2,p6)*izb(p6,b6m
c     &    )*zab2(p3,p2,p6,b6m)*zab2(b5p,p3,p4,p5)*prop34**(-1)*
c     &    prop126**(-1) - zb(p1,p4)*iza(p5,b5p)*izb(p2,p6)*izb(p6,b6m)*
c     &    zab2(p3,p2,p6,b6m)*zab2(b5p,p3,p4,p5)*prop34**(-1)*
c     &    prop126**(-1)*Qsum + 1.D0/2.D0*iza(p1,p5)*iza(p5,b5p)*izb(p2,
c     &    p6)*izb(p6,b6m)*zab2(p3,p2,p6,b6m)*zab2(b5p,p1,p5,p4)*
c     &    prop34**(-1) - 1.D0/2.D0*iza(p1,p5)*iza(p5,b5p)*izb(p2,p6)*
c     &    izb(p6,b6m)*zab2(p3,p2,p6,b6m)*zab2(b5p,p1,p5,p4)*
c     &    prop34**(-1)*Qsum**2
c      write(6,*) 'b5p,b6m,b(2,1)',b5p,b6m,b(2,1)
c      enddo
c      enddo
c!
c      write(6,*)
c      do b5p=p1,p2,p2-p1
c      do b6p=p1,p2,p2-p1
c      b(2,2)= + za(p1,b5p)*za(p2,p3)*iza(p1,p5)*iza(p4,p6)*iza(p5,b5p)*
c     &    iza(p6,b6p)*zab2(b6p,p4,p6,p1)*prop125**(-1) + za(p1,b5p)*za(
c     &    p2,p3)*iza(p1,p5)*iza(p4,p6)*iza(p5,b5p)*iza(p6,b6p)*zab2(b6p
c     &    ,p4,p6,p1)*prop125**(-1)*Qsum + za(p1,b6p)*za(p2,p3)*iza(p1,
c     &    p6)*iza(p4,p5)*iza(p5,b5p)*iza(p6,b6p)*zab2(b5p,p4,p5,p1)*
c     &    prop126**(-1) + za(p1,b6p)*za(p2,p3)*iza(p1,p6)*iza(p4,p5)*
c     &    iza(p5,b5p)*iza(p6,b6p)*zab2(b5p,p4,p5,p1)*prop126**(-1)*Qsum
c     &     + za(p2,p3)*za(p2,b5p)*zb(p1,p4)*iza(p2,p5)*iza(p5,b5p)*iza(
c     &    p6,b6p)*zab2(b6p,p3,p4,p6)*prop34**(-1)*prop125**(-1) - za(p2
c     &    ,p3)*za(p2,b5p)*zb(p1,p4)*iza(p2,p5)*iza(p5,b5p)*iza(p6,b6p)*
c     &    zab2(b6p,p3,p4,p6)*prop34**(-1)*prop125**(-1)*Qsum + 1.D0/2.D0
c     &    *za(p2,p3)*za(p2,b5p)*iza(p1,p6)*iza(p2,p5)*iza(p5,b5p)*iza(
c     &    p6,b6p)*zab2(b6p,p1,p6,p4)*prop34**(-1) - 1.D0/2.D0*za(p2,p3)
c     &    *za(p2,b5p)*iza(p1,p6)*iza(p2,p5)*iza(p5,b5p)*iza(p6,b6p)*
c     &    zab2(b6p,p1,p6,p4)*prop34**(-1)*Qsum**2 + za(p2,p3)*za(p2,b5p
c     &    )*iza(p2,p5)*iza(p4,p6)*iza(p5,b5p)*iza(p6,b6p)*zab2(b6p,p4,
c     &    p6,p1)*prop125**(-1)
c      b(2,2) = b(2,2) - za(p2,p3)*za(p2,b5p)*iza(p2,p5)*iza(p4,p6)*iza(
c     & p5,b5p)*iza(p6,b6p)*zab2(b6p,p4,p6,p1)*prop125**(-1)*Qsum + za(
c     &    p2,p3)*za(p2,b6p)*zb(p1,p4)*iza(p2,p6)*iza(p5,b5p)*iza(p6,b6p
c     &    )*zab2(b5p,p3,p4,p5)*prop34**(-1)*prop126**(-1) - za(p2,p3)*
c     &    za(p2,b6p)*zb(p1,p4)*iza(p2,p6)*iza(p5,b5p)*iza(p6,b6p)*zab2(
c     &    b5p,p3,p4,p5)*prop34**(-1)*prop126**(-1)*Qsum + 1.D0/2.D0*za(
c     &    p2,p3)*za(p2,b6p)*iza(p1,p5)*iza(p2,p6)*iza(p5,b5p)*iza(p6,
c     &    b6p)*zab2(b5p,p1,p5,p4)*prop34**(-1) - 1.D0/2.D0*za(p2,p3)*
c     &    za(p2,b6p)*iza(p1,p5)*iza(p2,p6)*iza(p5,b5p)*iza(p6,b6p)*
c     &    zab2(b5p,p1,p5,p4)*prop34**(-1)*Qsum**2 + za(p2,p3)*za(p2,b6p
c     &    )*iza(p2,p6)*iza(p4,p5)*iza(p5,b5p)*iza(p6,b6p)*zab2(b5p,p4,
c     &    p5,p1)*prop126**(-1) - za(p2,p3)*za(p2,b6p)*iza(p2,p6)*iza(p4
c     &    ,p5)*iza(p5,b5p)*iza(p6,b6p)*zab2(b5p,p4,p5,p1)*prop126**(-1)
c     &    *Qsum - 2*za(p2,p3)*za(b5p,b6p)*zb(p1,p4)*zb(p5,p6)*iza(p5,
c     &    b5p)*iza(p6,b6p)*prop34**(-1)*prop12**(-1) + 2*za(p2,p3)*zb(
c     &    p1,p4)*iza(p5,b5p)*iza(p6,b6p)*zab2(b5p,p1,p2,p5)*zab2(b6p,p3
c     &    ,p4,p6)*prop34**(-1)*prop12**(-1)*prop125**(-1)
c      b(2,2) = b(2,2) + 2*za(p2,p3)*zb(p1,p4)*iza(p5,b5p)*iza(p6,b6p)*
c     & zab2(b5p,p3,p4,p5)*zab2(b6p,p1,p2,p6)*prop34**(-1)*prop12**(-1)*
c     & prop126**(-1) + za(p2,p3)*iza(p1,p5)*iza(p4,p6)*iza(p6,b6p)*
c     &    zab2(b6p,p4,p6,p5)*prop125**(-1) + za(p2,p3)*iza(p1,p5)*iza(
c     &    p4,p6)*iza(p6,b6p)*zab2(b6p,p4,p6,p5)*prop125**(-1)*Qsum -
c     &    za(p2,p3)*iza(p1,p5)*iza(p5,b5p)*iza(p6,b6p)*zab2(b5p,p1,p5,
c     &    p4)*zab2(b6p,p3,p4,p6)*prop34**(-1)*prop125**(-1) - za(p2,p3)
c     &    *iza(p1,p5)*iza(p5,b5p)*iza(p6,b6p)*zab2(b5p,p1,p5,p4)*zab2(
c     &    b6p,p3,p4,p6)*prop34**(-1)*prop125**(-1)*Qsum + 1.D0/2.D0*za(
c     &    p2,p3)*iza(p1,p5)*iza(p5,b5p)*iza(p6,b6p)*zab2(b5p,p1,p5,p6)*
c     &    zab3(b6p,p1,p5,p6,p4)*prop34**(-1)*s156**(-1) + za(p2,p3)*
c     &    iza(p1,p5)*iza(p5,b5p)*iza(p6,b6p)*zab2(b5p,p1,p5,p6)*zab3(
c     &    b6p,p1,p5,p6,p4)*prop34**(-1)*Qsum*s156**(-1) + 1.D0/2.D0*za(
c     &    p2,p3)*iza(p1,p5)*iza(p5,b5p)*iza(p6,b6p)*zab2(b5p,p1,p5,p6)*
c     &    zab3(b6p,p1,p5,p6,p4)*prop34**(-1)*Qsum**2*s156**(-1) + za(p2
c     &    ,p3)*iza(p1,p6)*iza(p4,p5)*iza(p5,b5p)*zab2(b5p,p4,p5,p6)*
c     &    prop126**(-1)
c      b(2,2) = b(2,2) + za(p2,p3)*iza(p1,p6)*iza(p4,p5)*iza(p5,b5p)*
c     & zab2(b5p,p4,p5,p6)*prop126**(-1)*Qsum - za(p2,p3)*iza(p1,p6)*
c     &    iza(p5,b5p)*iza(p6,b6p)*zab2(b5p,p3,p4,p5)*zab2(b6p,p1,p6,p4)
c     &    *prop34**(-1)*prop126**(-1) - za(p2,p3)*iza(p1,p6)*iza(p5,b5p
c     &    )*iza(p6,b6p)*zab2(b5p,p3,p4,p5)*zab2(b6p,p1,p6,p4)*
c     &    prop34**(-1)*prop126**(-1)*Qsum + 1.D0/2.D0*za(p2,p3)*iza(p1,
c     &    p6)*iza(p5,b5p)*iza(p6,b6p)*zab2(b6p,p1,p6,p5)*zab3(b5p,p1,p5
c     &    ,p6,p4)*prop34**(-1)*s156**(-1) + za(p2,p3)*iza(p1,p6)*iza(p5
c     &    ,b5p)*iza(p6,b6p)*zab2(b6p,p1,p6,p5)*zab3(b5p,p1,p5,p6,p4)*
c     &    prop34**(-1)*Qsum*s156**(-1) + 1.D0/2.D0*za(p2,p3)*iza(p1,p6)
c     &    *iza(p5,b5p)*iza(p6,b6p)*zab2(b6p,p1,p6,p5)*zab3(b5p,p1,p5,p6
c     &    ,p4)*prop34**(-1)*Qsum**2*s156**(-1) + 2*za(p2,p3)*iza(p4,p5)
c     &    *iza(p5,b5p)*iza(p6,b6p)*zab2(b5p,p4,p5,p1)*zab2(b6p,p1,p2,p6
c     &    )*prop12**(-1)*prop126**(-1) + 2*za(p2,p3)*iza(p4,p5)*iza(p5,
c     &    b5p)*iza(p6,b6p)*zab2(b5p,p4,p5,p6)*zab2(b6p,p2,p3,p1)*
c     &    prop12**(-1)*s456**(-1)
c      b(2,2) = b(2,2) + 2*za(p2,p3)*iza(p4,p6)*iza(p5,b5p)*iza(p6,b6p)*
c     & zab2(b5p,p1,p2,p5)*zab2(b6p,p4,p6,p1)*prop12**(-1)*prop125**(-1)
c     &     + 2*za(p2,p3)*iza(p4,p6)*iza(p5,b5p)*iza(p6,b6p)*zab2(b5p,p2
c     &    ,p3,p1)*zab2(b6p,p4,p6,p5)*prop12**(-1)*s456**(-1) - za(p2,p5
c     &    )*za(p2,b6p)*za(p3,b5p)*zb(p1,p5)*zb(p4,p5)*iza(p2,p6)*iza(p5
c     &    ,b5p)*iza(p6,b6p)*prop34**(-1)*prop126**(-1) + za(p2,p5)*za(
c     &    p2,b6p)*za(p3,b5p)*zb(p1,p5)*zb(p4,p5)*iza(p2,p6)*iza(p5,b5p)
c     &    *iza(p6,b6p)*prop34**(-1)*prop126**(-1)*Qsum + 2*za(p2,p5)*
c     &    za(p3,p6)*za(b5p,b6p)*zb(p1,p5)*zb(p4,p6)*zb(p5,p6)*iza(p5,
c     &    b5p)*iza(p6,b6p)*prop34**(-1)*prop12**(-1)*prop125**(-1) + 2*
c     &    za(p2,p5)*za(p3,b5p)*zb(p1,p5)*zb(p4,p5)*iza(p5,b5p)*iza(p6,
c     &    b6p)*zab2(b6p,p3,p4,p6)*prop34**(-1)*prop12**(-1)*
c     &    prop125**(-1) - 2*za(p2,p5)*za(p3,b5p)*zb(p1,p5)*iza(p4,p6)*
c     &    iza(p5,b5p)*iza(p6,b6p)*zab2(b6p,p4,p6,p5)*prop12**(-1)*
c     &    prop125**(-1) + za(p2,p5)*za(p3,b5p)*zb(p4,p5)*iza(p1,p6)*
c     &    iza(p5,b5p)*iza(p6,b6p)*zab2(b6p,p1,p6,p5)*prop34**(-1)*
c     &    prop126**(-1)
c      b(2,2) = b(2,2) + za(p2,p5)*za(p3,b5p)*zb(p4,p5)*iza(p1,p6)*iza(
c     & p5,b5p)*iza(p6,b6p)*zab2(b6p,p1,p6,p5)*prop34**(-1)*
c     & prop126**(-1)*Qsum - 2*za(p2,p5)*za(p3,b6p)*zb(p1,p5)*zb(p4,p6)*
c     &    iza(p5,b5p)*iza(p6,b6p)*zab2(b5p,p3,p4,p5)*prop34**(-1)*
c     &    prop12**(-1)*prop125**(-1) - za(p2,p6)*za(p2,b5p)*za(p3,b6p)*
c     &    zb(p1,p6)*zb(p4,p6)*iza(p2,p5)*iza(p5,b5p)*iza(p6,b6p)*
c     &    prop34**(-1)*prop125**(-1) + za(p2,p6)*za(p2,b5p)*za(p3,b6p)*
c     &    zb(p1,p6)*zb(p4,p6)*iza(p2,p5)*iza(p5,b5p)*iza(p6,b6p)*
c     &    prop34**(-1)*prop125**(-1)*Qsum + 2*za(p2,p6)*za(p3,p5)*za(
c     &    b5p,b6p)*zb(p1,p6)*zb(p4,p5)*zb(p5,p6)*iza(p5,b5p)*iza(p6,b6p
c     &    )*prop34**(-1)*prop12**(-1)*prop126**(-1) - 2*za(p2,p6)*za(p3
c     &    ,b5p)*zb(p1,p6)*zb(p4,p5)*iza(p5,b5p)*iza(p6,b6p)*zab2(b6p,p3
c     &    ,p4,p6)*prop34**(-1)*prop12**(-1)*prop126**(-1) + 2*za(p2,p6)
c     &    *za(p3,b6p)*zb(p1,p6)*zb(p4,p6)*iza(p5,b5p)*iza(p6,b6p)*zab2(
c     &    b5p,p3,p4,p5)*prop34**(-1)*prop12**(-1)*prop126**(-1) - 2*za(
c     &    p2,p6)*za(p3,b6p)*zb(p1,p6)*iza(p4,p5)*iza(p5,b5p)*iza(p6,b6p
c     &    )*zab2(b5p,p4,p5,p6)*prop12**(-1)*prop126**(-1)
c      b(2,2) = b(2,2) + za(p2,p6)*za(p3,b6p)*zb(p4,p6)*iza(p1,p5)*iza(
c     & p5,b5p)*iza(p6,b6p)*zab2(b5p,p1,p5,p6)*prop34**(-1)*
c     & prop125**(-1) + za(p2,p6)*za(p3,b6p)*zb(p4,p6)*iza(p1,p5)*iza(p5
c     &    ,b5p)*iza(p6,b6p)*zab2(b5p,p1,p5,p6)*prop34**(-1)*
c     &    prop125**(-1)*Qsum + za(p2,b5p)*za(p2,b6p)*za(p3,p5)*zb(p1,p5
c     &    )*zb(p4,p5)*iza(p2,p6)*iza(p5,b5p)*iza(p6,b6p)*prop34**(-1)*
c     &    prop126**(-1) - za(p2,b5p)*za(p2,b6p)*za(p3,p5)*zb(p1,p5)*zb(
c     &    p4,p5)*iza(p2,p6)*iza(p5,b5p)*iza(p6,b6p)*prop34**(-1)*
c     &    prop126**(-1)*Qsum + za(p2,b5p)*za(p2,b6p)*za(p3,p6)*zb(p1,p6
c     &    )*zb(p4,p6)*iza(p2,p5)*iza(p5,b5p)*iza(p6,b6p)*prop34**(-1)*
c     &    prop125**(-1) - za(p2,b5p)*za(p2,b6p)*za(p3,p6)*zb(p1,p6)*zb(
c     &    p4,p6)*iza(p2,p5)*iza(p5,b5p)*iza(p6,b6p)*prop34**(-1)*
c     &    prop125**(-1)*Qsum - 1.D0/2.D0*za(p2,b5p)*za(p2,b6p)*zb(p1,p4
c     &    )*iza(p2,p5)*iza(p5,b5p)*iza(p6,b6p)*zab2(p3,p2,p5,p6)*
c     &    prop34**(-1)*s256**(-1) + za(p2,b5p)*za(p2,b6p)*zb(p1,p4)*
c     &    iza(p2,p5)*iza(p5,b5p)*iza(p6,b6p)*zab2(p3,p2,p5,p6)*
c     &    prop34**(-1)*Qsum*s256**(-1)
c      b(2,2) = b(2,2) - 1.D0/2.D0*za(p2,b5p)*za(p2,b6p)*zb(p1,p4)*iza(
c     & p2,p5)*iza(p5,b5p)*iza(p6,b6p)*zab2(p3,p2,p5,p6)*prop34**(-1)*
c     & Qsum**2*s256**(-1) - 1.D0/2.D0*za(p2,b5p)*za(p2,b6p)*zb(p1,p4)*
c     &    iza(p2,p6)*iza(p5,b5p)*iza(p6,b6p)*zab2(p3,p2,p6,p5)*
c     &    prop34**(-1)*s256**(-1) + za(p2,b5p)*za(p2,b6p)*zb(p1,p4)*
c     &    iza(p2,p6)*iza(p5,b5p)*iza(p6,b6p)*zab2(p3,p2,p6,p5)*
c     &    prop34**(-1)*Qsum*s256**(-1) - 1.D0/2.D0*za(p2,b5p)*za(p2,b6p
c     &    )*zb(p1,p4)*iza(p2,p6)*iza(p5,b5p)*iza(p6,b6p)*zab2(p3,p2,p6,
c     &    p5)*prop34**(-1)*Qsum**2*s256**(-1) + 2*za(p2,b5p)*za(p3,p5)*
c     &    zb(p1,p5)*zb(p4,p5)*iza(p5,b5p)*iza(p6,b6p)*zab2(b6p,p1,p2,p6
c     &    )*prop34**(-1)*prop12**(-1)*prop126**(-1) + 2*za(p2,b5p)*za(
c     &    p3,p5)*zb(p1,p5)*iza(p4,p6)*iza(p5,b5p)*iza(p6,b6p)*zab2(b6p,
c     &    p4,p6,p5)*prop12**(-1)*prop125**(-1) - za(p2,b5p)*za(p3,p5)*
c     &    zb(p4,p5)*iza(p1,p6)*iza(p5,b5p)*iza(p6,b6p)*zab2(b6p,p1,p6,
c     &    p5)*prop34**(-1)*prop126**(-1) - za(p2,b5p)*za(p3,p5)*zb(p4,
c     &    p5)*iza(p1,p6)*iza(p5,b5p)*iza(p6,b6p)*zab2(b6p,p1,p6,p5)*
c     &    prop34**(-1)*prop126**(-1)*Qsum
c      b(2,2) = b(2,2) - 2*za(p2,b5p)*za(p3,p6)*zb(p1,p5)*zb(p4,p6)*iza(
c     & p5,b5p)*iza(p6,b6p)*zab2(b6p,p1,p2,p6)*prop34**(-1)*prop12**(-1)
c     & *prop125**(-1) - za(p2,b5p)*za(p3,b6p)*zb(p1,p5)*zb(p4,p6)*iza(
c     &    p5,b5p)*iza(p6,b6p)*prop34**(-1)*prop12**(-1)*prop125**(-1)*
c     &    Mwsq**(-1)*s34*s12 - 2*za(p2,b5p)*za(p3,b6p)*zb(p1,p5)*zb(p4,
c     &    p6)*iza(p5,b5p)*iza(p6,b6p)*prop34**(-1)*prop12**(-1)*
c     &    prop125**(-1)*s56 + za(p2,b5p)*za(p3,b6p)*zb(p1,p5)*zb(p4,p6)
c     &    *iza(p5,b5p)*iza(p6,b6p)*prop34**(-1)*prop12**(-1)*
c     &    prop125**(-1)*s12 + za(p2,b5p)*za(p3,b6p)*zb(p1,p5)*zb(p4,p6)
c     &    *iza(p5,b5p)*iza(p6,b6p)*prop34**(-1)*prop12**(-1)*
c     &    prop125**(-1)*s34 - za(p2,b5p)*za(p3,b6p)*zb(p1,p5)*zb(p4,p6)
c     &    *iza(p5,b5p)*iza(p6,b6p)*prop34**(-1)*prop12**(-1)*
c     &    prop125**(-1)*s125 + za(p2,b5p)*za(p3,b6p)*zb(p1,p5)*zb(p4,p6
c     &    )*iza(p5,b5p)*iza(p6,b6p)*prop34**(-1)*prop12**(-1) + za(p2,
c     &    b5p)*za(p3,b6p)*zb(p1,p5)*zb(p4,p6)*iza(p5,b5p)*iza(p6,b6p)*
c     &    prop34**(-1)*prop125**(-1)*Mwsq**(-1)*s34
c      b(2,2) = b(2,2) - za(p2,b5p)*za(p3,b6p)*zb(p1,p5)*zb(p4,p6)*iza(
c     & p5,b5p)*iza(p6,b6p)*prop34**(-1)*prop125**(-1) + za(p2,b5p)*za(
c     &    p3,b6p)*zb(p1,p5)*zb(p4,p6)*iza(p5,b5p)*iza(p6,b6p)*
c     &    prop12**(-1)*prop125**(-1)*Mwsq**(-1)*s12 - za(p2,b5p)*za(p3,
c     &    b6p)*zb(p1,p5)*zb(p4,p6)*iza(p5,b5p)*iza(p6,b6p)*prop12**(-1)
c     &    *prop125**(-1) - za(p2,b5p)*za(p3,b6p)*zb(p1,p5)*zb(p4,p6)*
c     &    iza(p5,b5p)*iza(p6,b6p)*prop125**(-1)*Mwsq**(-1) - 2*za(p2,
c     &    b5p)*zb(p1,p5)*iza(p5,b5p)*iza(p6,b6p)*zab2(p3,p1,p2,p4)*
c     &    zab2(b6p,p3,p4,p6)*prop34**(-1)*prop12**(-1)*prop125**(-1) -
c     &    2*za(p2,b6p)*za(p3,p5)*zb(p1,p6)*zb(p4,p5)*iza(p5,b5p)*iza(p6
c     &    ,b6p)*zab2(b5p,p1,p2,p5)*prop34**(-1)*prop12**(-1)*
c     &    prop126**(-1) + 2*za(p2,b6p)*za(p3,p6)*zb(p1,p6)*zb(p4,p6)*
c     &    iza(p5,b5p)*iza(p6,b6p)*zab2(b5p,p1,p2,p5)*prop34**(-1)*
c     &    prop12**(-1)*prop125**(-1) + 2*za(p2,b6p)*za(p3,p6)*zb(p1,p6)
c     &    *iza(p4,p5)*iza(p5,b5p)*iza(p6,b6p)*zab2(b5p,p4,p5,p6)*
c     &    prop12**(-1)*prop126**(-1)
c      b(2,2) = b(2,2) - za(p2,b6p)*za(p3,p6)*zb(p4,p6)*iza(p1,p5)*iza(
c     & p5,b5p)*iza(p6,b6p)*zab2(b5p,p1,p5,p6)*prop34**(-1)*
c     & prop125**(-1) - za(p2,b6p)*za(p3,p6)*zb(p4,p6)*iza(p1,p5)*iza(p5
c     &    ,b5p)*iza(p6,b6p)*zab2(b5p,p1,p5,p6)*prop34**(-1)*
c     &    prop125**(-1)*Qsum - za(p2,b6p)*za(p3,b5p)*zb(p1,p6)*zb(p4,p5
c     &    )*iza(p5,b5p)*iza(p6,b6p)*prop34**(-1)*prop12**(-1)*
c     &    prop126**(-1)*Mwsq**(-1)*s34*s12 - 2*za(p2,b6p)*za(p3,b5p)*
c     &    zb(p1,p6)*zb(p4,p5)*iza(p5,b5p)*iza(p6,b6p)*prop34**(-1)*
c     &    prop12**(-1)*prop126**(-1)*s56 + za(p2,b6p)*za(p3,b5p)*zb(p1,
c     &    p6)*zb(p4,p5)*iza(p5,b5p)*iza(p6,b6p)*prop34**(-1)*
c     &    prop12**(-1)*prop126**(-1)*s12 + za(p2,b6p)*za(p3,b5p)*zb(p1,
c     &    p6)*zb(p4,p5)*iza(p5,b5p)*iza(p6,b6p)*prop34**(-1)*
c     &    prop12**(-1)*prop126**(-1)*s34 - za(p2,b6p)*za(p3,b5p)*zb(p1,
c     &    p6)*zb(p4,p5)*iza(p5,b5p)*iza(p6,b6p)*prop34**(-1)*
c     &    prop12**(-1)*prop126**(-1)*s126 + za(p2,b6p)*za(p3,b5p)*zb(p1
c     &    ,p6)*zb(p4,p5)*iza(p5,b5p)*iza(p6,b6p)*prop34**(-1)*
c     &    prop12**(-1)
c      b(2,2) = b(2,2) + za(p2,b6p)*za(p3,b5p)*zb(p1,p6)*zb(p4,p5)*iza(
c     & p5,b5p)*iza(p6,b6p)*prop34**(-1)*prop126**(-1)*Mwsq**(-1)*s34 -
c     &    za(p2,b6p)*za(p3,b5p)*zb(p1,p6)*zb(p4,p5)*iza(p5,b5p)*iza(p6,
c     &    b6p)*prop34**(-1)*prop126**(-1) + za(p2,b6p)*za(p3,b5p)*zb(p1
c     &    ,p6)*zb(p4,p5)*iza(p5,b5p)*iza(p6,b6p)*prop12**(-1)*
c     &    prop126**(-1)*Mwsq**(-1)*s12 - za(p2,b6p)*za(p3,b5p)*zb(p1,p6
c     &    )*zb(p4,p5)*iza(p5,b5p)*iza(p6,b6p)*prop12**(-1)*
c     &    prop126**(-1) - za(p2,b6p)*za(p3,b5p)*zb(p1,p6)*zb(p4,p5)*
c     &    iza(p5,b5p)*iza(p6,b6p)*prop126**(-1)*Mwsq**(-1) - 2*za(p2,
c     &    b6p)*zb(p1,p6)*iza(p5,b5p)*iza(p6,b6p)*zab2(p3,p1,p2,p4)*
c     &    zab2(b5p,p3,p4,p5)*prop34**(-1)*prop12**(-1)*prop126**(-1) -
c     &    2*za(p3,b5p)*zb(p4,p5)*iza(p5,b5p)*iza(p6,b6p)*zab2(p2,p3,p4,
c     &    p1)*zab2(b6p,p1,p2,p6)*prop34**(-1)*prop12**(-1)*
c     &    prop126**(-1) - 2*za(p3,b6p)*zb(p4,p6)*iza(p5,b5p)*iza(p6,b6p
c     &    )*zab2(p2,p3,p4,p1)*zab2(b5p,p1,p2,p5)*prop34**(-1)*
c     &    prop12**(-1)*prop125**(-1)
c      write(6,*) 'b5p,b6p,b(2,2)',b5p,b6p,b(2,2)
c      enddo
c      enddo

c      write(6,*)

      b(1,1)= - 2*za(p2,p3)*za(p2,p5)*zb(p1,p2)*zb(p1,p4)**2*izb(p1,p5)
     &    *izb(p1,p6)*izb(p4,p6)*prop12**(-1)*prop125**(-1) + 2*za(p2,
     &    p3)*za(p2,p5)*zb(p1,p2)*zb(p1,p4)*izb(p1,p5)*izb(p1,p6)*zab2(
     &    p6,p3,p4,p1)*prop34**(-1)*prop12**(-1)*prop125**(-1) - 2*za(
     &    p2,p3)*za(p2,p6)*zb(p1,p2)*zb(p1,p4)**2*izb(p1,p5)*izb(p1,p6)
     &    *izb(p4,p5)*prop12**(-1)*prop126**(-1) + 2*za(p2,p3)*za(p2,p6
     &    )*zb(p1,p2)*zb(p1,p4)*izb(p1,p5)*izb(p1,p6)*zab2(p5,p3,p4,p1)
     &    *prop34**(-1)*prop12**(-1)*prop126**(-1) + 2*za(p2,p3)*zb(p1,
     &    p4)**2*izb(p1,p5)*izb(p1,p6)*izb(p4,p5)*zab2(p6,p4,p5,p1)*
     &    prop12**(-1)*s456**(-1) + 2*za(p2,p3)*zb(p1,p4)**2*izb(p1,p5)
     &    *izb(p1,p6)*izb(p4,p6)*zab2(p5,p4,p6,p1)*prop12**(-1)*
     &    s456**(-1) + 2*za(p2,p5)*za(p3,p5)*zb(p1,p4)**2*izb(p1,p6)*
     &    izb(p4,p6)*prop12**(-1)*prop125**(-1) - 2*za(p2,p5)*za(p3,p5)
     &    *zb(p1,p4)*izb(p1,p6)*zab2(p6,p3,p4,p1)*prop34**(-1)*
     &    prop12**(-1)*prop125**(-1) + 2*za(p2,p5)*za(p3,p6)*zb(p1,p2)*
     &    zb(p1,p4)*izb(p1,p5)*izb(p1,p6)*zab2(p2,p3,p4,p1)*
     &    prop34**(-1)*prop12**(-1)*prop125**(-1)
      b(1,1) = b(1,1) + 2*za(p2,p5)*za(p3,p6)*zb(p1,p4)*izb(p1,p6)*
     & zab2(p5,p3,p4,p1)*prop34**(-1)*prop12**(-1)*prop125**(-1) + 2*
     &    za(p2,p6)*za(p3,p5)*zb(p1,p2)*zb(p1,p4)*izb(p1,p5)*izb(p1,p6)
     &    *zab2(p2,p3,p4,p1)*prop34**(-1)*prop12**(-1)*prop126**(-1) +
     &    2*za(p2,p6)*za(p3,p5)*zb(p1,p4)*izb(p1,p5)*zab2(p6,p3,p4,p1)*
     &    prop34**(-1)*prop12**(-1)*prop126**(-1) + 2*za(p2,p6)*za(p3,
     &    p6)*zb(p1,p4)**2*izb(p1,p5)*izb(p4,p5)*prop12**(-1)*
     &    prop126**(-1) - 2*za(p2,p6)*za(p3,p6)*zb(p1,p4)*izb(p1,p5)*
     &    zab2(p5,p3,p4,p1)*prop34**(-1)*prop12**(-1)*prop126**(-1) + 1.
     &    D0/2.D0*za(p3,p4)*zb(p1,p4)**2*izb(p1,p5)*izb(p1,p6)*izb(p2,
     &    p5)*zab2(p6,p2,p5,p1)*prop34**(-1)*s256**(-1) - za(p3,p4)*zb(
     &    p1,p4)**2*izb(p1,p5)*izb(p1,p6)*izb(p2,p5)*zab2(p6,p2,p5,p1)*
     &    prop34**(-1)*Qsum*s256**(-1) + 1.D0/2.D0*za(p3,p4)*zb(p1,p4)
     &    **2*izb(p1,p5)*izb(p1,p6)*izb(p2,p5)*zab2(p6,p2,p5,p1)*
     &    prop34**(-1)*Qsum**2*s256**(-1) + 1.D0/2.D0*za(p3,p4)*zb(p1,
     &    p4)**2*izb(p1,p5)*izb(p1,p6)*izb(p2,p6)*zab2(p5,p2,p6,p1)*
     &    prop34**(-1)*s256**(-1)
      b(1,1) = b(1,1) - za(p3,p4)*zb(p1,p4)**2*izb(p1,p5)*izb(p1,p6)*
     & izb(p2,p6)*zab2(p5,p2,p6,p1)*prop34**(-1)*Qsum*s256**(-1) + 1.D0/
     &    2.D0*za(p3,p4)*zb(p1,p4)**2*izb(p1,p5)*izb(p1,p6)*izb(p2,p6)*
     &    zab2(p5,p2,p6,p1)*prop34**(-1)*Qsum**2*s256**(-1) - za(p3,p5)
     &    *zb(p1,p4)*izb(p1,p6)*izb(p2,p6)*zab2(p5,p2,p6,p1)*
     &    prop34**(-1)*prop126**(-1) + za(p3,p5)*zb(p1,p4)*izb(p1,p6)*
     &    izb(p2,p6)*zab2(p5,p2,p6,p1)*prop34**(-1)*prop126**(-1)*Qsum
     &     - za(p3,p6)*zb(p1,p4)*izb(p1,p5)*izb(p2,p5)*zab2(p6,p2,p5,p1
     &    )*prop34**(-1)*prop125**(-1) + za(p3,p6)*zb(p1,p4)*izb(p1,p5)
     &    *izb(p2,p5)*zab2(p6,p2,p5,p1)*prop34**(-1)*prop125**(-1)*Qsum
     &     + zb(p1,p4)**2*izb(p1,p5)*izb(p1,p6)*izb(p2,p5)*izb(p4,p6)*
     &    zab2(p3,p2,p5,p1)*prop125**(-1) - zb(p1,p4)**2*izb(p1,p5)*
     &    izb(p1,p6)*izb(p2,p5)*izb(p4,p6)*zab2(p3,p2,p5,p1)*
     &    prop125**(-1)*Qsum + zb(p1,p4)**2*izb(p1,p5)*izb(p1,p6)*izb(
     &    p2,p6)*izb(p4,p5)*zab2(p3,p2,p6,p1)*prop126**(-1) - zb(p1,p4)
     &    **2*izb(p1,p5)*izb(p1,p6)*izb(p2,p6)*izb(p4,p5)*zab2(p3,p2,p6
     &    ,p1)*prop126**(-1)*Qsum
      b(1,1) = b(1,1) - zb(p1,p4)*izb(p1,p5)*izb(p1,p6)*izb(p2,p5)*
     & zab2(p3,p2,p5,p1)*zab2(p6,p3,p4,p1)*prop34**(-1)*prop125**(-1)
     &     + zb(p1,p4)*izb(p1,p5)*izb(p1,p6)*izb(p2,p5)*zab2(p3,p2,p5,
     &    p1)*zab2(p6,p3,p4,p1)*prop34**(-1)*prop125**(-1)*Qsum - zb(p1
     &    ,p4)*izb(p1,p5)*izb(p1,p6)*izb(p2,p6)*zab2(p3,p2,p6,p1)*zab2(
     &    p5,p3,p4,p1)*prop34**(-1)*prop126**(-1) + zb(p1,p4)*izb(p1,p5
     &    )*izb(p1,p6)*izb(p2,p6)*zab2(p3,p2,p6,p1)*zab2(p5,p3,p4,p1)*
     &    prop34**(-1)*prop126**(-1)*Qsum
      b(1,2)= - 2*za(p1,p2)*za(p2,p3)*zb(p1,p4)**2*zb(p1,p6)*iza(p2,p6)
     &    *izb(p1,p5)*izb(p4,p5)*prop12**(-1)*prop126**(-1) + 2*za(p1,
     &    p2)*za(p2,p3)*zb(p1,p4)*zb(p1,p6)*iza(p2,p6)*izb(p1,p5)*zab2(
     &    p5,p3,p4,p1)*prop34**(-1)*prop12**(-1)*prop126**(-1) + 2*za(
     &    p1,p2)*za(p3,p5)*zb(p1,p4)*zb(p1,p6)*iza(p2,p6)*izb(p1,p5)*
     &    zab2(p2,p3,p4,p1)*prop34**(-1)*prop12**(-1)*prop126**(-1) + 2
     &    *za(p2,p3)**2*zb(p1,p3)*zb(p1,p4)*zb(p4,p6)*iza(p2,p6)*izb(p1
     &    ,p5)*izb(p4,p5)*prop12**(-1)*s456**(-1) - 2*za(p2,p3)*za(p2,
     &    p5)*zb(p1,p2)*zb(p1,p4)*iza(p2,p6)*izb(p1,p5)*zab2(p2,p3,p4,
     &    p6)*prop34**(-1)*prop12**(-1)*prop125**(-1) - 2*za(p2,p3)*za(
     &    p2,p5)*zb(p1,p2)*zb(p4,p6)*iza(p2,p6)*izb(p1,p5)*zab2(p2,p3,
     &    p4,p1)*prop34**(-1)*prop12**(-1)*prop125**(-1) - 2*za(p2,p3)*
     &    za(p2,p5)*zb(p1,p2)*iza(p2,p6)*iza(p4,p6)*izb(p1,p5)*zab2(p2,
     &    p4,p6,p1)*prop12**(-1)*prop125**(-1) - 2*za(p2,p3)*za(p2,p5)*
     &    zb(p1,p4)*zb(p1,p6)*iza(p2,p6)*izb(p1,p5)*prop34**(-1)*
     &    prop12**(-1)
      b(1,2) = b(1,2) - 2*za(p2,p3)*za(p2,p5)*zb(p4,p6)*iza(p2,p6)*
     & zab2(p5,p3,p4,p1)*prop34**(-1)*prop12**(-1)*prop125**(-1) + za(
     &    p2,p3)*zb(p1,p2)*iza(p2,p6)*iza(p4,p6)*izb(p1,p5)*izb(p2,p5)*
     &    zab2(p2,p4,p6,p1)*prop125**(-1) - za(p2,p3)*zb(p1,p2)*iza(p2,
     &    p6)*iza(p4,p6)*izb(p1,p5)*izb(p2,p5)*zab2(p2,p4,p6,p1)*
     &    prop125**(-1)*Qsum - 2*za(p2,p3)*zb(p1,p4)*zb(p1,p6)*zb(p4,p6
     &    )*izb(p1,p5)*izb(p4,p5)*prop12**(-1)*prop126**(-1) - za(p2,p3
     &    )*zb(p1,p4)*iza(p1,p6)*iza(p2,p6)*izb(p1,p5)*izb(p4,p5)*zab2(
     &    p2,p1,p6,p4)*prop126**(-1) - za(p2,p3)*zb(p1,p4)*iza(p1,p6)*
     &    iza(p2,p6)*izb(p1,p5)*izb(p4,p5)*zab2(p2,p1,p6,p4)*
     &    prop126**(-1)*Qsum + za(p2,p3)*zb(p1,p6)*zb(p4,p6)*iza(p2,p6)
     &    *izb(p1,p5)*izb(p2,p5)*zab2(p6,p2,p5,p1)*prop34**(-1)*
     &    prop125**(-1) - za(p2,p3)*zb(p1,p6)*zb(p4,p6)*iza(p2,p6)*izb(
     &    p1,p5)*izb(p2,p5)*zab2(p6,p2,p5,p1)*prop34**(-1)*
     &    prop125**(-1)*Qsum + 2*za(p2,p3)*zb(p1,p6)*zb(p4,p6)*izb(p1,
     &    p5)*zab2(p5,p3,p4,p1)*prop34**(-1)*prop12**(-1)*prop126**(-1)
      b(1,2) = b(1,2) + 1.D0/2.D0*za(p2,p3)*zb(p1,p6)*iza(p1,p6)*izb(p1
     & ,p5)*zab2(p5,p1,p6,p4)*prop34**(-1)*s156**(-1) + za(p2,p3)*zb(p1
     &    ,p6)*iza(p1,p6)*izb(p1,p5)*zab2(p5,p1,p6,p4)*prop34**(-1)*
     &    Qsum*s156**(-1) + 1.D0/2.D0*za(p2,p3)*zb(p1,p6)*iza(p1,p6)*
     &    izb(p1,p5)*zab2(p5,p1,p6,p4)*prop34**(-1)*Qsum**2*s156**(-1)
     &     + za(p2,p3)*iza(p1,p6)*iza(p2,p6)*izb(p1,p5)*zab2(p2,p1,p6,
     &    p4)*zab2(p5,p3,p4,p1)*prop34**(-1)*prop126**(-1) + za(p2,p3)*
     &    iza(p1,p6)*iza(p2,p6)*izb(p1,p5)*zab2(p2,p1,p6,p4)*zab2(p5,p3
     &    ,p4,p1)*prop34**(-1)*prop126**(-1)*Qsum + 2*za(p2,p3)*iza(p2,
     &    p6)*iza(p4,p6)*izb(p1,p5)*zab2(p2,p4,p6,p1)*zab2(p5,p4,p6,p1)
     &    *prop12**(-1)*s456**(-1) + 2*za(p2,p5)**2*za(p3,p6)*zb(p1,p6)
     &    *zb(p4,p6)*iza(p2,p6)*prop34**(-1)*prop12**(-1)*prop125**(-1)
     &     + za(p2,p5)*za(p3,p5)*zb(p1,p4)*iza(p1,p6)*iza(p2,p6)*izb(p1
     &    ,p5)*zab2(p2,p1,p6,p5)*prop34**(-1)*prop126**(-1) + za(p2,p5)
     &    *za(p3,p5)*zb(p1,p4)*iza(p1,p6)*iza(p2,p6)*izb(p1,p5)*zab2(p2
     &    ,p1,p6,p5)*prop34**(-1)*prop126**(-1)*Qsum
      b(1,2) = b(1,2) + 2*za(p2,p5)*za(p3,p5)*zb(p1,p4)*iza(p2,p6)*
     & zab2(p2,p3,p4,p6)*prop34**(-1)*prop12**(-1)*prop125**(-1) + 2*
     &    za(p2,p5)*za(p3,p5)*zb(p1,p6)**2*zb(p4,p5)*izb(p1,p5)*
     &    prop34**(-1)*prop12**(-1)*prop126**(-1) - za(p2,p5)*za(p3,p5)
     &    *zb(p1,p6)*zb(p4,p5)*iza(p1,p6)*izb(p1,p5)*prop34**(-1)*
     &    prop126**(-1) - za(p2,p5)*za(p3,p5)*zb(p1,p6)*zb(p4,p5)*iza(
     &    p1,p6)*izb(p1,p5)*prop34**(-1)*prop126**(-1)*Qsum + 2*za(p2,
     &    p5)*za(p3,p5)*iza(p2,p6)*iza(p4,p6)*zab2(p2,p4,p6,p1)*
     &    prop12**(-1)*prop125**(-1) - za(p2,p5)*za(p3,p6)*zb(p1,p6)*
     &    zb(p4,p6)*iza(p2,p6)*izb(p2,p5)*prop34**(-1)*prop125**(-1) +
     &    za(p2,p5)*za(p3,p6)*zb(p1,p6)*zb(p4,p6)*iza(p2,p6)*izb(p2,p5)
     &    *prop34**(-1)*prop125**(-1)*Qsum + 1.D0/2.D0*za(p2,p5)*zb(p1,
     &    p4)*iza(p2,p6)*izb(p2,p5)*zab2(p3,p2,p5,p6)*prop34**(-1)*
     &    s256**(-1) - za(p2,p5)*zb(p1,p4)*iza(p2,p6)*izb(p2,p5)*zab2(
     &    p3,p2,p5,p6)*prop34**(-1)*Qsum*s256**(-1) + 1.D0/2.D0*za(p2,
     &    p5)*zb(p1,p4)*iza(p2,p6)*izb(p2,p5)*zab2(p3,p2,p5,p6)*
     &    prop34**(-1)*Qsum**2*s256**(-1)
      b(1,2) = b(1,2) - 2*za(p3,p5)*zb(p1,p4)*zb(p1,p6)*izb(p1,p5)*
     & zab2(p2,p3,p4,p6)*prop34**(-1)*prop12**(-1)*prop126**(-1) - za(
     &    p3,p5)*iza(p2,p6)*iza(p4,p6)*izb(p2,p5)*zab2(p2,p4,p6,p1)*
     &    prop125**(-1) + za(p3,p5)*iza(p2,p6)*iza(p4,p6)*izb(p2,p5)*
     &    zab2(p2,p4,p6,p1)*prop125**(-1)*Qsum + zb(p1,p4)*iza(p2,p6)*
     &    izb(p1,p5)*izb(p2,p5)*zab2(p2,p3,p4,p6)*zab2(p3,p2,p5,p1)*
     &    prop34**(-1)*prop125**(-1) - zb(p1,p4)*iza(p2,p6)*izb(p1,p5)*
     &    izb(p2,p5)*zab2(p2,p3,p4,p6)*zab2(p3,p2,p5,p1)*prop34**(-1)*
     &    prop125**(-1)*Qsum + 1.D0/2.D0*iza(p1,p6)*iza(p2,p6)*izb(p1,
     &    p5)*izb(p2,p5)*zab2(p2,p1,p6,p4)*zab2(p3,p2,p5,p1)*
     &    prop34**(-1) - 1.D0/2.D0*iza(p1,p6)*iza(p2,p6)*izb(p1,p5)*
     &    izb(p2,p5)*zab2(p2,p1,p6,p4)*zab2(p3,p2,p5,p1)*prop34**(-1)*
     &    Qsum**2
      b(2,1)= - 2*za(p1,p2)*za(p2,p3)*zb(p1,p4)**2*zb(p1,p5)*iza(p2,p5)
     &    *izb(p1,p6)*izb(p4,p6)*prop12**(-1)*prop125**(-1) + 2*za(p1,
     &    p2)*za(p2,p3)*zb(p1,p4)*zb(p1,p5)*iza(p2,p5)*izb(p1,p6)*zab2(
     &    p6,p3,p4,p1)*prop34**(-1)*prop12**(-1)*prop125**(-1) + 2*za(
     &    p1,p2)*za(p3,p6)*zb(p1,p4)*zb(p1,p5)*iza(p2,p5)*izb(p1,p6)*
     &    zab2(p2,p3,p4,p1)*prop34**(-1)*prop12**(-1)*prop125**(-1) + 2
     &    *za(p2,p3)**2*zb(p1,p3)*zb(p1,p4)*zb(p4,p5)*iza(p2,p5)*izb(p1
     &    ,p6)*izb(p4,p6)*prop12**(-1)*s456**(-1) - 2*za(p2,p3)*za(p2,
     &    p6)*zb(p1,p2)*zb(p1,p4)*iza(p2,p5)*izb(p1,p6)*zab2(p2,p3,p4,
     &    p5)*prop34**(-1)*prop12**(-1)*prop126**(-1) - 2*za(p2,p3)*za(
     &    p2,p6)*zb(p1,p2)*zb(p4,p5)*iza(p2,p5)*izb(p1,p6)*zab2(p2,p3,
     &    p4,p1)*prop34**(-1)*prop12**(-1)*prop126**(-1) - 2*za(p2,p3)*
     &    za(p2,p6)*zb(p1,p2)*iza(p2,p5)*iza(p4,p5)*izb(p1,p6)*zab2(p2,
     &    p4,p5,p1)*prop12**(-1)*prop126**(-1) - 2*za(p2,p3)*za(p2,p6)*
     &    zb(p1,p4)*zb(p1,p5)*iza(p2,p5)*izb(p1,p6)*prop34**(-1)*
     &    prop12**(-1)
      b(2,1) = b(2,1) - 2*za(p2,p3)*za(p2,p6)*zb(p4,p5)*iza(p2,p5)*
     & zab2(p6,p3,p4,p1)*prop34**(-1)*prop12**(-1)*prop126**(-1) + za(
     &    p2,p3)*zb(p1,p2)*iza(p2,p5)*iza(p4,p5)*izb(p1,p6)*izb(p2,p6)*
     &    zab2(p2,p4,p5,p1)*prop126**(-1) - za(p2,p3)*zb(p1,p2)*iza(p2,
     &    p5)*iza(p4,p5)*izb(p1,p6)*izb(p2,p6)*zab2(p2,p4,p5,p1)*
     &    prop126**(-1)*Qsum - 2*za(p2,p3)*zb(p1,p4)*zb(p1,p5)*zb(p4,p5
     &    )*izb(p1,p6)*izb(p4,p6)*prop12**(-1)*prop125**(-1) - za(p2,p3
     &    )*zb(p1,p4)*iza(p1,p5)*iza(p2,p5)*izb(p1,p6)*izb(p4,p6)*zab2(
     &    p2,p1,p5,p4)*prop125**(-1) - za(p2,p3)*zb(p1,p4)*iza(p1,p5)*
     &    iza(p2,p5)*izb(p1,p6)*izb(p4,p6)*zab2(p2,p1,p5,p4)*
     &    prop125**(-1)*Qsum + za(p2,p3)*zb(p1,p5)*zb(p4,p5)*iza(p2,p5)
     &    *izb(p1,p6)*izb(p2,p6)*zab2(p5,p2,p6,p1)*prop34**(-1)*
     &    prop126**(-1) - za(p2,p3)*zb(p1,p5)*zb(p4,p5)*iza(p2,p5)*izb(
     &    p1,p6)*izb(p2,p6)*zab2(p5,p2,p6,p1)*prop34**(-1)*
     &    prop126**(-1)*Qsum + 2*za(p2,p3)*zb(p1,p5)*zb(p4,p5)*izb(p1,
     &    p6)*zab2(p6,p3,p4,p1)*prop34**(-1)*prop12**(-1)*prop125**(-1)
      b(2,1) = b(2,1) + 1.D0/2.D0*za(p2,p3)*zb(p1,p5)*iza(p1,p5)*izb(p1
     & ,p6)*zab2(p6,p1,p5,p4)*prop34**(-1)*s156**(-1) + za(p2,p3)*zb(p1
     &    ,p5)*iza(p1,p5)*izb(p1,p6)*zab2(p6,p1,p5,p4)*prop34**(-1)*
     &    Qsum*s156**(-1) + 1.D0/2.D0*za(p2,p3)*zb(p1,p5)*iza(p1,p5)*
     &    izb(p1,p6)*zab2(p6,p1,p5,p4)*prop34**(-1)*Qsum**2*s156**(-1)
     &     + za(p2,p3)*iza(p1,p5)*iza(p2,p5)*izb(p1,p6)*zab2(p2,p1,p5,
     &    p4)*zab2(p6,p3,p4,p1)*prop34**(-1)*prop125**(-1) + za(p2,p3)*
     &    iza(p1,p5)*iza(p2,p5)*izb(p1,p6)*zab2(p2,p1,p5,p4)*zab2(p6,p3
     &    ,p4,p1)*prop34**(-1)*prop125**(-1)*Qsum + 2*za(p2,p3)*iza(p2,
     &    p5)*iza(p4,p5)*izb(p1,p6)*zab2(p2,p4,p5,p1)*zab2(p6,p4,p5,p1)
     &    *prop12**(-1)*s456**(-1) + 2*za(p2,p6)**2*za(p3,p5)*zb(p1,p5)
     &    *zb(p4,p5)*iza(p2,p5)*prop34**(-1)*prop12**(-1)*prop126**(-1)
     &     - za(p2,p6)*za(p3,p5)*zb(p1,p5)*zb(p4,p5)*iza(p2,p5)*izb(p2,
     &    p6)*prop34**(-1)*prop126**(-1) + za(p2,p6)*za(p3,p5)*zb(p1,p5
     &    )*zb(p4,p5)*iza(p2,p5)*izb(p2,p6)*prop34**(-1)*prop126**(-1)*
     &    Qsum
      b(2,1) = b(2,1) + za(p2,p6)*za(p3,p6)*zb(p1,p4)*iza(p1,p5)*iza(p2
     & ,p5)*izb(p1,p6)*zab2(p2,p1,p5,p6)*prop34**(-1)*prop125**(-1) +
     &    za(p2,p6)*za(p3,p6)*zb(p1,p4)*iza(p1,p5)*iza(p2,p5)*izb(p1,p6
     &    )*zab2(p2,p1,p5,p6)*prop34**(-1)*prop125**(-1)*Qsum + 2*za(p2
     &    ,p6)*za(p3,p6)*zb(p1,p4)*iza(p2,p5)*zab2(p2,p3,p4,p5)*
     &    prop34**(-1)*prop12**(-1)*prop126**(-1) + 2*za(p2,p6)*za(p3,
     &    p6)*zb(p1,p5)**2*zb(p4,p6)*izb(p1,p6)*prop34**(-1)*
     &    prop12**(-1)*prop125**(-1) - za(p2,p6)*za(p3,p6)*zb(p1,p5)*
     &    zb(p4,p6)*iza(p1,p5)*izb(p1,p6)*prop34**(-1)*prop125**(-1) -
     &    za(p2,p6)*za(p3,p6)*zb(p1,p5)*zb(p4,p6)*iza(p1,p5)*izb(p1,p6)
     &    *prop34**(-1)*prop125**(-1)*Qsum + 2*za(p2,p6)*za(p3,p6)*iza(
     &    p2,p5)*iza(p4,p5)*zab2(p2,p4,p5,p1)*prop12**(-1)*
     &    prop126**(-1) + 1.D0/2.D0*za(p2,p6)*zb(p1,p4)*iza(p2,p5)*izb(
     &    p2,p6)*zab2(p3,p2,p6,p5)*prop34**(-1)*s256**(-1) - za(p2,p6)*
     &    zb(p1,p4)*iza(p2,p5)*izb(p2,p6)*zab2(p3,p2,p6,p5)*
     &    prop34**(-1)*Qsum*s256**(-1)
      b(2,1) = b(2,1) + 1.D0/2.D0*za(p2,p6)*zb(p1,p4)*iza(p2,p5)*izb(p2
     & ,p6)*zab2(p3,p2,p6,p5)*prop34**(-1)*Qsum**2*s256**(-1) - 2*za(p3
     &    ,p6)*zb(p1,p4)*zb(p1,p5)*izb(p1,p6)*zab2(p2,p3,p4,p5)*
     &    prop34**(-1)*prop12**(-1)*prop125**(-1) - za(p3,p6)*iza(p2,p5
     &    )*iza(p4,p5)*izb(p2,p6)*zab2(p2,p4,p5,p1)*prop126**(-1) + za(
     &    p3,p6)*iza(p2,p5)*iza(p4,p5)*izb(p2,p6)*zab2(p2,p4,p5,p1)*
     &    prop126**(-1)*Qsum + zb(p1,p4)*iza(p2,p5)*izb(p1,p6)*izb(p2,
     &    p6)*zab2(p2,p3,p4,p5)*zab2(p3,p2,p6,p1)*prop34**(-1)*
     &    prop126**(-1) - zb(p1,p4)*iza(p2,p5)*izb(p1,p6)*izb(p2,p6)*
     &    zab2(p2,p3,p4,p5)*zab2(p3,p2,p6,p1)*prop34**(-1)*
     &    prop126**(-1)*Qsum + 1.D0/2.D0*iza(p1,p5)*iza(p2,p5)*izb(p1,
     &    p6)*izb(p2,p6)*zab2(p2,p1,p5,p4)*zab2(p3,p2,p6,p1)*
     &    prop34**(-1) - 1.D0/2.D0*iza(p1,p5)*iza(p2,p5)*izb(p1,p6)*
     &    izb(p2,p6)*zab2(p2,p1,p5,p4)*zab2(p3,p2,p6,p1)*prop34**(-1)*
     &    Qsum**2
      b(2,2)= - 2*za(p1,p2)*za(p2,p3)*zb(p1,p4)*zb(p1,p5)*iza(p2,p5)*
     &    iza(p2,p6)*zab2(p2,p3,p4,p6)*prop34**(-1)*prop12**(-1)*
     &    prop125**(-1) - 2*za(p1,p2)*za(p2,p3)*zb(p1,p4)*zb(p1,p6)*
     &    iza(p2,p5)*iza(p2,p6)*zab2(p2,p3,p4,p5)*prop34**(-1)*
     &    prop12**(-1)*prop126**(-1) - 2*za(p1,p2)*za(p2,p3)*zb(p1,p5)*
     &    zb(p4,p6)*iza(p2,p5)*iza(p2,p6)*zab2(p2,p3,p4,p1)*
     &    prop34**(-1)*prop12**(-1)*prop125**(-1) - 2*za(p1,p2)*za(p2,
     &    p3)*zb(p1,p5)*iza(p2,p5)*iza(p2,p6)*iza(p4,p6)*zab2(p2,p4,p6,
     &    p1)*prop12**(-1)*prop125**(-1) - 2*za(p1,p2)*za(p2,p3)*zb(p1,
     &    p6)*zb(p4,p5)*iza(p2,p5)*iza(p2,p6)*zab2(p2,p3,p4,p1)*
     &    prop34**(-1)*prop12**(-1)*prop126**(-1) - 2*za(p1,p2)*za(p2,
     &    p3)*zb(p1,p6)*iza(p2,p5)*iza(p2,p6)*iza(p4,p5)*zab2(p2,p4,p5,
     &    p1)*prop12**(-1)*prop126**(-1) + za(p1,p2)*za(p2,p3)*iza(p1,
     &    p5)*iza(p2,p5)*iza(p2,p6)*iza(p4,p6)*zab2(p2,p4,p6,p1)*
     &    prop125**(-1) + za(p1,p2)*za(p2,p3)*iza(p1,p5)*iza(p2,p5)*
     &    iza(p2,p6)*iza(p4,p6)*zab2(p2,p4,p6,p1)*prop125**(-1)*Qsum
      b(2,2) = b(2,2) + za(p1,p2)*za(p2,p3)*iza(p1,p6)*iza(p2,p5)*iza(
     & p2,p6)*iza(p4,p5)*zab2(p2,p4,p5,p1)*prop126**(-1) + za(p1,p2)*
     &    za(p2,p3)*iza(p1,p6)*iza(p2,p5)*iza(p2,p6)*iza(p4,p5)*zab2(p2
     &    ,p4,p5,p1)*prop126**(-1)*Qsum - 2*za(p2,p3)**2*zb(p1,p3)*iza(
     &    p2,p5)*iza(p2,p6)*iza(p4,p5)*zab2(p2,p4,p5,p6)*prop12**(-1)*
     &    s456**(-1) - 2*za(p2,p3)**2*zb(p1,p3)*iza(p2,p5)*iza(p2,p6)*
     &    iza(p4,p6)*zab2(p2,p4,p6,p5)*prop12**(-1)*s456**(-1) - 2*za(
     &    p2,p3)*zb(p1,p5)*zb(p4,p5)*iza(p2,p6)*zab2(p2,p3,p4,p6)*
     &    prop34**(-1)*prop12**(-1)*prop125**(-1) + 2*za(p2,p3)*zb(p1,
     &    p5)*zb(p4,p6)*iza(p2,p6)*zab2(p2,p3,p4,p5)*prop34**(-1)*
     &    prop12**(-1)*prop125**(-1) + 2*za(p2,p3)*zb(p1,p5)*iza(p2,p6)
     &    *iza(p4,p6)*zab2(p2,p4,p6,p5)*prop12**(-1)*prop125**(-1) + 2*
     &    za(p2,p3)*zb(p1,p6)*zb(p4,p5)*iza(p2,p5)*zab2(p2,p3,p4,p6)*
     &    prop34**(-1)*prop12**(-1)*prop126**(-1) - 2*za(p2,p3)*zb(p1,
     &    p6)*zb(p4,p6)*iza(p2,p5)*zab2(p2,p3,p4,p5)*prop34**(-1)*
     &    prop12**(-1)*prop126**(-1)
      b(2,2) = b(2,2) + 2*za(p2,p3)*zb(p1,p6)*iza(p2,p5)*iza(p4,p5)*
     & zab2(p2,p4,p5,p6)*prop12**(-1)*prop126**(-1) - za(p2,p3)*zb(p4,
     &    p5)*iza(p1,p6)*iza(p2,p6)*zab2(p2,p1,p6,p5)*prop34**(-1)*
     &    prop126**(-1) - za(p2,p3)*zb(p4,p5)*iza(p1,p6)*iza(p2,p6)*
     &    zab2(p2,p1,p6,p5)*prop34**(-1)*prop126**(-1)*Qsum - za(p2,p3)
     &    *zb(p4,p6)*iza(p1,p5)*iza(p2,p5)*zab2(p2,p1,p5,p6)*
     &    prop34**(-1)*prop125**(-1) - za(p2,p3)*zb(p4,p6)*iza(p1,p5)*
     &    iza(p2,p5)*zab2(p2,p1,p5,p6)*prop34**(-1)*prop125**(-1)*Qsum
     &     - za(p2,p3)*iza(p1,p5)*iza(p2,p5)*iza(p2,p6)*zab2(p2,p1,p5,
     &    p4)*zab2(p2,p3,p4,p6)*prop34**(-1)*prop125**(-1) - za(p2,p3)*
     &    iza(p1,p5)*iza(p2,p5)*iza(p2,p6)*zab2(p2,p1,p5,p4)*zab2(p2,p3
     &    ,p4,p6)*prop34**(-1)*prop125**(-1)*Qsum + 1.D0/2.D0*za(p2,p3)
     &    *iza(p1,p5)*iza(p2,p5)*iza(p2,p6)*zab2(p2,p1,p5,p6)*zab3(p2,
     &    p1,p5,p6,p4)*prop34**(-1)*s156**(-1) + za(p2,p3)*iza(p1,p5)*
     &    iza(p2,p5)*iza(p2,p6)*zab2(p2,p1,p5,p6)*zab3(p2,p1,p5,p6,p4)*
     &    prop34**(-1)*Qsum*s156**(-1)
      b(2,2) = b(2,2) + 1.D0/2.D0*za(p2,p3)*iza(p1,p5)*iza(p2,p5)*iza(
     & p2,p6)*zab2(p2,p1,p5,p6)*zab3(p2,p1,p5,p6,p4)*prop34**(-1)*
     & Qsum**2*s156**(-1) - za(p2,p3)*iza(p1,p5)*iza(p2,p6)*iza(p4,p6)*
     &    zab2(p2,p4,p6,p5)*prop125**(-1) - za(p2,p3)*iza(p1,p5)*iza(p2
     &    ,p6)*iza(p4,p6)*zab2(p2,p4,p6,p5)*prop125**(-1)*Qsum - za(p2,
     &    p3)*iza(p1,p6)*iza(p2,p5)*iza(p2,p6)*zab2(p2,p1,p6,p4)*zab2(
     &    p2,p3,p4,p5)*prop34**(-1)*prop126**(-1) - za(p2,p3)*iza(p1,p6
     &    )*iza(p2,p5)*iza(p2,p6)*zab2(p2,p1,p6,p4)*zab2(p2,p3,p4,p5)*
     &    prop34**(-1)*prop126**(-1)*Qsum + 1.D0/2.D0*za(p2,p3)*iza(p1,
     &    p6)*iza(p2,p5)*iza(p2,p6)*zab2(p2,p1,p6,p5)*zab3(p2,p1,p5,p6,
     &    p4)*prop34**(-1)*s156**(-1) + za(p2,p3)*iza(p1,p6)*iza(p2,p5)
     &    *iza(p2,p6)*zab2(p2,p1,p6,p5)*zab3(p2,p1,p5,p6,p4)*
     &    prop34**(-1)*Qsum*s156**(-1) + 1.D0/2.D0*za(p2,p3)*iza(p1,p6)
     &    *iza(p2,p5)*iza(p2,p6)*zab2(p2,p1,p6,p5)*zab3(p2,p1,p5,p6,p4)
     &    *prop34**(-1)*Qsum**2*s156**(-1) - za(p2,p3)*iza(p1,p6)*iza(
     &    p2,p5)*iza(p4,p5)*zab2(p2,p4,p5,p6)*prop126**(-1)
      b(2,2) = b(2,2) - za(p2,p3)*iza(p1,p6)*iza(p2,p5)*iza(p4,p5)*
     & zab2(p2,p4,p5,p6)*prop126**(-1)*Qsum

c      write(6,*) 'b(1,1)',b(1,1)
c      write(6,*) 'b(1,2)',b(1,2)
c      write(6,*) 'b(2,1)',b(2,1)
c      write(6,*) 'b(2,2)',b(2,2)
c      pause

      return
      end
