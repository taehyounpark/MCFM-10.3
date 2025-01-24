!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a7ZZn(p1,p2,p3,p4,p5,p6,p7,p,n,za,zb,a7n)
      implicit none
c---- Matrix element for ZZ radiation from line 12
c---- including Z Z interchange
c  q(-p1) +q~(-p2) --> e^-(p3)+e^+(p4)+mu^-(p5)+mu^+(p6)+g(p7)
c     overall factor of i*gs*e^2*2*rt2 removed





c        5-----<-- 6   3-----<--4
c            |Z            |Z
c   2 ----<--|-------------|----1
c                 0
c                 0
c                 0

c     with line 7 contracted with n
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer i,j,p1,p2,p3,p4,p5,p6,p7
      real(dp):: s3,n(4),p(mxpart,4),
     & s134,s234,s156,s256,s34,s56
      complex(dp):: a7n(2,2,2),vecm(mxpart,mxpart),zab2,iza,izb

      s3(p1,p2,p3)=s(p1,p2)+s(p1,p3)+s(p2,p3)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      iza(p1,p2)=cone/za(p1,p2)
      izb(p1,p2)=cone/zb(p1,p2)

      call checkndotp(p,n,p7)
      s34=s(p3,p4)
      s56=s(p5,p6)
      s134=s3(p1,p3,p4)
      s234=s3(p2,p3,p4)
      s156=s3(p1,p5,p6)
      s256=s3(p2,p5,p6)
      do i=1,7
      do j=i,7
      call ndveccur(i,j,n,p,vecm)
      enddo
      enddo

      a7n(1,1,1)= + s156**(-1)*s234**(-1) * ( za(p1,p5)*za(p2,p3)*zb(p1
     &    ,p6)*zb(p2,p4)*vecm(p2,p1) + za(p1,p5)*za(p2,p3)*zb(p1,p6)*
     &    zb(p3,p4)*vecm(p3,p1) - za(p2,p3)*za(p5,p6)*zb(p1,p6)*zb(p2,
     &    p4)*vecm(p2,p6) - za(p2,p3)*za(p5,p6)*zb(p1,p6)*zb(p3,p4)*
     &    vecm(p3,p6) )
      a7n(1,1,1) = a7n(1,1,1) + s256**(-1)*s134**(-1) * ( za(p1,p3)*za(
     &    p2,p5)*zb(p1,p4)*zb(p2,p6)*vecm(p2,p1) + za(p1,p3)*za(p2,p5)*
     &    zb(p1,p4)*zb(p5,p6)*vecm(p5,p1) - za(p2,p5)*za(p3,p4)*zb(p1,
     &    p4)*zb(p2,p6)*vecm(p2,p4) - za(p2,p5)*za(p3,p4)*zb(p1,p4)*zb(
     &    p5,p6)*vecm(p5,p4) )
      a7n(1,1,1) = a7n(1,1,1) + iza(p1,p7)*izb(p1,p7)*s256**(-1) * (
     &    za(p2,p5)*zb(p1,p4)*zab2(p3,p2,p5,p6)*vecm(p1,p1) - za(p2,p5)
     &    *zb(p4,p7)*zab2(p3,p2,p5,p6)*vecm(p7,p1) )
      a7n(1,1,1) = a7n(1,1,1) + iza(p1,p7)*izb(p1,p7)*s234**(-1) * (
     &    za(p2,p3)*zb(p1,p6)*zab2(p5,p2,p3,p4)*vecm(p1,p1) - za(p2,p3)
     &    *zb(p6,p7)*zab2(p5,p2,p3,p4)*vecm(p7,p1) )
      a7n(1,1,1) = a7n(1,1,1) + iza(p2,p7)*izb(p2,p7)*s156**(-1) * (
     &    za(p2,p3)*zb(p1,p6)*zab2(p5,p1,p6,p4)*vecm(p2,p2) - za(p3,p7)
     &    *zb(p1,p6)*zab2(p5,p1,p6,p4)*vecm(p2,p7) )
      a7n(1,1,1) = a7n(1,1,1) + iza(p2,p7)*izb(p2,p7)*s134**(-1) * (
     &    za(p2,p5)*zb(p1,p4)*zab2(p3,p1,p4,p6)*vecm(p2,p2) - za(p5,p7)
     &    *zb(p1,p4)*zab2(p3,p1,p4,p6)*vecm(p2,p7) )
      a7n(1,2,1)= + s156**(-1)*s234**(-1) * ( za(p1,p5)*za(p2,p4)*zb(p1
     &    ,p6)*zb(p2,p3)*vecm(p2,p1) - za(p1,p5)*za(p2,p4)*zb(p1,p6)*
     &    zb(p3,p4)*vecm(p4,p1) - za(p2,p4)*za(p5,p6)*zb(p1,p6)*zb(p2,
     &    p3)*vecm(p2,p6) + za(p2,p4)*za(p5,p6)*zb(p1,p6)*zb(p3,p4)*
     &    vecm(p4,p6) )
      a7n(1,2,1) = a7n(1,2,1) + s256**(-1)*s134**(-1) * ( za(p1,p4)*za(
     &    p2,p5)*zb(p1,p3)*zb(p2,p6)*vecm(p2,p1) + za(p1,p4)*za(p2,p5)*
     &    zb(p1,p3)*zb(p5,p6)*vecm(p5,p1) + za(p2,p5)*za(p3,p4)*zb(p1,
     &    p3)*zb(p2,p6)*vecm(p2,p3) + za(p2,p5)*za(p3,p4)*zb(p1,p3)*zb(
     &    p5,p6)*vecm(p5,p3) )
      a7n(1,2,1) = a7n(1,2,1) + iza(p1,p7)*izb(p1,p7)*s256**(-1) * (
     &    za(p2,p5)*zb(p1,p3)*zab2(p4,p2,p5,p6)*vecm(p1,p1) - za(p2,p5)
     &    *zb(p3,p7)*zab2(p4,p2,p5,p6)*vecm(p7,p1) )
      a7n(1,2,1) = a7n(1,2,1) + iza(p1,p7)*izb(p1,p7)*s234**(-1) * (
     &    za(p2,p4)*zb(p1,p6)*zab2(p5,p2,p4,p3)*vecm(p1,p1) - za(p2,p4)
     &    *zb(p6,p7)*zab2(p5,p2,p4,p3)*vecm(p7,p1) )
      a7n(1,2,1) = a7n(1,2,1) + iza(p2,p7)*izb(p2,p7)*s156**(-1) * (
     &    za(p2,p4)*zb(p1,p6)*zab2(p5,p1,p6,p3)*vecm(p2,p2) - za(p4,p7)
     &    *zb(p1,p6)*zab2(p5,p1,p6,p3)*vecm(p2,p7) )
      a7n(1,2,1) = a7n(1,2,1) + iza(p2,p7)*izb(p2,p7)*s134**(-1) * (
     &    za(p2,p5)*zb(p1,p3)*zab2(p4,p1,p3,p6)*vecm(p2,p2) - za(p5,p7)
     &    *zb(p1,p3)*zab2(p4,p1,p3,p6)*vecm(p2,p7) )
      a7n(2,1,1)= + s156**(-1)*s234**(-1) * ( za(p1,p5)*za(p2,p3)*zb(p1
     &    ,p6)*zb(p2,p4)*vecm(p1,p2) + za(p1,p5)*za(p2,p3)*zb(p2,p4)*
     &    zb(p5,p6)*vecm(p5,p2) - za(p1,p5)*za(p3,p4)*zb(p1,p6)*zb(p2,
     &    p4)*vecm(p1,p4) - za(p1,p5)*za(p3,p4)*zb(p2,p4)*zb(p5,p6)*
     &    vecm(p5,p4) )
      a7n(2,1,1) = a7n(2,1,1) + s256**(-1)*s134**(-1) * ( za(p1,p3)*za(
     &    p2,p5)*zb(p1,p4)*zb(p2,p6)*vecm(p1,p2) + za(p1,p3)*za(p2,p5)*
     &    zb(p2,p6)*zb(p3,p4)*vecm(p3,p2) - za(p1,p3)*za(p5,p6)*zb(p1,
     &    p4)*zb(p2,p6)*vecm(p1,p6) - za(p1,p3)*za(p5,p6)*zb(p2,p6)*zb(
     &    p3,p4)*vecm(p3,p6) )
      a7n(2,1,1) = a7n(2,1,1) + iza(p1,p7)*izb(p1,p7)*s256**(-1) * (
     &    za(p1,p3)*zb(p2,p6)*zab2(p5,p2,p6,p4)*vecm(p1,p1) - za(p3,p7)
     &    *zb(p2,p6)*zab2(p5,p2,p6,p4)*vecm(p1,p7) )
      a7n(2,1,1) = a7n(2,1,1) + iza(p1,p7)*izb(p1,p7)*s234**(-1) * (
     &    za(p1,p5)*zb(p2,p4)*zab2(p3,p2,p4,p6)*vecm(p1,p1) - za(p5,p7)
     &    *zb(p2,p4)*zab2(p3,p2,p4,p6)*vecm(p1,p7) )
      a7n(2,1,1) = a7n(2,1,1) + iza(p2,p7)*izb(p2,p7)*s156**(-1) * (
     &    za(p1,p5)*zb(p2,p4)*zab2(p3,p1,p5,p6)*vecm(p2,p2) - za(p1,p5)
     &    *zb(p4,p7)*zab2(p3,p1,p5,p6)*vecm(p7,p2) )
      a7n(2,1,1) = a7n(2,1,1) + iza(p2,p7)*izb(p2,p7)*s134**(-1) * (
     &    za(p1,p3)*zb(p2,p6)*zab2(p5,p1,p3,p4)*vecm(p2,p2) - za(p1,p3)
     &    *zb(p6,p7)*zab2(p5,p1,p3,p4)*vecm(p7,p2) )
      a7n(2,2,1)= + s156**(-1)*s234**(-1) * ( za(p1,p5)*za(p2,p4)*zb(p1
     &    ,p6)*zb(p2,p3)*vecm(p1,p2) + za(p1,p5)*za(p2,p4)*zb(p2,p3)*
     &    zb(p5,p6)*vecm(p5,p2) + za(p1,p5)*za(p3,p4)*zb(p1,p6)*zb(p2,
     &    p3)*vecm(p1,p3) + za(p1,p5)*za(p3,p4)*zb(p2,p3)*zb(p5,p6)*
     &    vecm(p5,p3) )
      a7n(2,2,1) = a7n(2,2,1) + s256**(-1)*s134**(-1) * ( za(p1,p4)*za(
     &    p2,p5)*zb(p1,p3)*zb(p2,p6)*vecm(p1,p2) - za(p1,p4)*za(p2,p5)*
     &    zb(p2,p6)*zb(p3,p4)*vecm(p4,p2) - za(p1,p4)*za(p5,p6)*zb(p1,
     &    p3)*zb(p2,p6)*vecm(p1,p6) + za(p1,p4)*za(p5,p6)*zb(p2,p6)*zb(
     &    p3,p4)*vecm(p4,p6) )
      a7n(2,2,1) = a7n(2,2,1) + iza(p1,p7)*izb(p1,p7)*s256**(-1) * (
     &    za(p1,p4)*zb(p2,p6)*zab2(p5,p2,p6,p3)*vecm(p1,p1) - za(p4,p7)
     &    *zb(p2,p6)*zab2(p5,p2,p6,p3)*vecm(p1,p7) )
      a7n(2,2,1) = a7n(2,2,1) + iza(p1,p7)*izb(p1,p7)*s234**(-1) * (
     &    za(p1,p5)*zb(p2,p3)*zab2(p4,p2,p3,p6)*vecm(p1,p1) - za(p5,p7)
     &    *zb(p2,p3)*zab2(p4,p2,p3,p6)*vecm(p1,p7) )
      a7n(2,2,1) = a7n(2,2,1) + iza(p2,p7)*izb(p2,p7)*s156**(-1) * (
     &    za(p1,p5)*zb(p2,p3)*zab2(p4,p1,p5,p6)*vecm(p2,p2) - za(p1,p5)
     &    *zb(p3,p7)*zab2(p4,p1,p5,p6)*vecm(p7,p2) )
      a7n(2,2,1) = a7n(2,2,1) + iza(p2,p7)*izb(p2,p7)*s134**(-1) * (
     &    za(p1,p4)*zb(p2,p6)*zab2(p5,p1,p4,p3)*vecm(p2,p2) - za(p1,p4)
     &    *zb(p6,p7)*zab2(p5,p1,p4,p3)*vecm(p7,p2) )
      a7n(1,1,2)= + s156**(-1)*s234**(-1) * ( za(p1,p6)*za(p2,p3)*zb(p1
     &    ,p5)*zb(p2,p4)*vecm(p2,p1) + za(p1,p6)*za(p2,p3)*zb(p1,p5)*
     &    zb(p3,p4)*vecm(p3,p1) + za(p2,p3)*za(p5,p6)*zb(p1,p5)*zb(p2,
     &    p4)*vecm(p2,p5) + za(p2,p3)*za(p5,p6)*zb(p1,p5)*zb(p3,p4)*
     &    vecm(p3,p5) )
      a7n(1,1,2) = a7n(1,1,2) + s256**(-1)*s134**(-1) * ( za(p1,p3)*za(
     &    p2,p6)*zb(p1,p4)*zb(p2,p5)*vecm(p2,p1) - za(p1,p3)*za(p2,p6)*
     &    zb(p1,p4)*zb(p5,p6)*vecm(p6,p1) - za(p2,p6)*za(p3,p4)*zb(p1,
     &    p4)*zb(p2,p5)*vecm(p2,p4) + za(p2,p6)*za(p3,p4)*zb(p1,p4)*zb(
     &    p5,p6)*vecm(p6,p4) )
      a7n(1,1,2) = a7n(1,1,2) + iza(p1,p7)*izb(p1,p7)*s256**(-1) * (
     &    za(p2,p6)*zb(p1,p4)*zab2(p3,p2,p6,p5)*vecm(p1,p1) - za(p2,p6)
     &    *zb(p4,p7)*zab2(p3,p2,p6,p5)*vecm(p7,p1) )
      a7n(1,1,2) = a7n(1,1,2) + iza(p1,p7)*izb(p1,p7)*s234**(-1) * (
     &    za(p2,p3)*zb(p1,p5)*zab2(p6,p2,p3,p4)*vecm(p1,p1) - za(p2,p3)
     &    *zb(p5,p7)*zab2(p6,p2,p3,p4)*vecm(p7,p1) )
      a7n(1,1,2) = a7n(1,1,2) + iza(p2,p7)*izb(p2,p7)*s156**(-1) * (
     &    za(p2,p3)*zb(p1,p5)*zab2(p6,p1,p5,p4)*vecm(p2,p2) - za(p3,p7)
     &    *zb(p1,p5)*zab2(p6,p1,p5,p4)*vecm(p2,p7) )
      a7n(1,1,2) = a7n(1,1,2) + iza(p2,p7)*izb(p2,p7)*s134**(-1) * (
     &    za(p2,p6)*zb(p1,p4)*zab2(p3,p1,p4,p5)*vecm(p2,p2) - za(p6,p7)
     &    *zb(p1,p4)*zab2(p3,p1,p4,p5)*vecm(p2,p7) )
      a7n(1,2,2)= + s156**(-1)*s234**(-1) * ( za(p1,p6)*za(p2,p4)*zb(p1
     &    ,p5)*zb(p2,p3)*vecm(p2,p1) - za(p1,p6)*za(p2,p4)*zb(p1,p5)*
     &    zb(p3,p4)*vecm(p4,p1) + za(p2,p4)*za(p5,p6)*zb(p1,p5)*zb(p2,
     &    p3)*vecm(p2,p5) - za(p2,p4)*za(p5,p6)*zb(p1,p5)*zb(p3,p4)*
     &    vecm(p4,p5) )
      a7n(1,2,2) = a7n(1,2,2) + s256**(-1)*s134**(-1) * ( za(p1,p4)*za(
     &    p2,p6)*zb(p1,p3)*zb(p2,p5)*vecm(p2,p1) - za(p1,p4)*za(p2,p6)*
     &    zb(p1,p3)*zb(p5,p6)*vecm(p6,p1) + za(p2,p6)*za(p3,p4)*zb(p1,
     &    p3)*zb(p2,p5)*vecm(p2,p3) - za(p2,p6)*za(p3,p4)*zb(p1,p3)*zb(
     &    p5,p6)*vecm(p6,p3) )
      a7n(1,2,2) = a7n(1,2,2) + iza(p1,p7)*izb(p1,p7)*s256**(-1) * (
     &    za(p2,p6)*zb(p1,p3)*zab2(p4,p2,p6,p5)*vecm(p1,p1) - za(p2,p6)
     &    *zb(p3,p7)*zab2(p4,p2,p6,p5)*vecm(p7,p1) )
      a7n(1,2,2) = a7n(1,2,2) + iza(p1,p7)*izb(p1,p7)*s234**(-1) * (
     &    za(p2,p4)*zb(p1,p5)*zab2(p6,p2,p4,p3)*vecm(p1,p1) - za(p2,p4)
     &    *zb(p5,p7)*zab2(p6,p2,p4,p3)*vecm(p7,p1) )
      a7n(1,2,2) = a7n(1,2,2) + iza(p2,p7)*izb(p2,p7)*s156**(-1) * (
     &    za(p2,p4)*zb(p1,p5)*zab2(p6,p1,p5,p3)*vecm(p2,p2) - za(p4,p7)
     &    *zb(p1,p5)*zab2(p6,p1,p5,p3)*vecm(p2,p7) )
      a7n(1,2,2) = a7n(1,2,2) + iza(p2,p7)*izb(p2,p7)*s134**(-1) * (
     &    za(p2,p6)*zb(p1,p3)*zab2(p4,p1,p3,p5)*vecm(p2,p2) - za(p6,p7)
     &    *zb(p1,p3)*zab2(p4,p1,p3,p5)*vecm(p2,p7) )
      a7n(2,1,2)= + s156**(-1)*s234**(-1) * ( za(p1,p6)*za(p2,p3)*zb(p1
     &    ,p5)*zb(p2,p4)*vecm(p1,p2) - za(p1,p6)*za(p2,p3)*zb(p2,p4)*
     &    zb(p5,p6)*vecm(p6,p2) - za(p1,p6)*za(p3,p4)*zb(p1,p5)*zb(p2,
     &    p4)*vecm(p1,p4) + za(p1,p6)*za(p3,p4)*zb(p2,p4)*zb(p5,p6)*
     &    vecm(p6,p4) )
      a7n(2,1,2) = a7n(2,1,2) + s256**(-1)*s134**(-1) * ( za(p1,p3)*za(
     &    p2,p6)*zb(p1,p4)*zb(p2,p5)*vecm(p1,p2) + za(p1,p3)*za(p2,p6)*
     &    zb(p2,p5)*zb(p3,p4)*vecm(p3,p2) + za(p1,p3)*za(p5,p6)*zb(p1,
     &    p4)*zb(p2,p5)*vecm(p1,p5) + za(p1,p3)*za(p5,p6)*zb(p2,p5)*zb(
     &    p3,p4)*vecm(p3,p5) )
      a7n(2,1,2) = a7n(2,1,2) + iza(p1,p7)*izb(p1,p7)*s256**(-1) * (
     &    za(p1,p3)*zb(p2,p5)*zab2(p6,p2,p5,p4)*vecm(p1,p1) - za(p3,p7)
     &    *zb(p2,p5)*zab2(p6,p2,p5,p4)*vecm(p1,p7) )
      a7n(2,1,2) = a7n(2,1,2) + iza(p1,p7)*izb(p1,p7)*s234**(-1) * (
     &    za(p1,p6)*zb(p2,p4)*zab2(p3,p2,p4,p5)*vecm(p1,p1) - za(p6,p7)
     &    *zb(p2,p4)*zab2(p3,p2,p4,p5)*vecm(p1,p7) )
      a7n(2,1,2) = a7n(2,1,2) + iza(p2,p7)*izb(p2,p7)*s156**(-1) * (
     &    za(p1,p6)*zb(p2,p4)*zab2(p3,p1,p6,p5)*vecm(p2,p2) - za(p1,p6)
     &    *zb(p4,p7)*zab2(p3,p1,p6,p5)*vecm(p7,p2) )
      a7n(2,1,2) = a7n(2,1,2) + iza(p2,p7)*izb(p2,p7)*s134**(-1) * (
     &    za(p1,p3)*zb(p2,p5)*zab2(p6,p1,p3,p4)*vecm(p2,p2) - za(p1,p3)
     &    *zb(p5,p7)*zab2(p6,p1,p3,p4)*vecm(p7,p2) )
      a7n(2,2,2)= + s156**(-1)*s234**(-1) * ( za(p1,p6)*za(p2,p4)*zb(p1
     &    ,p5)*zb(p2,p3)*vecm(p1,p2) - za(p1,p6)*za(p2,p4)*zb(p2,p3)*
     &    zb(p5,p6)*vecm(p6,p2) + za(p1,p6)*za(p3,p4)*zb(p1,p5)*zb(p2,
     &    p3)*vecm(p1,p3) - za(p1,p6)*za(p3,p4)*zb(p2,p3)*zb(p5,p6)*
     &    vecm(p6,p3) )
      a7n(2,2,2) = a7n(2,2,2) + s256**(-1)*s134**(-1) * ( za(p1,p4)*za(
     &    p2,p6)*zb(p1,p3)*zb(p2,p5)*vecm(p1,p2) - za(p1,p4)*za(p2,p6)*
     &    zb(p2,p5)*zb(p3,p4)*vecm(p4,p2) + za(p1,p4)*za(p5,p6)*zb(p1,
     &    p3)*zb(p2,p5)*vecm(p1,p5) - za(p1,p4)*za(p5,p6)*zb(p2,p5)*zb(
     &    p3,p4)*vecm(p4,p5) )
      a7n(2,2,2) = a7n(2,2,2) + iza(p1,p7)*izb(p1,p7)*s256**(-1) * (
     &    za(p1,p4)*zb(p2,p5)*zab2(p6,p2,p5,p3)*vecm(p1,p1) - za(p4,p7)
     &    *zb(p2,p5)*zab2(p6,p2,p5,p3)*vecm(p1,p7) )
      a7n(2,2,2) = a7n(2,2,2) + iza(p1,p7)*izb(p1,p7)*s234**(-1) * (
     &    za(p1,p6)*zb(p2,p3)*zab2(p4,p2,p3,p5)*vecm(p1,p1) - za(p6,p7)
     &    *zb(p2,p3)*zab2(p4,p2,p3,p5)*vecm(p1,p7) )
      a7n(2,2,2) = a7n(2,2,2) + iza(p2,p7)*izb(p2,p7)*s156**(-1) * (
     &    za(p1,p6)*zb(p2,p3)*zab2(p4,p1,p6,p5)*vecm(p2,p2) - za(p1,p6)
     &    *zb(p3,p7)*zab2(p4,p1,p6,p5)*vecm(p7,p2) )
      a7n(2,2,2) = a7n(2,2,2) + iza(p2,p7)*izb(p2,p7)*s134**(-1) * (
     &    za(p1,p4)*zb(p2,p5)*zab2(p6,p1,p4,p3)*vecm(p2,p2) - za(p1,p4)
     &    *zb(p5,p7)*zab2(p6,p1,p4,p3)*vecm(p7,p2) )
      a7n(:,:,:)=a7n(:,:,:)/(s34*s56)
      return
      end
