!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a8WZsr(p1,p2,p3,p4,p5,p6,p7,p8,asr)
      implicit none
c     Overall factor of  i_*gs^2*e^4*(T^A)_{i7,i1}*(T^A)_{i8,i2}
c     Tr{T^A T^B} = delta^{AB}
c     including QED type propagators
c     This is the part of the amplitude with two W's
c     Exceptionally the coupling of this need an extra factor of 1/(2*xw)
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer::p1,p2,p3,p4,p5,p6,p7,p8
      real(dp):: s3,s128,s278,s345,s346,s356,s456,s34,s28,s56,s3456
      complex(dp):: zab2,asr(4:7,2,2)

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p1,p3)

      s34=s(p3,p4)
      s56=s(p5,p6)
      s28=s(p2,p8)
      s128=s3(p1,p2,p8)
      s278=s3(p2,p7,p8)
      s345=s3(p3,p4,p5)
      s346=s3(p3,p4,p6)
      s356=s3(p3,p5,p6)
      s456=s3(p4,p5,p6)
      s3456=s(p3,p4)+s(p3,p5)+s(p3,p6)+s(p4,p5)+s(p4,p6)+s(p5,p6)

c     amplitude electron(4,h28,h56)
      asr(4,1,1)=(za(p7,p8)*zb(p4,p1)/s278
     & *(zb(p6,p3)*zab2(p3,p7,p8,p2)+zb(p6,p5)*zab2(p5,p7,p8,p2))
     & -zab2(p7,p3,p5,p6)*zb(p1,p2)*zab2(p8,p1,p2,p4)/s128)
     & *za(p3,p5)/(s356*s28*s3456*s56)

      asr(4,2,1)=(-za(p2,p7)*zb(p4,p1)/s278
     & *(zb(p6,p3)*zab2(p3,p2,p7,p8)+zb(p6,p5)*zab2(p5,p2,p7,p8))
     & -zab2(p7,p3,p5,p6)*zb(p1,p8)*zab2(p2,p1,p8,p4)/s128)
     & *za(p3,p5)/(s356*s28*s3456*s56)

      asr(4,1,2)=(za(p7,p8)*zb(p4,p1)/s278
     & *(zb(p5,p3)*zab2(p3,p7,p8,p2)+zb(p5,p6)*zab2(p6,p7,p8,p2))
     & -zab2(p7,p3,p6,p5)*zb(p1,p2)*zab2(p8,p1,p2,p4)/s128)
     & *za(p3,p6)/(s356*s28*s3456*s56)

      asr(4,2,2)=(-za(p2,p7)*zb(p4,p1)/s278
     & *(zb(p5,p3)*zab2(p3,p2,p7,p8)+zb(p5,p6)*zab2(p6,p2,p7,p8))
     & -zab2(p7,p3,p6,p5)*zb(p1,p8)*zab2(p2,p1,p8,p4)/s128)
     & *za(p3,p6)/(s356*s28*s3456*s56)

c     amplitude aneutrino(5,h28,h56)
      asr(5,1,1)=(za(p7,p3)*zb(p1,p2)/s128
     & *(za(p4,p5)*zab2(p8,p1,p2,p4)+za(p6,p5)*zab2(p8,p1,p2,p6))
     &- zab2(p5,p4,p6,p1)*za(p7,p8)*zab2(p3,p7,p8,p2)/s278)
     & *zb(p4,p6)/(s456*s28*s3456*s56)

      asr(5,2,1)=(za(p7,p3)*zb(p1,p8)/s128
     & *(za(p4,p5)*zab2(p2,p1,p8,p4)+za(p6,p5)*zab2(p2,p1,p8,p6))
     & +za(p2,p7)*zab2(p5,p4,p6,p1)*zab2(p3,p2,p7,p8)/s278)
     & *zb(p4,p6)/(s456*s28*s3456*s56)

      asr(5,1,2)=(za(p7,p3)*zb(p1,p2)/s128
     & *(zab2(p8,p1,p2,p4)*za(p4,p6)+zab2(p8,p1,p2,p5)*za(p5,p6))
     & -zab2(p6,p4,p5,p1)*za(p7,p8)*zab2(p3,p7,p8,p2)/s278)
     & *zb(p4,p5)/(s456*s28*s3456*s56)

      asr(5,2,2)=(za(p7,p3)*zb(p1,p8)/s128
     & *(zab2(p2,p1,p8,p4)*za(p4,p6)+zab2(p2,p1,p8,p5)*za(p5,p6))
     & +za(p2,p7)*zab2(p6,p4,p5,p1)*zab2(p3,p2,p7,p8)/s278)
     & *zb(p4,p5)/(s456*s28*s3456*s56)

c     amplitude a(6,h28,h56)
      asr(6,1,1)=(
     &   za(p4,p3)*za(p7,p5)*zb(p2,p1)*zab2(p8,p1,p2,p4)/s128
     &  +za(p6,p3)*za(p7,p5)*zb(p2,p1)*zab2(p8,p1,p2,p6)/s128
     &  +za(p7,p8)*zab2(p3,p4,p6,p1)*zab2(p5,p7,p8,p2)/s278)
     & *zb(p4,p6)/(s28*s34*s3456*s346)

      asr(6,2,1)=(za(p7,p5)*zb(p8,p1)/s128
     & *(zab2(p2,p1,p8,p4)*za(p4,p3)+zab2(p2,p1,p8,p6)*za(p6,p3))
     & + za(p7,p2)*zab2(p3,p4,p6,p1)*zab2(p5,p2,p7,p8)/s278)
     & *zb(p4,p6)/(s28*s34*s3456*s346)

      asr(6,:,2)=czip

      asr(7,1,1)=(za(p7,p8)*zb(p6,p1)/s278
     & *(zb(p4,p5)*zab2(p5,p7,p8,p2)+zb(p4,p3)*zab2(p3,p7,p8,p2))
     & -zab2(p7,p5,p3,p4)*zb(p1,p2)*zab2(p8,p1,p2,p6)/s128)
     & *za(p5,p3)/(s345*s28*s3456*s34)

      asr(7,2,1)=(-za(p2,p7)*zb(p6,p1)/s278
     & *(zb(p4,p5)*zab2(p5,p2,p7,p8)+zb(p4,p3)*zab2(p3,p2,p7,p8))
     & -zab2(p7,p5,p3,p4)*zb(p1,p8)*zab2(p2,p1,p8,p6)/s128)
     & *za(p5,p3)/(s345*s28*s3456*s34)

      asr(7,:,2)=czip
      return
      end

