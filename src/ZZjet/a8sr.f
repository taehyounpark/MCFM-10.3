!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a8sr(p1,p2,p3,p4,p5,p6,p7,p8,
     & A56AB,A56BA,A34AB,A34BA)
      implicit none
      include 'types.f'
      integer p1,p2,p3,p4,p5,p6,p7,p8,h1,h2
      complex(dp):: A56AB(2,2,2,2,2),A56BA(2,2,2,2,2),
     & A34AB(2,2,2,2,2),A34BA(2,2,2,2,2),temp(2,2,2,2,2)
c     a56(2,2,2,2,2) equiv f(polg1,polg2,polq,pol34,pol56)
c     result for the ZZ process supplementary diagrams
c     with a factor of 4*e^4*gs^2 removed
c     setting up color orders andcontribution 34<-->56
c     Amplitudes for the process
c     0 --> qbar(j1)+q(j2)+e^-(j3)+e^+(p4)+mu^-(p5)+mu^+(p6)+g(p7)+g(p8)
c     only s34,s56,propagators included.
c     pol=1=LH, pol=2=RH
      call a8srpart(p1,p2,p3,p4,p5,p6,p7,p8,A56AB)
      call a8srpart(p1,p2,p3,p4,p5,p6,p8,p7,A56BA)

      call a8srpart(p1,p2,p5,p6,p3,p4,p7,p8,temp)
      do h1=1,2
      do h2=1,2
      A34AB(:,:,:,h1,h2)= temp(:,:,:,h2,h1)
      enddo
      enddo

      call a8srpart(p1,p2,p5,p6,p3,p4,p8,p7,temp)
      do h1=1,2
      do h2=1,2
      A34BA(:,:,:,h1,h2)=temp(:,:,:,h2,h1)
      enddo
      enddo

      return
      end



      subroutine a8srpart(p1,p2,p3,p4,j5,j6,p7,p8,A56)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer::p1,p2,p3,p4,p5,p6,p7,p8
      integer::j5,j6,h56
      real(dp):: s356,s456,s3456,s178,s278,s56,s78,s3
      complex(dp):: A56(2,2,2,2,2),zab2,zab3,iza,izb
c     a56(2,2,2,2,2) equiv f(polg1,polg2,polq,pol34,pol56)
c     result for the ZZ process supplementary diagrams
c     with a factor of 4*e^4*gs^2 removed
c     one colour ordering only
c     contribution 34<-->56 needs to be calculated
c     Amplitudes for the process
c     0 --> qbar(j1)+q(j2)+e^-(j3)+e^+(p4)+mu^-(p5)+mu^+(p6)+g(p7)+g(p8)
c     only s34,s56,propagators included.
c     pol=1=LH, pol=2=RH

      s3(p1,p2,p3)=s(p1,p2)+s(p1,p3)+s(p2,p3)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zab3(p1,p2,p3,p4,p5)=
     & za(p1,p2)*zb(p2,p5)+za(p1,p3)*zb(p3,p5)+za(p1,p4)*zb(p4,p5)
      iza(p1,p2)=1d0/za(p1,p2)
      izb(p1,p2)=1d0/zb(p1,p2)
      s78=s(p7,p8)
      s56=s(j5,j6)
      s356=s3(p3,j5,j6)
      s456=s3(p4,j5,j6)
      s178=s3(p1,p7,p8)
      s278=s3(p2,p7,p8)
      s3456=s(p3,p4)+s(p3,j5)+s(p3,j6)+s(p4,j5)+s(p4,j6)+s(j5,j6)
c      do h56=1,2
c      if (h56 == 1) then
      p5=j5
      p6=j6
c      elseif (h56 == 2) then
c      p5=j6
c      p6=j5
c      endif
      h56=1
      A56(1,1,1,1,h56)= + s3456**(-1)*s278**(-1)*s356**(-1) * ( za(p2,
     &    p7)*za(p3,p5)*zb(p1,p4)*zb(p3,p6)*izb(p1,p8)*izb(p7,p8)*zab3(
     &    p3,p2,p7,p8,p1)*s56**(-1) - za(p2,p7)*za(p3,p5)*zb(p1,p4)*
     &    iza(p5,p6)*izb(p1,p8)*izb(p7,p8)*zab3(p5,p2,p7,p8,p1) + za(p2
     &    ,p8)*za(p3,p5)*zb(p1,p4)*zb(p3,p6)*izb(p1,p7)*izb(p7,p8)*
     &    zab3(p3,p2,p7,p8,p1)*s56**(-1) - za(p2,p8)*za(p3,p5)*zb(p1,p4
     &    )*iza(p5,p6)*izb(p1,p7)*izb(p7,p8)*zab3(p5,p2,p7,p8,p1) + za(
     &    p3,p5)*zb(p1,p4)*zb(p3,p6)*izb(p1,p7)*izb(p1,p8)*izb(p2,p7)*
     &    zab2(p8,p2,p7,p1)*zab3(p3,p2,p7,p8,p1)*s56**(-1) - za(p3,p5)*
     &    zb(p1,p4)*iza(p5,p6)*izb(p1,p7)*izb(p1,p8)*izb(p2,p7)*zab2(p8
     &    ,p2,p7,p1)*zab3(p5,p2,p7,p8,p1) )
      A56(1,1,1,1,h56) = A56(1,1,1,1,h56) + s3456**(-1)*s278**(-1)*
     & s456**(-1) * (  - za(p2,p7)*za(p4,p5)*zb(p1,p4)*zb(p4,p6)*izb(p1
     &    ,p8)*izb(p7,p8)*zab3(p3,p2,p7,p8,p1)*s56**(-1) - za(p2,p7)*
     &    zb(p1,p6)*zb(p4,p6)*izb(p1,p8)*izb(p5,p6)*izb(p7,p8)*zab3(p3,
     &    p2,p7,p8,p1) - za(p2,p8)*za(p4,p5)*zb(p1,p4)*zb(p4,p6)*izb(p1
     &    ,p7)*izb(p7,p8)*zab3(p3,p2,p7,p8,p1)*s56**(-1) - za(p2,p8)*
     &    zb(p1,p6)*zb(p4,p6)*izb(p1,p7)*izb(p5,p6)*izb(p7,p8)*zab3(p3,
     &    p2,p7,p8,p1) - za(p4,p5)*zb(p1,p4)*zb(p4,p6)*izb(p1,p7)*izb(
     &    p1,p8)*izb(p2,p7)*zab2(p8,p2,p7,p1)*zab3(p3,p2,p7,p8,p1)*
     &    s56**(-1) - zb(p1,p6)*zb(p4,p6)*izb(p1,p7)*izb(p1,p8)*izb(p2,
     &    p7)*izb(p5,p6)*zab2(p8,p2,p7,p1)*zab3(p3,p2,p7,p8,p1) )
      A56(1,1,2,1,h56)= + s3456**(-1)*s178**(-1)*s356**(-1) * (  - za(
     &    p1,p7)*za(p3,p5)*zb(p2,p4)*zb(p3,p6)*izb(p1,p8)*izb(p7,p8)*
     &    zab3(p3,p1,p7,p8,p1)*s56**(-1) + za(p1,p7)*za(p3,p5)*zb(p2,p4
     &    )*iza(p5,p6)*izb(p1,p8)*izb(p7,p8)*zab3(p5,p1,p7,p8,p1) - za(
     &    p1,p8)*za(p3,p5)*zb(p2,p4)*zb(p3,p6)*izb(p1,p7)*izb(p7,p8)*
     &    zab3(p3,p1,p7,p8,p1)*s56**(-1) + za(p1,p8)*za(p3,p5)*zb(p2,p4
     &    )*iza(p5,p6)*izb(p1,p7)*izb(p7,p8)*zab3(p5,p1,p7,p8,p1) + za(
     &    p3,p5)*zb(p2,p4)*zb(p3,p6)*izb(p1,p7)*izb(p1,p8)**2*zab2(p7,
     &    p1,p8,p1)*zab3(p3,p1,p7,p8,p1)*s56**(-1) - za(p3,p5)*zb(p2,p4
     &    )*iza(p5,p6)*izb(p1,p7)*izb(p1,p8)**2*zab2(p7,p1,p8,p1)*zab3(
     &    p5,p1,p7,p8,p1) )
      A56(1,1,2,1,h56) = A56(1,1,2,1,h56) + s3456**(-1)*s178**(-1)*
     & s456**(-1) * ( za(p1,p7)*za(p4,p5)*zb(p2,p4)*zb(p4,p6)*izb(p1,p8
     &    )*izb(p7,p8)*zab3(p3,p1,p7,p8,p1)*s56**(-1) + za(p1,p7)*zb(p2
     &    ,p6)*zb(p4,p6)*izb(p1,p8)*izb(p5,p6)*izb(p7,p8)*zab3(p3,p1,p7
     &    ,p8,p1) + za(p1,p8)*za(p4,p5)*zb(p2,p4)*zb(p4,p6)*izb(p1,p7)*
     &    izb(p7,p8)*zab3(p3,p1,p7,p8,p1)*s56**(-1) + za(p1,p8)*zb(p2,
     &    p6)*zb(p4,p6)*izb(p1,p7)*izb(p5,p6)*izb(p7,p8)*zab3(p3,p1,p7,
     &    p8,p1) - za(p4,p5)*zb(p2,p4)*zb(p4,p6)*izb(p1,p7)*izb(p1,p8)
     &    **2*zab2(p7,p1,p8,p1)*zab3(p3,p1,p7,p8,p1)*s56**(-1) - zb(p2,
     &    p6)*zb(p4,p6)*izb(p1,p7)*izb(p1,p8)**2*izb(p5,p6)*zab2(p7,p1,
     &    p8,p1)*zab3(p3,p1,p7,p8,p1) )
      A56(1,1,2,1,h56) = A56(1,1,2,1,h56) + s3456**(-1)*s278**(-1)*
     & s356**(-1) * (  - za(p1,p3)*za(p3,p5)*zb(p1,p2)**2*zb(p3,p6)*
     &    izb(p1,p7)*izb(p1,p8)*izb(p2,p7)*zab2(p8,p2,p7,p4)*s56**(-1)
     &     - za(p1,p3)*za(p3,p5)*zb(p1,p2)*zb(p3,p6)*izb(p1,p7)*izb(p7,
     &    p8)*zab2(p8,p2,p7,p4)*s56**(-1) - za(p1,p3)*za(p3,p5)*zb(p1,
     &    p2)*zb(p3,p6)*izb(p1,p8)*izb(p7,p8)*zab2(p7,p2,p8,p4)*
     &    s56**(-1) + za(p1,p5)*za(p3,p5)*zb(p1,p2)**2*iza(p5,p6)*izb(
     &    p1,p7)*izb(p1,p8)*izb(p2,p7)*zab2(p8,p2,p7,p4) + za(p1,p5)*
     &    za(p3,p5)*zb(p1,p2)*iza(p5,p6)*izb(p1,p7)*izb(p7,p8)*zab2(p8,
     &    p2,p7,p4) + za(p1,p5)*za(p3,p5)*zb(p1,p2)*iza(p5,p6)*izb(p1,
     &    p8)*izb(p7,p8)*zab2(p7,p2,p8,p4) )
      A56(1,1,2,1,h56) = A56(1,1,2,1,h56) + s3456**(-1)*s278**(-1)*
     & s456**(-1) * ( za(p1,p3)*za(p4,p5)*zb(p1,p2)**2*zb(p4,p6)*izb(p1
     &    ,p7)*izb(p1,p8)*izb(p2,p7)*zab2(p8,p2,p7,p4)*s56**(-1) + za(
     &    p1,p3)*za(p4,p5)*zb(p1,p2)*zb(p4,p6)*izb(p1,p7)*izb(p7,p8)*
     &    zab2(p8,p2,p7,p4)*s56**(-1) + za(p1,p3)*za(p4,p5)*zb(p1,p2)*
     &    zb(p4,p6)*izb(p1,p8)*izb(p7,p8)*zab2(p7,p2,p8,p4)*s56**(-1)
     &     + za(p1,p3)*zb(p1,p2)**2*zb(p4,p6)*izb(p1,p7)*izb(p1,p8)*
     &    izb(p2,p7)*izb(p5,p6)*zab2(p8,p2,p7,p6) + za(p1,p3)*zb(p1,p2)
     &    *zb(p4,p6)*izb(p1,p7)*izb(p5,p6)*izb(p7,p8)*zab2(p8,p2,p7,p6)
     &     + za(p1,p3)*zb(p1,p2)*zb(p4,p6)*izb(p1,p8)*izb(p5,p6)*izb(p7
     &    ,p8)*zab2(p7,p2,p8,p6) )
      A56(1,1,2,1,h56) = A56(1,1,2,1,h56) + s3456**(-1)*s356**(-1) * (
     &    za(p3,p5)*zb(p1,p2)*zb(p2,p4)*zb(p3,p6)*izb(p1,p7)*izb(p1,p8)
     &    **2*izb(p2,p7)*zab2(p3,p1,p8,p1)*s56**(-1) - za(p3,p5)*zb(p1,
     &    p2)*zb(p2,p4)*iza(p5,p6)*izb(p1,p7)*izb(p1,p8)**2*izb(p2,p7)*
     &    zab2(p5,p1,p8,p1) )
      A56(1,1,2,1,h56) = A56(1,1,2,1,h56) + s3456**(-1)*s456**(-1) * (
     &     - za(p4,p5)*zb(p1,p2)*zb(p2,p4)*zb(p4,p6)*izb(p1,p7)*izb(p1,
     &    p8)**2*izb(p2,p7)*zab2(p3,p1,p8,p1)*s56**(-1) - zb(p1,p2)*zb(
     &    p2,p6)*zb(p4,p6)*izb(p1,p7)*izb(p1,p8)**2*izb(p2,p7)*izb(p5,
     &    p6)*zab2(p3,p1,p8,p1) )
      A56(1,1,1,2,h56)= + s3456**(-1)*s278**(-1)*s356**(-1) * ( za(p2,
     &    p7)*za(p3,p5)*zb(p1,p3)*zb(p3,p6)*izb(p1,p8)*izb(p7,p8)*zab3(
     &    p4,p2,p7,p8,p1)*s56**(-1) + za(p2,p7)*zb(p1,p6)*zb(p3,p6)*
     &    izb(p1,p8)*izb(p5,p6)*izb(p7,p8)*zab3(p4,p2,p7,p8,p1) + za(p2
     &    ,p8)*za(p3,p5)*zb(p1,p3)*zb(p3,p6)*izb(p1,p7)*izb(p7,p8)*
     &    zab3(p4,p2,p7,p8,p1)*s56**(-1) + za(p2,p8)*zb(p1,p6)*zb(p3,p6
     &    )*izb(p1,p7)*izb(p5,p6)*izb(p7,p8)*zab3(p4,p2,p7,p8,p1) + za(
     &    p3,p5)*zb(p1,p3)*zb(p3,p6)*izb(p1,p7)*izb(p1,p8)*izb(p2,p7)*
     &    zab2(p8,p2,p7,p1)*zab3(p4,p2,p7,p8,p1)*s56**(-1) + zb(p1,p6)*
     &    zb(p3,p6)*izb(p1,p7)*izb(p1,p8)*izb(p2,p7)*izb(p5,p6)*zab2(p8
     &    ,p2,p7,p1)*zab3(p4,p2,p7,p8,p1) )
      A56(1,1,1,2,h56) = A56(1,1,1,2,h56) + s3456**(-1)*s278**(-1)*
     & s456**(-1) * (  - za(p2,p7)*za(p4,p5)*zb(p1,p3)*zb(p4,p6)*izb(p1
     &    ,p8)*izb(p7,p8)*zab3(p4,p2,p7,p8,p1)*s56**(-1) + za(p2,p7)*
     &    za(p4,p5)*zb(p1,p3)*iza(p5,p6)*izb(p1,p8)*izb(p7,p8)*zab3(p5,
     &    p2,p7,p8,p1) - za(p2,p8)*za(p4,p5)*zb(p1,p3)*zb(p4,p6)*izb(p1
     &    ,p7)*izb(p7,p8)*zab3(p4,p2,p7,p8,p1)*s56**(-1) + za(p2,p8)*
     &    za(p4,p5)*zb(p1,p3)*iza(p5,p6)*izb(p1,p7)*izb(p7,p8)*zab3(p5,
     &    p2,p7,p8,p1) - za(p4,p5)*zb(p1,p3)*zb(p4,p6)*izb(p1,p7)*izb(
     &    p1,p8)*izb(p2,p7)*zab2(p8,p2,p7,p1)*zab3(p4,p2,p7,p8,p1)*
     &    s56**(-1) + za(p4,p5)*zb(p1,p3)*iza(p5,p6)*izb(p1,p7)*izb(p1,
     &    p8)*izb(p2,p7)*zab2(p8,p2,p7,p1)*zab3(p5,p2,p7,p8,p1) )
      A56(1,1,2,2,h56)= + s3456**(-1)*s178**(-1)*s356**(-1) * (  - za(
     &    p1,p7)*za(p3,p5)*zb(p2,p3)*zb(p3,p6)*izb(p2,p8)*izb(p7,p8)*
     &    zab3(p4,p1,p7,p8,p2)*s56**(-1) - za(p1,p7)*zb(p2,p6)*zb(p3,p6
     &    )*izb(p2,p8)*izb(p5,p6)*izb(p7,p8)*zab3(p4,p1,p7,p8,p2) - za(
     &    p1,p8)*za(p3,p5)*zb(p2,p3)*zb(p3,p6)*izb(p2,p7)*izb(p7,p8)*
     &    zab3(p4,p1,p7,p8,p2)*s56**(-1) - za(p1,p8)*zb(p2,p6)*zb(p3,p6
     &    )*izb(p2,p7)*izb(p5,p6)*izb(p7,p8)*zab3(p4,p1,p7,p8,p2) + za(
     &    p3,p5)*zb(p2,p3)*zb(p3,p6)*izb(p1,p8)*izb(p2,p7)*izb(p2,p8)*
     &    zab2(p7,p1,p8,p2)*zab3(p4,p1,p7,p8,p2)*s56**(-1) + zb(p2,p6)*
     &    zb(p3,p6)*izb(p1,p8)*izb(p2,p7)*izb(p2,p8)*izb(p5,p6)*zab2(p7
     &    ,p1,p8,p2)*zab3(p4,p1,p7,p8,p2) )
      A56(1,1,2,2,h56) = A56(1,1,2,2,h56) + s3456**(-1)*s178**(-1)*
     & s456**(-1) * ( za(p1,p7)*za(p4,p5)*zb(p2,p3)*zb(p4,p6)*izb(p2,p8
     &    )*izb(p7,p8)*zab3(p4,p1,p7,p8,p2)*s56**(-1) - za(p1,p7)*za(p4
     &    ,p5)*zb(p2,p3)*iza(p5,p6)*izb(p2,p8)*izb(p7,p8)*zab3(p5,p1,p7
     &    ,p8,p2) + za(p1,p8)*za(p4,p5)*zb(p2,p3)*zb(p4,p6)*izb(p2,p7)*
     &    izb(p7,p8)*zab3(p4,p1,p7,p8,p2)*s56**(-1) - za(p1,p8)*za(p4,
     &    p5)*zb(p2,p3)*iza(p5,p6)*izb(p2,p7)*izb(p7,p8)*zab3(p5,p1,p7,
     &    p8,p2) - za(p4,p5)*zb(p2,p3)*zb(p4,p6)*izb(p1,p8)*izb(p2,p7)*
     &    izb(p2,p8)*zab2(p7,p1,p8,p2)*zab3(p4,p1,p7,p8,p2)*s56**(-1)
     &     + za(p4,p5)*zb(p2,p3)*iza(p5,p6)*izb(p1,p8)*izb(p2,p7)*izb(
     &    p2,p8)*zab2(p7,p1,p8,p2)*zab3(p5,p1,p7,p8,p2) )

      A56(1,2,1,1,h56)= + s3456**(-1)*s178**(-1)*s356**(-1) * ( za(p1,
     &    p7)*za(p2,p3)*za(p3,p5)*zb(p1,p8)*zb(p3,p6)*iza(p1,p8)*zab2(
     &    p7,p1,p8,p4)*s56**(-1)*s78**(-1) - za(p1,p7)*za(p2,p5)*za(p3,
     &    p5)*zb(p1,p8)*iza(p1,p8)*iza(p5,p6)*zab2(p7,p1,p8,p4)*
     &    s78**(-1) )
      A56(1,2,1,1,h56) = A56(1,2,1,1,h56) + s3456**(-1)*s178**(-1)*
     & s456**(-1) * (  - za(p1,p7)*za(p2,p3)*za(p4,p5)*zb(p1,p8)*zb(p4,
     &    p6)*iza(p1,p8)*zab2(p7,p1,p8,p4)*s56**(-1)*s78**(-1) - za(p1,
     &    p7)*za(p2,p3)*zb(p1,p8)*zb(p4,p6)*iza(p1,p8)*izb(p5,p6)*zab2(
     &    p7,p1,p8,p6)*s78**(-1) )
      A56(1,2,1,1,h56) = A56(1,2,1,1,h56) + s3456**(-1)*s278**(-1)*
     & s356**(-1) * ( za(p2,p7)*za(p3,p5)*zb(p1,p4)*zb(p2,p8)*zb(p3,p6)
     &    *izb(p2,p7)*zab2(p3,p2,p7,p8)*s56**(-1)*s78**(-1) - za(p2,p7)
     &    *za(p3,p5)*zb(p1,p4)*zb(p2,p8)*iza(p5,p6)*izb(p2,p7)*zab2(p5,
     &    p2,p7,p8)*s78**(-1) )
      A56(1,2,1,1,h56) = A56(1,2,1,1,h56) + s3456**(-1)*s278**(-1)*
     & s456**(-1) * (  - za(p2,p7)*za(p4,p5)*zb(p1,p4)*zb(p2,p8)*zb(p4,
     &    p6)*izb(p2,p7)*zab2(p3,p2,p7,p8)*s56**(-1)*s78**(-1) - za(p2,
     &    p7)*zb(p1,p6)*zb(p2,p8)*zb(p4,p6)*izb(p2,p7)*izb(p5,p6)*zab2(
     &    p3,p2,p7,p8)*s78**(-1) )
      A56(1,2,1,1,h56) = A56(1,2,1,1,h56) + s3456**(-1)*s356**(-1) * (
     &     - za(p3,p5)*zb(p3,p6)*iza(p1,p8)*izb(p2,p7)*zab2(p3,p2,p7,p8
     &    )*zab2(p7,p1,p8,p4)*s56**(-1)*s78**(-1) + za(p3,p5)*iza(p1,p8
     &    )*iza(p5,p6)*izb(p2,p7)*zab2(p5,p2,p7,p8)*zab2(p7,p1,p8,p4)*
     &    s78**(-1) )
      A56(1,2,1,1,h56) = A56(1,2,1,1,h56) + s3456**(-1)*s456**(-1) * (
     &    za(p4,p5)*zb(p4,p6)*iza(p1,p8)*izb(p2,p7)*zab2(p3,p2,p7,p8)*
     &    zab2(p7,p1,p8,p4)*s56**(-1)*s78**(-1) + zb(p4,p6)*iza(p1,p8)*
     &    izb(p2,p7)*izb(p5,p6)*zab2(p3,p2,p7,p8)*zab2(p7,p1,p8,p6)*
     &    s78**(-1) )
      A56(1,2,2,1,h56)= + s3456**(-1)*s178**(-1)*s356**(-1) * ( za(p1,
     &    p7)**2*za(p3,p5)*zb(p2,p4)*zb(p3,p6)*iza(p1,p8)*zab2(p3,p1,p7
     &    ,p8)*s56**(-1)*s78**(-1) - za(p1,p7)**2*za(p3,p5)*zb(p2,p4)*
     &    iza(p1,p8)*iza(p5,p6)*zab2(p5,p1,p7,p8)*s78**(-1) )
      A56(1,2,2,1,h56) = A56(1,2,2,1,h56) + s3456**(-1)*s178**(-1)*
     & s456**(-1) * (  - za(p1,p7)**2*za(p4,p5)*zb(p2,p4)*zb(p4,p6)*
     &    iza(p1,p8)*zab2(p3,p1,p7,p8)*s56**(-1)*s78**(-1) - za(p1,p7)
     &    **2*zb(p2,p6)*zb(p4,p6)*iza(p1,p8)*izb(p5,p6)*zab2(p3,p1,p7,
     &    p8)*s78**(-1) )
      A56(1,2,2,1,h56) = A56(1,2,2,1,h56) + s3456**(-1)*s278**(-1)*
     & s356**(-1) * ( za(p1,p3)*za(p3,p5)*zb(p2,p8)**2*zb(p3,p6)*izb(p2
     &    ,p7)*zab2(p7,p2,p8,p4)*s56**(-1)*s78**(-1) - za(p1,p5)*za(p3,
     &    p5)*zb(p2,p8)**2*iza(p5,p6)*izb(p2,p7)*zab2(p7,p2,p8,p4)*
     &    s78**(-1) )
      A56(1,2,2,1,h56) = A56(1,2,2,1,h56) + s3456**(-1)*s278**(-1)*
     & s456**(-1) * (  - za(p1,p3)*za(p4,p5)*zb(p2,p8)**2*zb(p4,p6)*
     &    izb(p2,p7)*zab2(p7,p2,p8,p4)*s56**(-1)*s78**(-1) - za(p1,p3)*
     &    zb(p2,p8)**2*zb(p4,p6)*izb(p2,p7)*izb(p5,p6)*zab2(p7,p2,p8,p6
     &    )*s78**(-1) )
      A56(1,2,2,1,h56) = A56(1,2,2,1,h56) + s3456**(-1)*s356**(-1) * (
     &     - za(p1,p3)*za(p1,p7)*za(p3,p5)*zb(p2,p4)*zb(p2,p8)*zb(p3,p6
     &    )*iza(p1,p8)*izb(p2,p7)*s56**(-1)*s78**(-1) + za(p1,p5)*za(p1
     &    ,p7)*za(p3,p5)*zb(p2,p4)*zb(p2,p8)*iza(p1,p8)*iza(p5,p6)*izb(
     &    p2,p7)*s78**(-1) )
      A56(1,2,2,1,h56) = A56(1,2,2,1,h56) + s3456**(-1)*s456**(-1) * (
     &    za(p1,p3)*za(p1,p7)*za(p4,p5)*zb(p2,p4)*zb(p2,p8)*zb(p4,p6)*
     &    iza(p1,p8)*izb(p2,p7)*s56**(-1)*s78**(-1) + za(p1,p3)*za(p1,
     &    p7)*zb(p2,p6)*zb(p2,p8)*zb(p4,p6)*iza(p1,p8)*izb(p2,p7)*izb(
     &    p5,p6)*s78**(-1) )
      A56(1,2,1,2,h56)= + s3456**(-1)*s178**(-1)*s356**(-1) * ( za(p1,
     &    p7)*za(p2,p4)*za(p3,p5)*zb(p1,p8)*zb(p3,p6)*iza(p1,p8)*zab2(
     &    p7,p1,p8,p3)*s56**(-1)*s78**(-1) + za(p1,p7)*za(p2,p4)*zb(p1,
     &    p8)*zb(p3,p6)*iza(p1,p8)*izb(p5,p6)*zab2(p7,p1,p8,p6)*
     &    s78**(-1) )
      A56(1,2,1,2,h56) = A56(1,2,1,2,h56) + s3456**(-1)*s178**(-1)*
     & s456**(-1) * (  - za(p1,p7)*za(p2,p4)*za(p4,p5)*zb(p1,p8)*zb(p4,
     &    p6)*iza(p1,p8)*zab2(p7,p1,p8,p3)*s56**(-1)*s78**(-1) + za(p1,
     &    p7)*za(p2,p5)*za(p4,p5)*zb(p1,p8)*iza(p1,p8)*iza(p5,p6)*zab2(
     &    p7,p1,p8,p3)*s78**(-1) )
      A56(1,2,1,2,h56) = A56(1,2,1,2,h56) + s3456**(-1)*s278**(-1)*
     & s356**(-1) * ( za(p2,p7)*za(p3,p5)*zb(p1,p3)*zb(p2,p8)*zb(p3,p6)
     &    *izb(p2,p7)*zab2(p4,p2,p7,p8)*s56**(-1)*s78**(-1) + za(p2,p7)
     &    *zb(p1,p6)*zb(p2,p8)*zb(p3,p6)*izb(p2,p7)*izb(p5,p6)*zab2(p4,
     &    p2,p7,p8)*s78**(-1) )
      A56(1,2,1,2,h56) = A56(1,2,1,2,h56) + s3456**(-1)*s278**(-1)*
     & s456**(-1) * (  - za(p2,p7)*za(p4,p5)*zb(p1,p3)*zb(p2,p8)*zb(p4,
     &    p6)*izb(p2,p7)*zab2(p4,p2,p7,p8)*s56**(-1)*s78**(-1) + za(p2,
     &    p7)*za(p4,p5)*zb(p1,p3)*zb(p2,p8)*iza(p5,p6)*izb(p2,p7)*zab2(
     &    p5,p2,p7,p8)*s78**(-1) )
      A56(1,2,1,2,h56) = A56(1,2,1,2,h56) + s3456**(-1)*s356**(-1) * (
     &     - za(p3,p5)*zb(p3,p6)*iza(p1,p8)*izb(p2,p7)*zab2(p4,p2,p7,p8
     &    )*zab2(p7,p1,p8,p3)*s56**(-1)*s78**(-1) - zb(p3,p6)*iza(p1,p8
     &    )*izb(p2,p7)*izb(p5,p6)*zab2(p4,p2,p7,p8)*zab2(p7,p1,p8,p6)*
     &    s78**(-1) )
      A56(1,2,1,2,h56) = A56(1,2,1,2,h56) + s3456**(-1)*s456**(-1) * (
     &    za(p4,p5)*zb(p4,p6)*iza(p1,p8)*izb(p2,p7)*zab2(p4,p2,p7,p8)*
     &    zab2(p7,p1,p8,p3)*s56**(-1)*s78**(-1) - za(p4,p5)*iza(p1,p8)*
     &    iza(p5,p6)*izb(p2,p7)*zab2(p5,p2,p7,p8)*zab2(p7,p1,p8,p3)*
     &    s78**(-1) )
      A56(1,2,2,2,h56)= + s3456**(-1)*s178**(-1)*s356**(-1) * ( za(p1,
     &    p7)**2*za(p3,p5)*zb(p2,p3)*zb(p3,p6)*iza(p1,p8)*zab2(p4,p1,p7
     &    ,p8)*s56**(-1)*s78**(-1) + za(p1,p7)**2*zb(p2,p6)*zb(p3,p6)*
     &    iza(p1,p8)*izb(p5,p6)*zab2(p4,p1,p7,p8)*s78**(-1) )
      A56(1,2,2,2,h56) = A56(1,2,2,2,h56) + s3456**(-1)*s178**(-1)*
     & s456**(-1) * (  - za(p1,p7)**2*za(p4,p5)*zb(p2,p3)*zb(p4,p6)*
     &    iza(p1,p8)*zab2(p4,p1,p7,p8)*s56**(-1)*s78**(-1) + za(p1,p7)
     &    **2*za(p4,p5)*zb(p2,p3)*iza(p1,p8)*iza(p5,p6)*zab2(p5,p1,p7,
     &    p8)*s78**(-1) )
      A56(1,2,2,2,h56) = A56(1,2,2,2,h56) + s3456**(-1)*s278**(-1)*
     & s356**(-1) * ( za(p1,p4)*za(p3,p5)*zb(p2,p8)**2*zb(p3,p6)*izb(p2
     &    ,p7)*zab2(p7,p2,p8,p3)*s56**(-1)*s78**(-1) + za(p1,p4)*zb(p2,
     &    p8)**2*zb(p3,p6)*izb(p2,p7)*izb(p5,p6)*zab2(p7,p2,p8,p6)*
     &    s78**(-1) )
      A56(1,2,2,2,h56) = A56(1,2,2,2,h56) + s3456**(-1)*s278**(-1)*
     & s456**(-1) * (  - za(p1,p4)*za(p4,p5)*zb(p2,p8)**2*zb(p4,p6)*
     &    izb(p2,p7)*zab2(p7,p2,p8,p3)*s56**(-1)*s78**(-1) + za(p1,p5)*
     &    za(p4,p5)*zb(p2,p8)**2*iza(p5,p6)*izb(p2,p7)*zab2(p7,p2,p8,p3
     &    )*s78**(-1) )
      A56(1,2,2,2,h56) = A56(1,2,2,2,h56) + s3456**(-1)*s356**(-1) * (
     &     - za(p1,p4)*za(p1,p7)*za(p3,p5)*zb(p2,p3)*zb(p2,p8)*zb(p3,p6
     &    )*iza(p1,p8)*izb(p2,p7)*s56**(-1)*s78**(-1) - za(p1,p4)*za(p1
     &    ,p7)*zb(p2,p6)*zb(p2,p8)*zb(p3,p6)*iza(p1,p8)*izb(p2,p7)*izb(
     &    p5,p6)*s78**(-1) )
      A56(1,2,2,2,h56) = A56(1,2,2,2,h56) + s3456**(-1)*s456**(-1) * (
     &    za(p1,p4)*za(p1,p7)*za(p4,p5)*zb(p2,p3)*zb(p2,p8)*zb(p4,p6)*
     &    iza(p1,p8)*izb(p2,p7)*s56**(-1)*s78**(-1) - za(p1,p5)*za(p1,
     &    p7)*za(p4,p5)*zb(p2,p3)*zb(p2,p8)*iza(p1,p8)*iza(p5,p6)*izb(
     &    p2,p7)*s78**(-1) )

      A56(2,1,1,1,h56)= + s3456**(-1)*s178**(-1)*s356**(-1) * ( za(p2,
     &    p3)*za(p3,p5)*zb(p1,p7)**2*zb(p3,p6)*izb(p1,p8)*zab2(p8,p1,p7
     &    ,p4)*s56**(-1)*s78**(-1) - za(p2,p5)*za(p3,p5)*zb(p1,p7)**2*
     &    iza(p5,p6)*izb(p1,p8)*zab2(p8,p1,p7,p4)*s78**(-1) )
      A56(2,1,1,1,h56) = A56(2,1,1,1,h56) + s3456**(-1)*s178**(-1)*
     & s456**(-1) * (  - za(p2,p3)*za(p4,p5)*zb(p1,p7)**2*zb(p4,p6)*
     &    izb(p1,p8)*zab2(p8,p1,p7,p4)*s56**(-1)*s78**(-1) - za(p2,p3)*
     &    zb(p1,p7)**2*zb(p4,p6)*izb(p1,p8)*izb(p5,p6)*zab2(p8,p1,p7,p6
     &    )*s78**(-1) )
      A56(2,1,1,1,h56) = A56(2,1,1,1,h56) + s3456**(-1)*s278**(-1)*
     & s356**(-1) * ( za(p2,p8)**2*za(p3,p5)*zb(p1,p4)*zb(p3,p6)*iza(p2
     &    ,p7)*zab2(p3,p2,p8,p7)*s56**(-1)*s78**(-1) - za(p2,p8)**2*za(
     &    p3,p5)*zb(p1,p4)*iza(p2,p7)*iza(p5,p6)*zab2(p5,p2,p8,p7)*
     &    s78**(-1) )
      A56(2,1,1,1,h56) = A56(2,1,1,1,h56) + s3456**(-1)*s278**(-1)*
     & s456**(-1) * (  - za(p2,p8)**2*za(p4,p5)*zb(p1,p4)*zb(p4,p6)*
     &    iza(p2,p7)*zab2(p3,p2,p8,p7)*s56**(-1)*s78**(-1) - za(p2,p8)
     &    **2*zb(p1,p6)*zb(p4,p6)*iza(p2,p7)*izb(p5,p6)*zab2(p3,p2,p8,
     &    p7)*s78**(-1) )
      A56(2,1,1,1,h56) = A56(2,1,1,1,h56) + s3456**(-1)*s356**(-1) * (
     &     - za(p2,p3)*za(p2,p8)*za(p3,p5)*zb(p1,p4)*zb(p1,p7)*zb(p3,p6
     &    )*iza(p2,p7)*izb(p1,p8)*s56**(-1)*s78**(-1) + za(p2,p5)*za(p2
     &    ,p8)*za(p3,p5)*zb(p1,p4)*zb(p1,p7)*iza(p2,p7)*iza(p5,p6)*izb(
     &    p1,p8)*s78**(-1) )
      A56(2,1,1,1,h56) = A56(2,1,1,1,h56) + s3456**(-1)*s456**(-1) * (
     &    za(p2,p3)*za(p2,p8)*za(p4,p5)*zb(p1,p4)*zb(p1,p7)*zb(p4,p6)*
     &    iza(p2,p7)*izb(p1,p8)*s56**(-1)*s78**(-1) + za(p2,p3)*za(p2,
     &    p8)*zb(p1,p6)*zb(p1,p7)*zb(p4,p6)*iza(p2,p7)*izb(p1,p8)*izb(
     &    p5,p6)*s78**(-1) )
      A56(2,1,2,1,h56)= + s3456**(-1)*s178**(-1)*s356**(-1) * ( za(p1,
     &    p8)*za(p3,p5)*zb(p1,p7)*zb(p2,p4)*zb(p3,p6)*izb(p1,p8)*zab2(
     &    p3,p1,p8,p7)*s56**(-1)*s78**(-1) - za(p1,p8)*za(p3,p5)*zb(p1,
     &    p7)*zb(p2,p4)*iza(p5,p6)*izb(p1,p8)*zab2(p5,p1,p8,p7)*
     &    s78**(-1) )
      A56(2,1,2,1,h56) = A56(2,1,2,1,h56) + s3456**(-1)*s178**(-1)*
     & s456**(-1) * (  - za(p1,p8)*za(p4,p5)*zb(p1,p7)*zb(p2,p4)*zb(p4,
     &    p6)*izb(p1,p8)*zab2(p3,p1,p8,p7)*s56**(-1)*s78**(-1) - za(p1,
     &    p8)*zb(p1,p7)*zb(p2,p6)*zb(p4,p6)*izb(p1,p8)*izb(p5,p6)*zab2(
     &    p3,p1,p8,p7)*s78**(-1) )
      A56(2,1,2,1,h56) = A56(2,1,2,1,h56) + s3456**(-1)*s278**(-1)*
     & s356**(-1) * ( za(p1,p3)*za(p2,p8)*za(p3,p5)*zb(p2,p7)*zb(p3,p6)
     &    *iza(p2,p7)*zab2(p8,p2,p7,p4)*s56**(-1)*s78**(-1) - za(p1,p5)
     &    *za(p2,p8)*za(p3,p5)*zb(p2,p7)*iza(p2,p7)*iza(p5,p6)*zab2(p8,
     &    p2,p7,p4)*s78**(-1) )
      A56(2,1,2,1,h56) = A56(2,1,2,1,h56) + s3456**(-1)*s278**(-1)*
     & s456**(-1) * (  - za(p1,p3)*za(p2,p8)*za(p4,p5)*zb(p2,p7)*zb(p4,
     &    p6)*iza(p2,p7)*zab2(p8,p2,p7,p4)*s56**(-1)*s78**(-1) - za(p1,
     &    p3)*za(p2,p8)*zb(p2,p7)*zb(p4,p6)*iza(p2,p7)*izb(p5,p6)*zab2(
     &    p8,p2,p7,p6)*s78**(-1) )
      A56(2,1,2,1,h56) = A56(2,1,2,1,h56) + s3456**(-1)*s356**(-1) * (
     &     - za(p3,p5)*zb(p3,p6)*iza(p2,p7)*izb(p1,p8)*zab2(p3,p1,p8,p7
     &    )*zab2(p8,p2,p7,p4)*s56**(-1)*s78**(-1) + za(p3,p5)*iza(p2,p7
     &    )*iza(p5,p6)*izb(p1,p8)*zab2(p5,p1,p8,p7)*zab2(p8,p2,p7,p4)*
     &    s78**(-1) )
      A56(2,1,2,1,h56) = A56(2,1,2,1,h56) + s3456**(-1)*s456**(-1) * (
     &    za(p4,p5)*zb(p4,p6)*iza(p2,p7)*izb(p1,p8)*zab2(p3,p1,p8,p7)*
     &    zab2(p8,p2,p7,p4)*s56**(-1)*s78**(-1) + zb(p4,p6)*iza(p2,p7)*
     &    izb(p1,p8)*izb(p5,p6)*zab2(p3,p1,p8,p7)*zab2(p8,p2,p7,p6)*
     &    s78**(-1) )
      A56(2,1,1,2,h56)= + s3456**(-1)*s178**(-1)*s356**(-1) * ( za(p2,
     &    p4)*za(p3,p5)*zb(p1,p7)**2*zb(p3,p6)*izb(p1,p8)*zab2(p8,p1,p7
     &    ,p3)*s56**(-1)*s78**(-1) + za(p2,p4)*zb(p1,p7)**2*zb(p3,p6)*
     &    izb(p1,p8)*izb(p5,p6)*zab2(p8,p1,p7,p6)*s78**(-1) )
      A56(2,1,1,2,h56) = A56(2,1,1,2,h56) + s3456**(-1)*s178**(-1)*
     & s456**(-1) * (  - za(p2,p4)*za(p4,p5)*zb(p1,p7)**2*zb(p4,p6)*
     &    izb(p1,p8)*zab2(p8,p1,p7,p3)*s56**(-1)*s78**(-1) + za(p2,p5)*
     &    za(p4,p5)*zb(p1,p7)**2*iza(p5,p6)*izb(p1,p8)*zab2(p8,p1,p7,p3
     &    )*s78**(-1) )
      A56(2,1,1,2,h56) = A56(2,1,1,2,h56) + s3456**(-1)*s278**(-1)*
     & s356**(-1) * ( za(p2,p8)**2*za(p3,p5)*zb(p1,p3)*zb(p3,p6)*iza(p2
     &    ,p7)*zab2(p4,p2,p8,p7)*s56**(-1)*s78**(-1) + za(p2,p8)**2*zb(
     &    p1,p6)*zb(p3,p6)*iza(p2,p7)*izb(p5,p6)*zab2(p4,p2,p8,p7)*
     &    s78**(-1) )
      A56(2,1,1,2,h56) = A56(2,1,1,2,h56) + s3456**(-1)*s278**(-1)*
     & s456**(-1) * (  - za(p2,p8)**2*za(p4,p5)*zb(p1,p3)*zb(p4,p6)*
     &    iza(p2,p7)*zab2(p4,p2,p8,p7)*s56**(-1)*s78**(-1) + za(p2,p8)
     &    **2*za(p4,p5)*zb(p1,p3)*iza(p2,p7)*iza(p5,p6)*zab2(p5,p2,p8,
     &    p7)*s78**(-1) )
      A56(2,1,1,2,h56) = A56(2,1,1,2,h56) + s3456**(-1)*s356**(-1) * (
     &     - za(p2,p4)*za(p2,p8)*za(p3,p5)*zb(p1,p3)*zb(p1,p7)*zb(p3,p6
     &    )*iza(p2,p7)*izb(p1,p8)*s56**(-1)*s78**(-1) - za(p2,p4)*za(p2
     &    ,p8)*zb(p1,p6)*zb(p1,p7)*zb(p3,p6)*iza(p2,p7)*izb(p1,p8)*izb(
     &    p5,p6)*s78**(-1) )
      A56(2,1,1,2,h56) = A56(2,1,1,2,h56) + s3456**(-1)*s456**(-1) * (
     &    za(p2,p4)*za(p2,p8)*za(p4,p5)*zb(p1,p3)*zb(p1,p7)*zb(p4,p6)*
     &    iza(p2,p7)*izb(p1,p8)*s56**(-1)*s78**(-1) - za(p2,p5)*za(p2,
     &    p8)*za(p4,p5)*zb(p1,p3)*zb(p1,p7)*iza(p2,p7)*iza(p5,p6)*izb(
     &    p1,p8)*s78**(-1) )
      A56(2,1,2,2,h56)= + s3456**(-1)*s178**(-1)*s356**(-1) * ( za(p1,
     &    p8)*za(p3,p5)*zb(p1,p7)*zb(p2,p3)*zb(p3,p6)*izb(p1,p8)*zab2(
     &    p4,p1,p8,p7)*s56**(-1)*s78**(-1) + za(p1,p8)*zb(p1,p7)*zb(p2,
     &    p6)*zb(p3,p6)*izb(p1,p8)*izb(p5,p6)*zab2(p4,p1,p8,p7)*
     &    s78**(-1) )
      A56(2,1,2,2,h56) = A56(2,1,2,2,h56) + s3456**(-1)*s178**(-1)*
     & s456**(-1) * (  - za(p1,p8)*za(p4,p5)*zb(p1,p7)*zb(p2,p3)*zb(p4,
     &    p6)*izb(p1,p8)*zab2(p4,p1,p8,p7)*s56**(-1)*s78**(-1) + za(p1,
     &    p8)*za(p4,p5)*zb(p1,p7)*zb(p2,p3)*iza(p5,p6)*izb(p1,p8)*zab2(
     &    p5,p1,p8,p7)*s78**(-1) )
      A56(2,1,2,2,h56) = A56(2,1,2,2,h56) + s3456**(-1)*s278**(-1)*
     & s356**(-1) * ( za(p1,p4)*za(p2,p8)*za(p3,p5)*zb(p2,p7)*zb(p3,p6)
     &    *iza(p2,p7)*zab2(p8,p2,p7,p3)*s56**(-1)*s78**(-1) + za(p1,p4)
     &    *za(p2,p8)*zb(p2,p7)*zb(p3,p6)*iza(p2,p7)*izb(p5,p6)*zab2(p8,
     &    p2,p7,p6)*s78**(-1) )
      A56(2,1,2,2,h56) = A56(2,1,2,2,h56) + s3456**(-1)*s278**(-1)*
     & s456**(-1) * (  - za(p1,p4)*za(p2,p8)*za(p4,p5)*zb(p2,p7)*zb(p4,
     &    p6)*iza(p2,p7)*zab2(p8,p2,p7,p3)*s56**(-1)*s78**(-1) + za(p1,
     &    p5)*za(p2,p8)*za(p4,p5)*zb(p2,p7)*iza(p2,p7)*iza(p5,p6)*zab2(
     &    p8,p2,p7,p3)*s78**(-1) )
      A56(2,1,2,2,h56) = A56(2,1,2,2,h56) + s3456**(-1)*s356**(-1) * (
     &     - za(p3,p5)*zb(p3,p6)*iza(p2,p7)*izb(p1,p8)*zab2(p4,p1,p8,p7
     &    )*zab2(p8,p2,p7,p3)*s56**(-1)*s78**(-1) - zb(p3,p6)*iza(p2,p7
     &    )*izb(p1,p8)*izb(p5,p6)*zab2(p4,p1,p8,p7)*zab2(p8,p2,p7,p6)*
     &    s78**(-1) )
      A56(2,1,2,2,h56) = A56(2,1,2,2,h56) + s3456**(-1)*s456**(-1) * (
     &    za(p4,p5)*zb(p4,p6)*iza(p2,p7)*izb(p1,p8)*zab2(p4,p1,p8,p7)*
     &    zab2(p8,p2,p7,p3)*s56**(-1)*s78**(-1) - za(p4,p5)*iza(p2,p7)*
     &    iza(p5,p6)*izb(p1,p8)*zab2(p5,p1,p8,p7)*zab2(p8,p2,p7,p3)*
     &    s78**(-1) )

      A56(2,2,1,1,h56)= + s3456**(-1)*s178**(-1)*s356**(-1) * (  - za(
     &    p2,p3)*za(p3,p5)*zb(p1,p7)*zb(p3,p6)*iza(p2,p8)*iza(p7,p8)*
     &    zab3(p2,p1,p7,p8,p4)*s56**(-1) - za(p2,p3)*za(p3,p5)*zb(p1,p8
     &    )*zb(p3,p6)*iza(p2,p7)*iza(p7,p8)*zab3(p2,p1,p7,p8,p4)*
     &    s56**(-1) + za(p2,p3)*za(p3,p5)*zb(p3,p6)*iza(p1,p8)*iza(p2,
     &    p7)*iza(p2,p8)*zab2(p2,p1,p8,p7)*zab3(p2,p1,p7,p8,p4)*
     &    s56**(-1) + za(p2,p5)*za(p3,p5)*zb(p1,p7)*iza(p2,p8)*iza(p5,
     &    p6)*iza(p7,p8)*zab3(p2,p1,p7,p8,p4) + za(p2,p5)*za(p3,p5)*zb(
     &    p1,p8)*iza(p2,p7)*iza(p5,p6)*iza(p7,p8)*zab3(p2,p1,p7,p8,p4)
     &     - za(p2,p5)*za(p3,p5)*iza(p1,p8)*iza(p2,p7)*iza(p2,p8)*iza(
     &    p5,p6)*zab2(p2,p1,p8,p7)*zab3(p2,p1,p7,p8,p4) )
      A56(2,2,1,1,h56) = A56(2,2,1,1,h56) + s3456**(-1)*s178**(-1)*
     & s456**(-1) * ( za(p2,p3)*za(p4,p5)*zb(p1,p7)*zb(p4,p6)*iza(p2,p8
     &    )*iza(p7,p8)*zab3(p2,p1,p7,p8,p4)*s56**(-1) + za(p2,p3)*za(p4
     &    ,p5)*zb(p1,p8)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*zab3(p2,p1,p7,
     &    p8,p4)*s56**(-1) - za(p2,p3)*za(p4,p5)*zb(p4,p6)*iza(p1,p8)*
     &    iza(p2,p7)*iza(p2,p8)*zab2(p2,p1,p8,p7)*zab3(p2,p1,p7,p8,p4)*
     &    s56**(-1) + za(p2,p3)*zb(p1,p7)*zb(p4,p6)*iza(p2,p8)*iza(p7,
     &    p8)*izb(p5,p6)*zab3(p2,p1,p7,p8,p6) + za(p2,p3)*zb(p1,p8)*zb(
     &    p4,p6)*iza(p2,p7)*iza(p7,p8)*izb(p5,p6)*zab3(p2,p1,p7,p8,p6)
     &     - za(p2,p3)*zb(p4,p6)*iza(p1,p8)*iza(p2,p7)*iza(p2,p8)*izb(
     &    p5,p6)*zab2(p2,p1,p8,p7)*zab3(p2,p1,p7,p8,p6) )
      A56(2,2,2,1,h56)= + s3456**(-1)*s278**(-1)*s356**(-1) * ( za(p1,
     &    p3)*za(p3,p5)*zb(p2,p7)*zb(p3,p6)*iza(p1,p8)*iza(p7,p8)*zab3(
     &    p1,p2,p7,p8,p4)*s56**(-1) + za(p1,p3)*za(p3,p5)*zb(p2,p8)*zb(
     &    p3,p6)*iza(p1,p7)*iza(p7,p8)*zab3(p1,p2,p7,p8,p4)*s56**(-1)
     &     + za(p1,p3)*za(p3,p5)*zb(p3,p6)*iza(p1,p7)*iza(p1,p8)*iza(p2
     &    ,p7)*zab2(p1,p2,p7,p8)*zab3(p1,p2,p7,p8,p4)*s56**(-1) - za(p1
     &    ,p5)*za(p3,p5)*zb(p2,p7)*iza(p1,p8)*iza(p5,p6)*iza(p7,p8)*
     &    zab3(p1,p2,p7,p8,p4) - za(p1,p5)*za(p3,p5)*zb(p2,p8)*iza(p1,
     &    p7)*iza(p5,p6)*iza(p7,p8)*zab3(p1,p2,p7,p8,p4) - za(p1,p5)*
     &    za(p3,p5)*iza(p1,p7)*iza(p1,p8)*iza(p2,p7)*iza(p5,p6)*zab2(p1
     &    ,p2,p7,p8)*zab3(p1,p2,p7,p8,p4) )
      A56(2,2,2,1,h56) = A56(2,2,2,1,h56) + s3456**(-1)*s278**(-1)*
     & s456**(-1) * (  - za(p1,p3)*za(p4,p5)*zb(p2,p7)*zb(p4,p6)*iza(p1
     &    ,p8)*iza(p7,p8)*zab3(p1,p2,p7,p8,p4)*s56**(-1) - za(p1,p3)*
     &    za(p4,p5)*zb(p2,p8)*zb(p4,p6)*iza(p1,p7)*iza(p7,p8)*zab3(p1,
     &    p2,p7,p8,p4)*s56**(-1) - za(p1,p3)*za(p4,p5)*zb(p4,p6)*iza(p1
     &    ,p7)*iza(p1,p8)*iza(p2,p7)*zab2(p1,p2,p7,p8)*zab3(p1,p2,p7,p8
     &    ,p4)*s56**(-1) - za(p1,p3)*zb(p2,p7)*zb(p4,p6)*iza(p1,p8)*
     &    iza(p7,p8)*izb(p5,p6)*zab3(p1,p2,p7,p8,p6) - za(p1,p3)*zb(p2,
     &    p8)*zb(p4,p6)*iza(p1,p7)*iza(p7,p8)*izb(p5,p6)*zab3(p1,p2,p7,
     &    p8,p6) - za(p1,p3)*zb(p4,p6)*iza(p1,p7)*iza(p1,p8)*iza(p2,p7)
     &    *izb(p5,p6)*zab2(p1,p2,p7,p8)*zab3(p1,p2,p7,p8,p6) )
      A56(2,2,1,2,h56)= + s3456**(-1)*s178**(-1)*s356**(-1) * (  - za(
     &    p2,p4)*za(p3,p5)*zb(p1,p7)*zb(p3,p6)*iza(p2,p8)*iza(p7,p8)*
     &    zab3(p2,p1,p7,p8,p3)*s56**(-1) - za(p2,p4)*za(p3,p5)*zb(p1,p8
     &    )*zb(p3,p6)*iza(p2,p7)*iza(p7,p8)*zab3(p2,p1,p7,p8,p3)*
     &    s56**(-1) + za(p2,p4)*za(p3,p5)*zb(p3,p6)*iza(p1,p8)*iza(p2,
     &    p7)*iza(p2,p8)*zab2(p2,p1,p8,p7)*zab3(p2,p1,p7,p8,p3)*
     &    s56**(-1) - za(p2,p4)*zb(p1,p7)*zb(p3,p6)*iza(p2,p8)*iza(p7,
     &    p8)*izb(p5,p6)*zab3(p2,p1,p7,p8,p6) - za(p2,p4)*zb(p1,p8)*zb(
     &    p3,p6)*iza(p2,p7)*iza(p7,p8)*izb(p5,p6)*zab3(p2,p1,p7,p8,p6)
     &     + za(p2,p4)*zb(p3,p6)*iza(p1,p8)*iza(p2,p7)*iza(p2,p8)*izb(
     &    p5,p6)*zab2(p2,p1,p8,p7)*zab3(p2,p1,p7,p8,p6) )
      A56(2,2,1,2,h56) = A56(2,2,1,2,h56) + s3456**(-1)*s178**(-1)*
     & s456**(-1) * ( za(p2,p4)*za(p4,p5)*zb(p1,p7)*zb(p4,p6)*iza(p2,p8
     &    )*iza(p7,p8)*zab3(p2,p1,p7,p8,p3)*s56**(-1) + za(p2,p4)*za(p4
     &    ,p5)*zb(p1,p8)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*zab3(p2,p1,p7,
     &    p8,p3)*s56**(-1) - za(p2,p4)*za(p4,p5)*zb(p4,p6)*iza(p1,p8)*
     &    iza(p2,p7)*iza(p2,p8)*zab2(p2,p1,p8,p7)*zab3(p2,p1,p7,p8,p3)*
     &    s56**(-1) - za(p2,p5)*za(p4,p5)*zb(p1,p7)*iza(p2,p8)*iza(p5,
     &    p6)*iza(p7,p8)*zab3(p2,p1,p7,p8,p3) - za(p2,p5)*za(p4,p5)*zb(
     &    p1,p8)*iza(p2,p7)*iza(p5,p6)*iza(p7,p8)*zab3(p2,p1,p7,p8,p3)
     &     + za(p2,p5)*za(p4,p5)*iza(p1,p8)*iza(p2,p7)*iza(p2,p8)*iza(
     &    p5,p6)*zab2(p2,p1,p8,p7)*zab3(p2,p1,p7,p8,p3) )
      A56(2,2,2,2,h56)= + s3456**(-1)*s278**(-1)*s356**(-1) * ( za(p1,
     &    p4)*za(p3,p5)*zb(p2,p7)*zb(p3,p6)*iza(p1,p8)*iza(p7,p8)*zab3(
     &    p1,p2,p7,p8,p3)*s56**(-1) + za(p1,p4)*za(p3,p5)*zb(p2,p8)*zb(
     &    p3,p6)*iza(p1,p7)*iza(p7,p8)*zab3(p1,p2,p7,p8,p3)*s56**(-1)
     &     + za(p1,p4)*za(p3,p5)*zb(p3,p6)*iza(p1,p7)*iza(p1,p8)*iza(p2
     &    ,p7)*zab2(p1,p2,p7,p8)*zab3(p1,p2,p7,p8,p3)*s56**(-1) + za(p1
     &    ,p4)*zb(p2,p7)*zb(p3,p6)*iza(p1,p8)*iza(p7,p8)*izb(p5,p6)*
     &    zab3(p1,p2,p7,p8,p6) + za(p1,p4)*zb(p2,p8)*zb(p3,p6)*iza(p1,
     &    p7)*iza(p7,p8)*izb(p5,p6)*zab3(p1,p2,p7,p8,p6) + za(p1,p4)*
     &    zb(p3,p6)*iza(p1,p7)*iza(p1,p8)*iza(p2,p7)*izb(p5,p6)*zab2(p1
     &    ,p2,p7,p8)*zab3(p1,p2,p7,p8,p6) )
      A56(2,2,2,2,h56) = A56(2,2,2,2,h56) + s3456**(-1)*s278**(-1)*
     & s456**(-1) * (  - za(p1,p4)*za(p4,p5)*zb(p2,p7)*zb(p4,p6)*iza(p1
     &    ,p8)*iza(p7,p8)*zab3(p1,p2,p7,p8,p3)*s56**(-1) - za(p1,p4)*
     &    za(p4,p5)*zb(p2,p8)*zb(p4,p6)*iza(p1,p7)*iza(p7,p8)*zab3(p1,
     &    p2,p7,p8,p3)*s56**(-1) - za(p1,p4)*za(p4,p5)*zb(p4,p6)*iza(p1
     &    ,p7)*iza(p1,p8)*iza(p2,p7)*zab2(p1,p2,p7,p8)*zab3(p1,p2,p7,p8
     &    ,p3)*s56**(-1) + za(p1,p5)*za(p4,p5)*zb(p2,p7)*iza(p1,p8)*
     &    iza(p5,p6)*iza(p7,p8)*zab3(p1,p2,p7,p8,p3) + za(p1,p5)*za(p4,
     &    p5)*zb(p2,p8)*iza(p1,p7)*iza(p5,p6)*iza(p7,p8)*zab3(p1,p2,p7,
     &    p8,p3) + za(p1,p5)*za(p4,p5)*iza(p1,p7)*iza(p1,p8)*iza(p2,p7)
     &    *iza(p5,p6)*zab2(p1,p2,p7,p8)*zab3(p1,p2,p7,p8,p3) )


      h56=2
      A56(1,1,1,1,h56)= + s3456**(-1)*s278**(-1)*s356**(-1) * ( za(p2,
     &    p7)*za(p3,p6)*zb(p1,p4)*zb(p3,p5)*izb(p1,p8)*izb(p7,p8)*zab3(
     &    p3,p2,p7,p8,p1)*s56**(-1) + za(p2,p7)*za(p3,p6)*zb(p1,p4)*
     &    iza(p5,p6)*izb(p1,p8)*izb(p7,p8)*zab3(p6,p2,p7,p8,p1) + za(p2
     &    ,p8)*za(p3,p6)*zb(p1,p4)*zb(p3,p5)*izb(p1,p7)*izb(p7,p8)*
     &    zab3(p3,p2,p7,p8,p1)*s56**(-1) + za(p2,p8)*za(p3,p6)*zb(p1,p4
     &    )*iza(p5,p6)*izb(p1,p7)*izb(p7,p8)*zab3(p6,p2,p7,p8,p1) + za(
     &    p3,p6)*zb(p1,p4)*zb(p3,p5)*izb(p1,p7)*izb(p1,p8)*izb(p2,p7)*
     &    zab2(p8,p2,p7,p1)*zab3(p3,p2,p7,p8,p1)*s56**(-1) + za(p3,p6)*
     &    zb(p1,p4)*iza(p5,p6)*izb(p1,p7)*izb(p1,p8)*izb(p2,p7)*zab2(p8
     &    ,p2,p7,p1)*zab3(p6,p2,p7,p8,p1) )
      A56(1,1,1,1,h56) = A56(1,1,1,1,h56) + s3456**(-1)*s278**(-1)*
     & s456**(-1) * (  - za(p2,p7)*za(p4,p6)*zb(p1,p4)*zb(p4,p5)*izb(p1
     &    ,p8)*izb(p7,p8)*zab3(p3,p2,p7,p8,p1)*s56**(-1) + za(p2,p7)*
     &    zb(p1,p5)*zb(p4,p5)*izb(p1,p8)*izb(p5,p6)*izb(p7,p8)*zab3(p3,
     &    p2,p7,p8,p1) - za(p2,p8)*za(p4,p6)*zb(p1,p4)*zb(p4,p5)*izb(p1
     &    ,p7)*izb(p7,p8)*zab3(p3,p2,p7,p8,p1)*s56**(-1) + za(p2,p8)*
     &    zb(p1,p5)*zb(p4,p5)*izb(p1,p7)*izb(p5,p6)*izb(p7,p8)*zab3(p3,
     &    p2,p7,p8,p1) - za(p4,p6)*zb(p1,p4)*zb(p4,p5)*izb(p1,p7)*izb(
     &    p1,p8)*izb(p2,p7)*zab2(p8,p2,p7,p1)*zab3(p3,p2,p7,p8,p1)*
     &    s56**(-1) + zb(p1,p5)*zb(p4,p5)*izb(p1,p7)*izb(p1,p8)*izb(p2,
     &    p7)*izb(p5,p6)*zab2(p8,p2,p7,p1)*zab3(p3,p2,p7,p8,p1) )
      A56(1,1,2,1,h56)= + s3456**(-1)*s178**(-1)*s356**(-1) * (  - za(
     &    p1,p7)*za(p3,p6)*zb(p2,p4)*zb(p3,p5)*izb(p1,p8)*izb(p7,p8)*
     &    zab3(p3,p1,p7,p8,p1)*s56**(-1) - za(p1,p7)*za(p3,p6)*zb(p2,p4
     &    )*iza(p5,p6)*izb(p1,p8)*izb(p7,p8)*zab3(p6,p1,p7,p8,p1) - za(
     &    p1,p8)*za(p3,p6)*zb(p2,p4)*zb(p3,p5)*izb(p1,p7)*izb(p7,p8)*
     &    zab3(p3,p1,p7,p8,p1)*s56**(-1) - za(p1,p8)*za(p3,p6)*zb(p2,p4
     &    )*iza(p5,p6)*izb(p1,p7)*izb(p7,p8)*zab3(p6,p1,p7,p8,p1) + za(
     &    p3,p6)*zb(p2,p4)*zb(p3,p5)*izb(p1,p7)*izb(p1,p8)**2*zab2(p7,
     &    p1,p8,p1)*zab3(p3,p1,p7,p8,p1)*s56**(-1) + za(p3,p6)*zb(p2,p4
     &    )*iza(p5,p6)*izb(p1,p7)*izb(p1,p8)**2*zab2(p7,p1,p8,p1)*zab3(
     &    p6,p1,p7,p8,p1) )
      A56(1,1,2,1,h56) = A56(1,1,2,1,h56) + s3456**(-1)*s178**(-1)*
     & s456**(-1) * ( za(p1,p7)*za(p4,p6)*zb(p2,p4)*zb(p4,p5)*izb(p1,p8
     &    )*izb(p7,p8)*zab3(p3,p1,p7,p8,p1)*s56**(-1) - za(p1,p7)*zb(p2
     &    ,p5)*zb(p4,p5)*izb(p1,p8)*izb(p5,p6)*izb(p7,p8)*zab3(p3,p1,p7
     &    ,p8,p1) + za(p1,p8)*za(p4,p6)*zb(p2,p4)*zb(p4,p5)*izb(p1,p7)*
     &    izb(p7,p8)*zab3(p3,p1,p7,p8,p1)*s56**(-1) - za(p1,p8)*zb(p2,
     &    p5)*zb(p4,p5)*izb(p1,p7)*izb(p5,p6)*izb(p7,p8)*zab3(p3,p1,p7,
     &    p8,p1) - za(p4,p6)*zb(p2,p4)*zb(p4,p5)*izb(p1,p7)*izb(p1,p8)
     &    **2*zab2(p7,p1,p8,p1)*zab3(p3,p1,p7,p8,p1)*s56**(-1) + zb(p2,
     &    p5)*zb(p4,p5)*izb(p1,p7)*izb(p1,p8)**2*izb(p5,p6)*zab2(p7,p1,
     &    p8,p1)*zab3(p3,p1,p7,p8,p1) )
      A56(1,1,2,1,h56) = A56(1,1,2,1,h56) + s3456**(-1)*s278**(-1)*
     & s356**(-1) * (  - za(p1,p3)*za(p3,p6)*zb(p1,p2)**2*zb(p3,p5)*
     &    izb(p1,p7)*izb(p1,p8)*izb(p2,p7)*zab2(p8,p2,p7,p4)*s56**(-1)
     &     - za(p1,p3)*za(p3,p6)*zb(p1,p2)*zb(p3,p5)*izb(p1,p7)*izb(p7,
     &    p8)*zab2(p8,p2,p7,p4)*s56**(-1) - za(p1,p3)*za(p3,p6)*zb(p1,
     &    p2)*zb(p3,p5)*izb(p1,p8)*izb(p7,p8)*zab2(p7,p2,p8,p4)*
     &    s56**(-1) - za(p1,p6)*za(p3,p6)*zb(p1,p2)**2*iza(p5,p6)*izb(
     &    p1,p7)*izb(p1,p8)*izb(p2,p7)*zab2(p8,p2,p7,p4) - za(p1,p6)*
     &    za(p3,p6)*zb(p1,p2)*iza(p5,p6)*izb(p1,p7)*izb(p7,p8)*zab2(p8,
     &    p2,p7,p4) - za(p1,p6)*za(p3,p6)*zb(p1,p2)*iza(p5,p6)*izb(p1,
     &    p8)*izb(p7,p8)*zab2(p7,p2,p8,p4) )
      A56(1,1,2,1,h56) = A56(1,1,2,1,h56) + s3456**(-1)*s278**(-1)*
     & s456**(-1) * ( za(p1,p3)*za(p4,p6)*zb(p1,p2)**2*zb(p4,p5)*izb(p1
     &    ,p7)*izb(p1,p8)*izb(p2,p7)*zab2(p8,p2,p7,p4)*s56**(-1) + za(
     &    p1,p3)*za(p4,p6)*zb(p1,p2)*zb(p4,p5)*izb(p1,p7)*izb(p7,p8)*
     &    zab2(p8,p2,p7,p4)*s56**(-1) + za(p1,p3)*za(p4,p6)*zb(p1,p2)*
     &    zb(p4,p5)*izb(p1,p8)*izb(p7,p8)*zab2(p7,p2,p8,p4)*s56**(-1)
     &     - za(p1,p3)*zb(p1,p2)**2*zb(p4,p5)*izb(p1,p7)*izb(p1,p8)*
     &    izb(p2,p7)*izb(p5,p6)*zab2(p8,p2,p7,p5) - za(p1,p3)*zb(p1,p2)
     &    *zb(p4,p5)*izb(p1,p7)*izb(p5,p6)*izb(p7,p8)*zab2(p8,p2,p7,p5)
     &     - za(p1,p3)*zb(p1,p2)*zb(p4,p5)*izb(p1,p8)*izb(p5,p6)*izb(p7
     &    ,p8)*zab2(p7,p2,p8,p5) )
      A56(1,1,2,1,h56) = A56(1,1,2,1,h56) + s3456**(-1)*s356**(-1) * (
     &    za(p3,p6)*zb(p1,p2)*zb(p2,p4)*zb(p3,p5)*izb(p1,p7)*izb(p1,p8)
     &    **2*izb(p2,p7)*zab2(p3,p1,p8,p1)*s56**(-1) + za(p3,p6)*zb(p1,
     &    p2)*zb(p2,p4)*iza(p5,p6)*izb(p1,p7)*izb(p1,p8)**2*izb(p2,p7)*
     &    zab2(p6,p1,p8,p1) )
      A56(1,1,2,1,h56) = A56(1,1,2,1,h56) + s3456**(-1)*s456**(-1) * (
     &     - za(p4,p6)*zb(p1,p2)*zb(p2,p4)*zb(p4,p5)*izb(p1,p7)*izb(p1,
     &    p8)**2*izb(p2,p7)*zab2(p3,p1,p8,p1)*s56**(-1) + zb(p1,p2)*zb(
     &    p2,p5)*zb(p4,p5)*izb(p1,p7)*izb(p1,p8)**2*izb(p2,p7)*izb(p5,
     &    p6)*zab2(p3,p1,p8,p1) )
      A56(1,1,1,2,h56)= + s3456**(-1)*s278**(-1)*s356**(-1) * ( za(p2,
     &    p7)*za(p3,p6)*zb(p1,p3)*zb(p3,p5)*izb(p1,p8)*izb(p7,p8)*zab3(
     &    p4,p2,p7,p8,p1)*s56**(-1) - za(p2,p7)*zb(p1,p5)*zb(p3,p5)*
     &    izb(p1,p8)*izb(p5,p6)*izb(p7,p8)*zab3(p4,p2,p7,p8,p1) + za(p2
     &    ,p8)*za(p3,p6)*zb(p1,p3)*zb(p3,p5)*izb(p1,p7)*izb(p7,p8)*
     &    zab3(p4,p2,p7,p8,p1)*s56**(-1) - za(p2,p8)*zb(p1,p5)*zb(p3,p5
     &    )*izb(p1,p7)*izb(p5,p6)*izb(p7,p8)*zab3(p4,p2,p7,p8,p1) + za(
     &    p3,p6)*zb(p1,p3)*zb(p3,p5)*izb(p1,p7)*izb(p1,p8)*izb(p2,p7)*
     &    zab2(p8,p2,p7,p1)*zab3(p4,p2,p7,p8,p1)*s56**(-1) - zb(p1,p5)*
     &    zb(p3,p5)*izb(p1,p7)*izb(p1,p8)*izb(p2,p7)*izb(p5,p6)*zab2(p8
     &    ,p2,p7,p1)*zab3(p4,p2,p7,p8,p1) )
      A56(1,1,1,2,h56) = A56(1,1,1,2,h56) + s3456**(-1)*s278**(-1)*
     & s456**(-1) * (  - za(p2,p7)*za(p4,p6)*zb(p1,p3)*zb(p4,p5)*izb(p1
     &    ,p8)*izb(p7,p8)*zab3(p4,p2,p7,p8,p1)*s56**(-1) - za(p2,p7)*
     &    za(p4,p6)*zb(p1,p3)*iza(p5,p6)*izb(p1,p8)*izb(p7,p8)*zab3(p6,
     &    p2,p7,p8,p1) - za(p2,p8)*za(p4,p6)*zb(p1,p3)*zb(p4,p5)*izb(p1
     &    ,p7)*izb(p7,p8)*zab3(p4,p2,p7,p8,p1)*s56**(-1) - za(p2,p8)*
     &    za(p4,p6)*zb(p1,p3)*iza(p5,p6)*izb(p1,p7)*izb(p7,p8)*zab3(p6,
     &    p2,p7,p8,p1) - za(p4,p6)*zb(p1,p3)*zb(p4,p5)*izb(p1,p7)*izb(
     &    p1,p8)*izb(p2,p7)*zab2(p8,p2,p7,p1)*zab3(p4,p2,p7,p8,p1)*
     &    s56**(-1) - za(p4,p6)*zb(p1,p3)*iza(p5,p6)*izb(p1,p7)*izb(p1,
     &    p8)*izb(p2,p7)*zab2(p8,p2,p7,p1)*zab3(p6,p2,p7,p8,p1) )
      A56(1,1,2,2,h56)= + s3456**(-1)*s178**(-1)*s356**(-1) * (  - za(
     &    p1,p7)*za(p3,p6)*zb(p2,p3)*zb(p3,p5)*izb(p2,p8)*izb(p7,p8)*
     &    zab3(p4,p1,p7,p8,p2)*s56**(-1) + za(p1,p7)*zb(p2,p5)*zb(p3,p5
     &    )*izb(p2,p8)*izb(p5,p6)*izb(p7,p8)*zab3(p4,p1,p7,p8,p2) - za(
     &    p1,p8)*za(p3,p6)*zb(p2,p3)*zb(p3,p5)*izb(p2,p7)*izb(p7,p8)*
     &    zab3(p4,p1,p7,p8,p2)*s56**(-1) + za(p1,p8)*zb(p2,p5)*zb(p3,p5
     &    )*izb(p2,p7)*izb(p5,p6)*izb(p7,p8)*zab3(p4,p1,p7,p8,p2) + za(
     &    p3,p6)*zb(p2,p3)*zb(p3,p5)*izb(p1,p8)*izb(p2,p7)*izb(p2,p8)*
     &    zab2(p7,p1,p8,p2)*zab3(p4,p1,p7,p8,p2)*s56**(-1) - zb(p2,p5)*
     &    zb(p3,p5)*izb(p1,p8)*izb(p2,p7)*izb(p2,p8)*izb(p5,p6)*zab2(p7
     &    ,p1,p8,p2)*zab3(p4,p1,p7,p8,p2) )
      A56(1,1,2,2,h56) = A56(1,1,2,2,h56) + s3456**(-1)*s178**(-1)*
     & s456**(-1) * ( za(p1,p7)*za(p4,p6)*zb(p2,p3)*zb(p4,p5)*izb(p2,p8
     &    )*izb(p7,p8)*zab3(p4,p1,p7,p8,p2)*s56**(-1) + za(p1,p7)*za(p4
     &    ,p6)*zb(p2,p3)*iza(p5,p6)*izb(p2,p8)*izb(p7,p8)*zab3(p6,p1,p7
     &    ,p8,p2) + za(p1,p8)*za(p4,p6)*zb(p2,p3)*zb(p4,p5)*izb(p2,p7)*
     &    izb(p7,p8)*zab3(p4,p1,p7,p8,p2)*s56**(-1) + za(p1,p8)*za(p4,
     &    p6)*zb(p2,p3)*iza(p5,p6)*izb(p2,p7)*izb(p7,p8)*zab3(p6,p1,p7,
     &    p8,p2) - za(p4,p6)*zb(p2,p3)*zb(p4,p5)*izb(p1,p8)*izb(p2,p7)*
     &    izb(p2,p8)*zab2(p7,p1,p8,p2)*zab3(p4,p1,p7,p8,p2)*s56**(-1)
     &     - za(p4,p6)*zb(p2,p3)*iza(p5,p6)*izb(p1,p8)*izb(p2,p7)*izb(
     &    p2,p8)*zab2(p7,p1,p8,p2)*zab3(p6,p1,p7,p8,p2) )

      A56(1,2,1,1,h56)= + s3456**(-1)*s178**(-1)*s356**(-1) * ( za(p1,
     &    p7)*za(p2,p3)*za(p3,p6)*zb(p1,p8)*zb(p3,p5)*iza(p1,p8)*zab2(
     &    p7,p1,p8,p4)*s56**(-1)*s78**(-1) + za(p1,p7)*za(p2,p6)*za(p3,
     &    p6)*zb(p1,p8)*iza(p1,p8)*iza(p5,p6)*zab2(p7,p1,p8,p4)*
     &    s78**(-1) )
      A56(1,2,1,1,h56) = A56(1,2,1,1,h56) + s3456**(-1)*s178**(-1)*
     & s456**(-1) * (  - za(p1,p7)*za(p2,p3)*za(p4,p6)*zb(p1,p8)*zb(p4,
     &    p5)*iza(p1,p8)*zab2(p7,p1,p8,p4)*s56**(-1)*s78**(-1) + za(p1,
     &    p7)*za(p2,p3)*zb(p1,p8)*zb(p4,p5)*iza(p1,p8)*izb(p5,p6)*zab2(
     &    p7,p1,p8,p5)*s78**(-1) )
      A56(1,2,1,1,h56) = A56(1,2,1,1,h56) + s3456**(-1)*s278**(-1)*
     & s356**(-1) * ( za(p2,p7)*za(p3,p6)*zb(p1,p4)*zb(p2,p8)*zb(p3,p5)
     &    *izb(p2,p7)*zab2(p3,p2,p7,p8)*s56**(-1)*s78**(-1) + za(p2,p7)
     &    *za(p3,p6)*zb(p1,p4)*zb(p2,p8)*iza(p5,p6)*izb(p2,p7)*zab2(p6,
     &    p2,p7,p8)*s78**(-1) )
      A56(1,2,1,1,h56) = A56(1,2,1,1,h56) + s3456**(-1)*s278**(-1)*
     & s456**(-1) * (  - za(p2,p7)*za(p4,p6)*zb(p1,p4)*zb(p2,p8)*zb(p4,
     &    p5)*izb(p2,p7)*zab2(p3,p2,p7,p8)*s56**(-1)*s78**(-1) + za(p2,
     &    p7)*zb(p1,p5)*zb(p2,p8)*zb(p4,p5)*izb(p2,p7)*izb(p5,p6)*zab2(
     &    p3,p2,p7,p8)*s78**(-1) )
      A56(1,2,1,1,h56) = A56(1,2,1,1,h56) + s3456**(-1)*s356**(-1) * (
     &     - za(p3,p6)*zb(p3,p5)*iza(p1,p8)*izb(p2,p7)*zab2(p3,p2,p7,p8
     &    )*zab2(p7,p1,p8,p4)*s56**(-1)*s78**(-1) - za(p3,p6)*iza(p1,p8
     &    )*iza(p5,p6)*izb(p2,p7)*zab2(p6,p2,p7,p8)*zab2(p7,p1,p8,p4)*
     &    s78**(-1) )
      A56(1,2,1,1,h56) = A56(1,2,1,1,h56) + s3456**(-1)*s456**(-1) * (
     &    za(p4,p6)*zb(p4,p5)*iza(p1,p8)*izb(p2,p7)*zab2(p3,p2,p7,p8)*
     &    zab2(p7,p1,p8,p4)*s56**(-1)*s78**(-1) - zb(p4,p5)*iza(p1,p8)*
     &    izb(p2,p7)*izb(p5,p6)*zab2(p3,p2,p7,p8)*zab2(p7,p1,p8,p5)*
     &    s78**(-1) )
      A56(1,2,2,1,h56)= + s3456**(-1)*s178**(-1)*s356**(-1) * ( za(p1,
     &    p7)**2*za(p3,p6)*zb(p2,p4)*zb(p3,p5)*iza(p1,p8)*zab2(p3,p1,p7
     &    ,p8)*s56**(-1)*s78**(-1) + za(p1,p7)**2*za(p3,p6)*zb(p2,p4)*
     &    iza(p1,p8)*iza(p5,p6)*zab2(p6,p1,p7,p8)*s78**(-1) )
      A56(1,2,2,1,h56) = A56(1,2,2,1,h56) + s3456**(-1)*s178**(-1)*
     & s456**(-1) * (  - za(p1,p7)**2*za(p4,p6)*zb(p2,p4)*zb(p4,p5)*
     &    iza(p1,p8)*zab2(p3,p1,p7,p8)*s56**(-1)*s78**(-1) + za(p1,p7)
     &    **2*zb(p2,p5)*zb(p4,p5)*iza(p1,p8)*izb(p5,p6)*zab2(p3,p1,p7,
     &    p8)*s78**(-1) )
      A56(1,2,2,1,h56) = A56(1,2,2,1,h56) + s3456**(-1)*s278**(-1)*
     & s356**(-1) * ( za(p1,p3)*za(p3,p6)*zb(p2,p8)**2*zb(p3,p5)*izb(p2
     &    ,p7)*zab2(p7,p2,p8,p4)*s56**(-1)*s78**(-1) + za(p1,p6)*za(p3,
     &    p6)*zb(p2,p8)**2*iza(p5,p6)*izb(p2,p7)*zab2(p7,p2,p8,p4)*
     &    s78**(-1) )
      A56(1,2,2,1,h56) = A56(1,2,2,1,h56) + s3456**(-1)*s278**(-1)*
     & s456**(-1) * (  - za(p1,p3)*za(p4,p6)*zb(p2,p8)**2*zb(p4,p5)*
     &    izb(p2,p7)*zab2(p7,p2,p8,p4)*s56**(-1)*s78**(-1) + za(p1,p3)*
     &    zb(p2,p8)**2*zb(p4,p5)*izb(p2,p7)*izb(p5,p6)*zab2(p7,p2,p8,p5
     &    )*s78**(-1) )
      A56(1,2,2,1,h56) = A56(1,2,2,1,h56) + s3456**(-1)*s356**(-1) * (
     &     - za(p1,p3)*za(p1,p7)*za(p3,p6)*zb(p2,p4)*zb(p2,p8)*zb(p3,p5
     &    )*iza(p1,p8)*izb(p2,p7)*s56**(-1)*s78**(-1) - za(p1,p6)*za(p1
     &    ,p7)*za(p3,p6)*zb(p2,p4)*zb(p2,p8)*iza(p1,p8)*iza(p5,p6)*izb(
     &    p2,p7)*s78**(-1) )
      A56(1,2,2,1,h56) = A56(1,2,2,1,h56) + s3456**(-1)*s456**(-1) * (
     &    za(p1,p3)*za(p1,p7)*za(p4,p6)*zb(p2,p4)*zb(p2,p8)*zb(p4,p5)*
     &    iza(p1,p8)*izb(p2,p7)*s56**(-1)*s78**(-1) - za(p1,p3)*za(p1,
     &    p7)*zb(p2,p5)*zb(p2,p8)*zb(p4,p5)*iza(p1,p8)*izb(p2,p7)*izb(
     &    p5,p6)*s78**(-1) )
      A56(1,2,1,2,h56)= + s3456**(-1)*s178**(-1)*s356**(-1) * ( za(p1,
     &    p7)*za(p2,p4)*za(p3,p6)*zb(p1,p8)*zb(p3,p5)*iza(p1,p8)*zab2(
     &    p7,p1,p8,p3)*s56**(-1)*s78**(-1) - za(p1,p7)*za(p2,p4)*zb(p1,
     &    p8)*zb(p3,p5)*iza(p1,p8)*izb(p5,p6)*zab2(p7,p1,p8,p5)*
     &    s78**(-1) )
      A56(1,2,1,2,h56) = A56(1,2,1,2,h56) + s3456**(-1)*s178**(-1)*
     & s456**(-1) * (  - za(p1,p7)*za(p2,p4)*za(p4,p6)*zb(p1,p8)*zb(p4,
     &    p5)*iza(p1,p8)*zab2(p7,p1,p8,p3)*s56**(-1)*s78**(-1) - za(p1,
     &    p7)*za(p2,p6)*za(p4,p6)*zb(p1,p8)*iza(p1,p8)*iza(p5,p6)*zab2(
     &    p7,p1,p8,p3)*s78**(-1) )
      A56(1,2,1,2,h56) = A56(1,2,1,2,h56) + s3456**(-1)*s278**(-1)*
     & s356**(-1) * ( za(p2,p7)*za(p3,p6)*zb(p1,p3)*zb(p2,p8)*zb(p3,p5)
     &    *izb(p2,p7)*zab2(p4,p2,p7,p8)*s56**(-1)*s78**(-1) - za(p2,p7)
     &    *zb(p1,p5)*zb(p2,p8)*zb(p3,p5)*izb(p2,p7)*izb(p5,p6)*zab2(p4,
     &    p2,p7,p8)*s78**(-1) )
      A56(1,2,1,2,h56) = A56(1,2,1,2,h56) + s3456**(-1)*s278**(-1)*
     & s456**(-1) * (  - za(p2,p7)*za(p4,p6)*zb(p1,p3)*zb(p2,p8)*zb(p4,
     &    p5)*izb(p2,p7)*zab2(p4,p2,p7,p8)*s56**(-1)*s78**(-1) - za(p2,
     &    p7)*za(p4,p6)*zb(p1,p3)*zb(p2,p8)*iza(p5,p6)*izb(p2,p7)*zab2(
     &    p6,p2,p7,p8)*s78**(-1) )
      A56(1,2,1,2,h56) = A56(1,2,1,2,h56) + s3456**(-1)*s356**(-1) * (
     &     - za(p3,p6)*zb(p3,p5)*iza(p1,p8)*izb(p2,p7)*zab2(p4,p2,p7,p8
     &    )*zab2(p7,p1,p8,p3)*s56**(-1)*s78**(-1) + zb(p3,p5)*iza(p1,p8
     &    )*izb(p2,p7)*izb(p5,p6)*zab2(p4,p2,p7,p8)*zab2(p7,p1,p8,p5)*
     &    s78**(-1) )
      A56(1,2,1,2,h56) = A56(1,2,1,2,h56) + s3456**(-1)*s456**(-1) * (
     &    za(p4,p6)*zb(p4,p5)*iza(p1,p8)*izb(p2,p7)*zab2(p4,p2,p7,p8)*
     &    zab2(p7,p1,p8,p3)*s56**(-1)*s78**(-1) + za(p4,p6)*iza(p1,p8)*
     &    iza(p5,p6)*izb(p2,p7)*zab2(p6,p2,p7,p8)*zab2(p7,p1,p8,p3)*
     &    s78**(-1) )
      A56(1,2,2,2,h56)= + s3456**(-1)*s178**(-1)*s356**(-1) * ( za(p1,
     &    p7)**2*za(p3,p6)*zb(p2,p3)*zb(p3,p5)*iza(p1,p8)*zab2(p4,p1,p7
     &    ,p8)*s56**(-1)*s78**(-1) - za(p1,p7)**2*zb(p2,p5)*zb(p3,p5)*
     &    iza(p1,p8)*izb(p5,p6)*zab2(p4,p1,p7,p8)*s78**(-1) )
      A56(1,2,2,2,h56) = A56(1,2,2,2,h56) + s3456**(-1)*s178**(-1)*
     & s456**(-1) * (  - za(p1,p7)**2*za(p4,p6)*zb(p2,p3)*zb(p4,p5)*
     &    iza(p1,p8)*zab2(p4,p1,p7,p8)*s56**(-1)*s78**(-1) - za(p1,p7)
     &    **2*za(p4,p6)*zb(p2,p3)*iza(p1,p8)*iza(p5,p6)*zab2(p6,p1,p7,
     &    p8)*s78**(-1) )
      A56(1,2,2,2,h56) = A56(1,2,2,2,h56) + s3456**(-1)*s278**(-1)*
     & s356**(-1) * ( za(p1,p4)*za(p3,p6)*zb(p2,p8)**2*zb(p3,p5)*izb(p2
     &    ,p7)*zab2(p7,p2,p8,p3)*s56**(-1)*s78**(-1) - za(p1,p4)*zb(p2,
     &    p8)**2*zb(p3,p5)*izb(p2,p7)*izb(p5,p6)*zab2(p7,p2,p8,p5)*
     &    s78**(-1) )
      A56(1,2,2,2,h56) = A56(1,2,2,2,h56) + s3456**(-1)*s278**(-1)*
     & s456**(-1) * (  - za(p1,p4)*za(p4,p6)*zb(p2,p8)**2*zb(p4,p5)*
     &    izb(p2,p7)*zab2(p7,p2,p8,p3)*s56**(-1)*s78**(-1) - za(p1,p6)*
     &    za(p4,p6)*zb(p2,p8)**2*iza(p5,p6)*izb(p2,p7)*zab2(p7,p2,p8,p3
     &    )*s78**(-1) )
      A56(1,2,2,2,h56) = A56(1,2,2,2,h56) + s3456**(-1)*s356**(-1) * (
     &     - za(p1,p4)*za(p1,p7)*za(p3,p6)*zb(p2,p3)*zb(p2,p8)*zb(p3,p5
     &    )*iza(p1,p8)*izb(p2,p7)*s56**(-1)*s78**(-1) + za(p1,p4)*za(p1
     &    ,p7)*zb(p2,p5)*zb(p2,p8)*zb(p3,p5)*iza(p1,p8)*izb(p2,p7)*izb(
     &    p5,p6)*s78**(-1) )
      A56(1,2,2,2,h56) = A56(1,2,2,2,h56) + s3456**(-1)*s456**(-1) * (
     &    za(p1,p4)*za(p1,p7)*za(p4,p6)*zb(p2,p3)*zb(p2,p8)*zb(p4,p5)*
     &    iza(p1,p8)*izb(p2,p7)*s56**(-1)*s78**(-1) + za(p1,p6)*za(p1,
     &    p7)*za(p4,p6)*zb(p2,p3)*zb(p2,p8)*iza(p1,p8)*iza(p5,p6)*izb(
     &    p2,p7)*s78**(-1) )

      A56(2,1,1,1,h56)= + s3456**(-1)*s178**(-1)*s356**(-1) * ( za(p2,
     &    p3)*za(p3,p6)*zb(p1,p7)**2*zb(p3,p5)*izb(p1,p8)*zab2(p8,p1,p7
     &    ,p4)*s56**(-1)*s78**(-1) + za(p2,p6)*za(p3,p6)*zb(p1,p7)**2*
     &    iza(p5,p6)*izb(p1,p8)*zab2(p8,p1,p7,p4)*s78**(-1) )
      A56(2,1,1,1,h56) = A56(2,1,1,1,h56) + s3456**(-1)*s178**(-1)*
     & s456**(-1) * (  - za(p2,p3)*za(p4,p6)*zb(p1,p7)**2*zb(p4,p5)*
     &    izb(p1,p8)*zab2(p8,p1,p7,p4)*s56**(-1)*s78**(-1) + za(p2,p3)*
     &    zb(p1,p7)**2*zb(p4,p5)*izb(p1,p8)*izb(p5,p6)*zab2(p8,p1,p7,p5
     &    )*s78**(-1) )
      A56(2,1,1,1,h56) = A56(2,1,1,1,h56) + s3456**(-1)*s278**(-1)*
     & s356**(-1) * ( za(p2,p8)**2*za(p3,p6)*zb(p1,p4)*zb(p3,p5)*iza(p2
     &    ,p7)*zab2(p3,p2,p8,p7)*s56**(-1)*s78**(-1) + za(p2,p8)**2*za(
     &    p3,p6)*zb(p1,p4)*iza(p2,p7)*iza(p5,p6)*zab2(p6,p2,p8,p7)*
     &    s78**(-1) )
      A56(2,1,1,1,h56) = A56(2,1,1,1,h56) + s3456**(-1)*s278**(-1)*
     & s456**(-1) * (  - za(p2,p8)**2*za(p4,p6)*zb(p1,p4)*zb(p4,p5)*
     &    iza(p2,p7)*zab2(p3,p2,p8,p7)*s56**(-1)*s78**(-1) + za(p2,p8)
     &    **2*zb(p1,p5)*zb(p4,p5)*iza(p2,p7)*izb(p5,p6)*zab2(p3,p2,p8,
     &    p7)*s78**(-1) )
      A56(2,1,1,1,h56) = A56(2,1,1,1,h56) + s3456**(-1)*s356**(-1) * (
     &     - za(p2,p3)*za(p2,p8)*za(p3,p6)*zb(p1,p4)*zb(p1,p7)*zb(p3,p5
     &    )*iza(p2,p7)*izb(p1,p8)*s56**(-1)*s78**(-1) - za(p2,p6)*za(p2
     &    ,p8)*za(p3,p6)*zb(p1,p4)*zb(p1,p7)*iza(p2,p7)*iza(p5,p6)*izb(
     &    p1,p8)*s78**(-1) )
      A56(2,1,1,1,h56) = A56(2,1,1,1,h56) + s3456**(-1)*s456**(-1) * (
     &    za(p2,p3)*za(p2,p8)*za(p4,p6)*zb(p1,p4)*zb(p1,p7)*zb(p4,p5)*
     &    iza(p2,p7)*izb(p1,p8)*s56**(-1)*s78**(-1) - za(p2,p3)*za(p2,
     &    p8)*zb(p1,p5)*zb(p1,p7)*zb(p4,p5)*iza(p2,p7)*izb(p1,p8)*izb(
     &    p5,p6)*s78**(-1) )
      A56(2,1,2,1,h56)= + s3456**(-1)*s178**(-1)*s356**(-1) * ( za(p1,
     &    p8)*za(p3,p6)*zb(p1,p7)*zb(p2,p4)*zb(p3,p5)*izb(p1,p8)*zab2(
     &    p3,p1,p8,p7)*s56**(-1)*s78**(-1) + za(p1,p8)*za(p3,p6)*zb(p1,
     &    p7)*zb(p2,p4)*iza(p5,p6)*izb(p1,p8)*zab2(p6,p1,p8,p7)*
     &    s78**(-1) )
      A56(2,1,2,1,h56) = A56(2,1,2,1,h56) + s3456**(-1)*s178**(-1)*
     & s456**(-1) * (  - za(p1,p8)*za(p4,p6)*zb(p1,p7)*zb(p2,p4)*zb(p4,
     &    p5)*izb(p1,p8)*zab2(p3,p1,p8,p7)*s56**(-1)*s78**(-1) + za(p1,
     &    p8)*zb(p1,p7)*zb(p2,p5)*zb(p4,p5)*izb(p1,p8)*izb(p5,p6)*zab2(
     &    p3,p1,p8,p7)*s78**(-1) )
      A56(2,1,2,1,h56) = A56(2,1,2,1,h56) + s3456**(-1)*s278**(-1)*
     & s356**(-1) * ( za(p1,p3)*za(p2,p8)*za(p3,p6)*zb(p2,p7)*zb(p3,p5)
     &    *iza(p2,p7)*zab2(p8,p2,p7,p4)*s56**(-1)*s78**(-1) + za(p1,p6)
     &    *za(p2,p8)*za(p3,p6)*zb(p2,p7)*iza(p2,p7)*iza(p5,p6)*zab2(p8,
     &    p2,p7,p4)*s78**(-1) )
      A56(2,1,2,1,h56) = A56(2,1,2,1,h56) + s3456**(-1)*s278**(-1)*
     & s456**(-1) * (  - za(p1,p3)*za(p2,p8)*za(p4,p6)*zb(p2,p7)*zb(p4,
     &    p5)*iza(p2,p7)*zab2(p8,p2,p7,p4)*s56**(-1)*s78**(-1) + za(p1,
     &    p3)*za(p2,p8)*zb(p2,p7)*zb(p4,p5)*iza(p2,p7)*izb(p5,p6)*zab2(
     &    p8,p2,p7,p5)*s78**(-1) )
      A56(2,1,2,1,h56) = A56(2,1,2,1,h56) + s3456**(-1)*s356**(-1) * (
     &     - za(p3,p6)*zb(p3,p5)*iza(p2,p7)*izb(p1,p8)*zab2(p3,p1,p8,p7
     &    )*zab2(p8,p2,p7,p4)*s56**(-1)*s78**(-1) - za(p3,p6)*iza(p2,p7
     &    )*iza(p5,p6)*izb(p1,p8)*zab2(p6,p1,p8,p7)*zab2(p8,p2,p7,p4)*
     &    s78**(-1) )
      A56(2,1,2,1,h56) = A56(2,1,2,1,h56) + s3456**(-1)*s456**(-1) * (
     &    za(p4,p6)*zb(p4,p5)*iza(p2,p7)*izb(p1,p8)*zab2(p3,p1,p8,p7)*
     &    zab2(p8,p2,p7,p4)*s56**(-1)*s78**(-1) - zb(p4,p5)*iza(p2,p7)*
     &    izb(p1,p8)*izb(p5,p6)*zab2(p3,p1,p8,p7)*zab2(p8,p2,p7,p5)*
     &    s78**(-1) )
      A56(2,1,1,2,h56)= + s3456**(-1)*s178**(-1)*s356**(-1) * ( za(p2,
     &    p4)*za(p3,p6)*zb(p1,p7)**2*zb(p3,p5)*izb(p1,p8)*zab2(p8,p1,p7
     &    ,p3)*s56**(-1)*s78**(-1) - za(p2,p4)*zb(p1,p7)**2*zb(p3,p5)*
     &    izb(p1,p8)*izb(p5,p6)*zab2(p8,p1,p7,p5)*s78**(-1) )
      A56(2,1,1,2,h56) = A56(2,1,1,2,h56) + s3456**(-1)*s178**(-1)*
     & s456**(-1) * (  - za(p2,p4)*za(p4,p6)*zb(p1,p7)**2*zb(p4,p5)*
     &    izb(p1,p8)*zab2(p8,p1,p7,p3)*s56**(-1)*s78**(-1) - za(p2,p6)*
     &    za(p4,p6)*zb(p1,p7)**2*iza(p5,p6)*izb(p1,p8)*zab2(p8,p1,p7,p3
     &    )*s78**(-1) )
      A56(2,1,1,2,h56) = A56(2,1,1,2,h56) + s3456**(-1)*s278**(-1)*
     & s356**(-1) * ( za(p2,p8)**2*za(p3,p6)*zb(p1,p3)*zb(p3,p5)*iza(p2
     &    ,p7)*zab2(p4,p2,p8,p7)*s56**(-1)*s78**(-1) - za(p2,p8)**2*zb(
     &    p1,p5)*zb(p3,p5)*iza(p2,p7)*izb(p5,p6)*zab2(p4,p2,p8,p7)*
     &    s78**(-1) )
      A56(2,1,1,2,h56) = A56(2,1,1,2,h56) + s3456**(-1)*s278**(-1)*
     & s456**(-1) * (  - za(p2,p8)**2*za(p4,p6)*zb(p1,p3)*zb(p4,p5)*
     &    iza(p2,p7)*zab2(p4,p2,p8,p7)*s56**(-1)*s78**(-1) - za(p2,p8)
     &    **2*za(p4,p6)*zb(p1,p3)*iza(p2,p7)*iza(p5,p6)*zab2(p6,p2,p8,
     &    p7)*s78**(-1) )
      A56(2,1,1,2,h56) = A56(2,1,1,2,h56) + s3456**(-1)*s356**(-1) * (
     &     - za(p2,p4)*za(p2,p8)*za(p3,p6)*zb(p1,p3)*zb(p1,p7)*zb(p3,p5
     &    )*iza(p2,p7)*izb(p1,p8)*s56**(-1)*s78**(-1) + za(p2,p4)*za(p2
     &    ,p8)*zb(p1,p5)*zb(p1,p7)*zb(p3,p5)*iza(p2,p7)*izb(p1,p8)*izb(
     &    p5,p6)*s78**(-1) )
      A56(2,1,1,2,h56) = A56(2,1,1,2,h56) + s3456**(-1)*s456**(-1) * (
     &    za(p2,p4)*za(p2,p8)*za(p4,p6)*zb(p1,p3)*zb(p1,p7)*zb(p4,p5)*
     &    iza(p2,p7)*izb(p1,p8)*s56**(-1)*s78**(-1) + za(p2,p6)*za(p2,
     &    p8)*za(p4,p6)*zb(p1,p3)*zb(p1,p7)*iza(p2,p7)*iza(p5,p6)*izb(
     &    p1,p8)*s78**(-1) )
      A56(2,1,2,2,h56)= + s3456**(-1)*s178**(-1)*s356**(-1) * ( za(p1,
     &    p8)*za(p3,p6)*zb(p1,p7)*zb(p2,p3)*zb(p3,p5)*izb(p1,p8)*zab2(
     &    p4,p1,p8,p7)*s56**(-1)*s78**(-1) - za(p1,p8)*zb(p1,p7)*zb(p2,
     &    p5)*zb(p3,p5)*izb(p1,p8)*izb(p5,p6)*zab2(p4,p1,p8,p7)*
     &    s78**(-1) )
      A56(2,1,2,2,h56) = A56(2,1,2,2,h56) + s3456**(-1)*s178**(-1)*
     & s456**(-1) * (  - za(p1,p8)*za(p4,p6)*zb(p1,p7)*zb(p2,p3)*zb(p4,
     &    p5)*izb(p1,p8)*zab2(p4,p1,p8,p7)*s56**(-1)*s78**(-1) - za(p1,
     &    p8)*za(p4,p6)*zb(p1,p7)*zb(p2,p3)*iza(p5,p6)*izb(p1,p8)*zab2(
     &    p6,p1,p8,p7)*s78**(-1) )
      A56(2,1,2,2,h56) = A56(2,1,2,2,h56) + s3456**(-1)*s278**(-1)*
     & s356**(-1) * ( za(p1,p4)*za(p2,p8)*za(p3,p6)*zb(p2,p7)*zb(p3,p5)
     &    *iza(p2,p7)*zab2(p8,p2,p7,p3)*s56**(-1)*s78**(-1) - za(p1,p4)
     &    *za(p2,p8)*zb(p2,p7)*zb(p3,p5)*iza(p2,p7)*izb(p5,p6)*zab2(p8,
     &    p2,p7,p5)*s78**(-1) )
      A56(2,1,2,2,h56) = A56(2,1,2,2,h56) + s3456**(-1)*s278**(-1)*
     & s456**(-1) * (  - za(p1,p4)*za(p2,p8)*za(p4,p6)*zb(p2,p7)*zb(p4,
     &    p5)*iza(p2,p7)*zab2(p8,p2,p7,p3)*s56**(-1)*s78**(-1) - za(p1,
     &    p6)*za(p2,p8)*za(p4,p6)*zb(p2,p7)*iza(p2,p7)*iza(p5,p6)*zab2(
     &    p8,p2,p7,p3)*s78**(-1) )
      A56(2,1,2,2,h56) = A56(2,1,2,2,h56) + s3456**(-1)*s356**(-1) * (
     &     - za(p3,p6)*zb(p3,p5)*iza(p2,p7)*izb(p1,p8)*zab2(p4,p1,p8,p7
     &    )*zab2(p8,p2,p7,p3)*s56**(-1)*s78**(-1) + zb(p3,p5)*iza(p2,p7
     &    )*izb(p1,p8)*izb(p5,p6)*zab2(p4,p1,p8,p7)*zab2(p8,p2,p7,p5)*
     &    s78**(-1) )
      A56(2,1,2,2,h56) = A56(2,1,2,2,h56) + s3456**(-1)*s456**(-1) * (
     &    za(p4,p6)*zb(p4,p5)*iza(p2,p7)*izb(p1,p8)*zab2(p4,p1,p8,p7)*
     &    zab2(p8,p2,p7,p3)*s56**(-1)*s78**(-1) + za(p4,p6)*iza(p2,p7)*
     &    iza(p5,p6)*izb(p1,p8)*zab2(p6,p1,p8,p7)*zab2(p8,p2,p7,p3)*
     &    s78**(-1) )

      A56(2,2,1,1,h56)= + s3456**(-1)*s178**(-1)*s356**(-1) * (  - za(
     &    p2,p3)*za(p3,p6)*zb(p1,p7)*zb(p3,p5)*iza(p2,p8)*iza(p7,p8)*
     &    zab3(p2,p1,p7,p8,p4)*s56**(-1) - za(p2,p3)*za(p3,p6)*zb(p1,p8
     &    )*zb(p3,p5)*iza(p2,p7)*iza(p7,p8)*zab3(p2,p1,p7,p8,p4)*
     &    s56**(-1) + za(p2,p3)*za(p3,p6)*zb(p3,p5)*iza(p1,p8)*iza(p2,
     &    p7)*iza(p2,p8)*zab2(p2,p1,p8,p7)*zab3(p2,p1,p7,p8,p4)*
     &    s56**(-1) - za(p2,p6)*za(p3,p6)*zb(p1,p7)*iza(p2,p8)*iza(p5,
     &    p6)*iza(p7,p8)*zab3(p2,p1,p7,p8,p4) - za(p2,p6)*za(p3,p6)*zb(
     &    p1,p8)*iza(p2,p7)*iza(p5,p6)*iza(p7,p8)*zab3(p2,p1,p7,p8,p4)
     &     + za(p2,p6)*za(p3,p6)*iza(p1,p8)*iza(p2,p7)*iza(p2,p8)*iza(
     &    p5,p6)*zab2(p2,p1,p8,p7)*zab3(p2,p1,p7,p8,p4) )
      A56(2,2,1,1,h56) = A56(2,2,1,1,h56) + s3456**(-1)*s178**(-1)*
     & s456**(-1) * ( za(p2,p3)*za(p4,p6)*zb(p1,p7)*zb(p4,p5)*iza(p2,p8
     &    )*iza(p7,p8)*zab3(p2,p1,p7,p8,p4)*s56**(-1) + za(p2,p3)*za(p4
     &    ,p6)*zb(p1,p8)*zb(p4,p5)*iza(p2,p7)*iza(p7,p8)*zab3(p2,p1,p7,
     &    p8,p4)*s56**(-1) - za(p2,p3)*za(p4,p6)*zb(p4,p5)*iza(p1,p8)*
     &    iza(p2,p7)*iza(p2,p8)*zab2(p2,p1,p8,p7)*zab3(p2,p1,p7,p8,p4)*
     &    s56**(-1) - za(p2,p3)*zb(p1,p7)*zb(p4,p5)*iza(p2,p8)*iza(p7,
     &    p8)*izb(p5,p6)*zab3(p2,p1,p7,p8,p5) - za(p2,p3)*zb(p1,p8)*zb(
     &    p4,p5)*iza(p2,p7)*iza(p7,p8)*izb(p5,p6)*zab3(p2,p1,p7,p8,p5)
     &     + za(p2,p3)*zb(p4,p5)*iza(p1,p8)*iza(p2,p7)*iza(p2,p8)*izb(
     &    p5,p6)*zab2(p2,p1,p8,p7)*zab3(p2,p1,p7,p8,p5) )
      A56(2,2,2,1,h56)= + s3456**(-1)*s278**(-1)*s356**(-1) * ( za(p1,
     &    p3)*za(p3,p6)*zb(p2,p7)*zb(p3,p5)*iza(p1,p8)*iza(p7,p8)*zab3(
     &    p1,p2,p7,p8,p4)*s56**(-1) + za(p1,p3)*za(p3,p6)*zb(p2,p8)*zb(
     &    p3,p5)*iza(p1,p7)*iza(p7,p8)*zab3(p1,p2,p7,p8,p4)*s56**(-1)
     &     + za(p1,p3)*za(p3,p6)*zb(p3,p5)*iza(p1,p7)*iza(p1,p8)*iza(p2
     &    ,p7)*zab2(p1,p2,p7,p8)*zab3(p1,p2,p7,p8,p4)*s56**(-1) + za(p1
     &    ,p6)*za(p3,p6)*zb(p2,p7)*iza(p1,p8)*iza(p5,p6)*iza(p7,p8)*
     &    zab3(p1,p2,p7,p8,p4) + za(p1,p6)*za(p3,p6)*zb(p2,p8)*iza(p1,
     &    p7)*iza(p5,p6)*iza(p7,p8)*zab3(p1,p2,p7,p8,p4) + za(p1,p6)*
     &    za(p3,p6)*iza(p1,p7)*iza(p1,p8)*iza(p2,p7)*iza(p5,p6)*zab2(p1
     &    ,p2,p7,p8)*zab3(p1,p2,p7,p8,p4) )
      A56(2,2,2,1,h56) = A56(2,2,2,1,h56) + s3456**(-1)*s278**(-1)*
     & s456**(-1) * (  - za(p1,p3)*za(p4,p6)*zb(p2,p7)*zb(p4,p5)*iza(p1
     &    ,p8)*iza(p7,p8)*zab3(p1,p2,p7,p8,p4)*s56**(-1) - za(p1,p3)*
     &    za(p4,p6)*zb(p2,p8)*zb(p4,p5)*iza(p1,p7)*iza(p7,p8)*zab3(p1,
     &    p2,p7,p8,p4)*s56**(-1) - za(p1,p3)*za(p4,p6)*zb(p4,p5)*iza(p1
     &    ,p7)*iza(p1,p8)*iza(p2,p7)*zab2(p1,p2,p7,p8)*zab3(p1,p2,p7,p8
     &    ,p4)*s56**(-1) + za(p1,p3)*zb(p2,p7)*zb(p4,p5)*iza(p1,p8)*
     &    iza(p7,p8)*izb(p5,p6)*zab3(p1,p2,p7,p8,p5) + za(p1,p3)*zb(p2,
     &    p8)*zb(p4,p5)*iza(p1,p7)*iza(p7,p8)*izb(p5,p6)*zab3(p1,p2,p7,
     &    p8,p5) + za(p1,p3)*zb(p4,p5)*iza(p1,p7)*iza(p1,p8)*iza(p2,p7)
     &    *izb(p5,p6)*zab2(p1,p2,p7,p8)*zab3(p1,p2,p7,p8,p5) )
      A56(2,2,1,2,h56)= + s3456**(-1)*s178**(-1)*s356**(-1) * (  - za(
     &    p2,p4)*za(p3,p6)*zb(p1,p7)*zb(p3,p5)*iza(p2,p8)*iza(p7,p8)*
     &    zab3(p2,p1,p7,p8,p3)*s56**(-1) - za(p2,p4)*za(p3,p6)*zb(p1,p8
     &    )*zb(p3,p5)*iza(p2,p7)*iza(p7,p8)*zab3(p2,p1,p7,p8,p3)*
     &    s56**(-1) + za(p2,p4)*za(p3,p6)*zb(p3,p5)*iza(p1,p8)*iza(p2,
     &    p7)*iza(p2,p8)*zab2(p2,p1,p8,p7)*zab3(p2,p1,p7,p8,p3)*
     &    s56**(-1) + za(p2,p4)*zb(p1,p7)*zb(p3,p5)*iza(p2,p8)*iza(p7,
     &    p8)*izb(p5,p6)*zab3(p2,p1,p7,p8,p5) + za(p2,p4)*zb(p1,p8)*zb(
     &    p3,p5)*iza(p2,p7)*iza(p7,p8)*izb(p5,p6)*zab3(p2,p1,p7,p8,p5)
     &     - za(p2,p4)*zb(p3,p5)*iza(p1,p8)*iza(p2,p7)*iza(p2,p8)*izb(
     &    p5,p6)*zab2(p2,p1,p8,p7)*zab3(p2,p1,p7,p8,p5) )
      A56(2,2,1,2,h56) = A56(2,2,1,2,h56) + s3456**(-1)*s178**(-1)*
     & s456**(-1) * ( za(p2,p4)*za(p4,p6)*zb(p1,p7)*zb(p4,p5)*iza(p2,p8
     &    )*iza(p7,p8)*zab3(p2,p1,p7,p8,p3)*s56**(-1) + za(p2,p4)*za(p4
     &    ,p6)*zb(p1,p8)*zb(p4,p5)*iza(p2,p7)*iza(p7,p8)*zab3(p2,p1,p7,
     &    p8,p3)*s56**(-1) - za(p2,p4)*za(p4,p6)*zb(p4,p5)*iza(p1,p8)*
     &    iza(p2,p7)*iza(p2,p8)*zab2(p2,p1,p8,p7)*zab3(p2,p1,p7,p8,p3)*
     &    s56**(-1) + za(p2,p6)*za(p4,p6)*zb(p1,p7)*iza(p2,p8)*iza(p5,
     &    p6)*iza(p7,p8)*zab3(p2,p1,p7,p8,p3) + za(p2,p6)*za(p4,p6)*zb(
     &    p1,p8)*iza(p2,p7)*iza(p5,p6)*iza(p7,p8)*zab3(p2,p1,p7,p8,p3)
     &     - za(p2,p6)*za(p4,p6)*iza(p1,p8)*iza(p2,p7)*iza(p2,p8)*iza(
     &    p5,p6)*zab2(p2,p1,p8,p7)*zab3(p2,p1,p7,p8,p3) )
      A56(2,2,2,2,h56)= + s3456**(-1)*s278**(-1)*s356**(-1) * ( za(p1,
     &    p4)*za(p3,p6)*zb(p2,p7)*zb(p3,p5)*iza(p1,p8)*iza(p7,p8)*zab3(
     &    p1,p2,p7,p8,p3)*s56**(-1) + za(p1,p4)*za(p3,p6)*zb(p2,p8)*zb(
     &    p3,p5)*iza(p1,p7)*iza(p7,p8)*zab3(p1,p2,p7,p8,p3)*s56**(-1)
     &     + za(p1,p4)*za(p3,p6)*zb(p3,p5)*iza(p1,p7)*iza(p1,p8)*iza(p2
     &    ,p7)*zab2(p1,p2,p7,p8)*zab3(p1,p2,p7,p8,p3)*s56**(-1) - za(p1
     &    ,p4)*zb(p2,p7)*zb(p3,p5)*iza(p1,p8)*iza(p7,p8)*izb(p5,p6)*
     &    zab3(p1,p2,p7,p8,p5) - za(p1,p4)*zb(p2,p8)*zb(p3,p5)*iza(p1,
     &    p7)*iza(p7,p8)*izb(p5,p6)*zab3(p1,p2,p7,p8,p5) - za(p1,p4)*
     &    zb(p3,p5)*iza(p1,p7)*iza(p1,p8)*iza(p2,p7)*izb(p5,p6)*zab2(p1
     &    ,p2,p7,p8)*zab3(p1,p2,p7,p8,p5) )
      A56(2,2,2,2,h56) = A56(2,2,2,2,h56) + s3456**(-1)*s278**(-1)*
     & s456**(-1) * (  - za(p1,p4)*za(p4,p6)*zb(p2,p7)*zb(p4,p5)*iza(p1
     &    ,p8)*iza(p7,p8)*zab3(p1,p2,p7,p8,p3)*s56**(-1) - za(p1,p4)*
     &    za(p4,p6)*zb(p2,p8)*zb(p4,p5)*iza(p1,p7)*iza(p7,p8)*zab3(p1,
     &    p2,p7,p8,p3)*s56**(-1) - za(p1,p4)*za(p4,p6)*zb(p4,p5)*iza(p1
     &    ,p7)*iza(p1,p8)*iza(p2,p7)*zab2(p1,p2,p7,p8)*zab3(p1,p2,p7,p8
     &    ,p3)*s56**(-1) - za(p1,p6)*za(p4,p6)*zb(p2,p7)*iza(p1,p8)*
     &    iza(p5,p6)*iza(p7,p8)*zab3(p1,p2,p7,p8,p3) - za(p1,p6)*za(p4,
     &    p6)*zb(p2,p8)*iza(p1,p7)*iza(p5,p6)*iza(p7,p8)*zab3(p1,p2,p7,
     &    p8,p3) - za(p1,p6)*za(p4,p6)*iza(p1,p7)*iza(p1,p8)*iza(p2,p7)
     &    *iza(p5,p6)*zab2(p1,p2,p7,p8)*zab3(p1,p2,p7,p8,p3) )

c     enddo

      return
      end
