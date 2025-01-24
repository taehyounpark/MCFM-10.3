!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a7ZZsrn(p1,p2,p3,p4,p5,p6,p7,
     & p,n,za,zb,A34,A56)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer p1,p2,p3,p4,p5,p6,p7,h34,h56
      complex(dp)::A34(2,2,2),A56(2,2,2),temp(2,2,2)
      real(dp):: n(4),p(mxpart,4)
c     result for the ZZ process supplementary diagrams
c     with a factor of 2*rt2*e^4*gs removed
c     setting up color orders andcontribution 34<-->56
c     Amplitudes for the process
c     0 --> qbar(j1)+q(j2)+e^-(j3)+e^+(p4)+mu^-(p5)+mu^+(p6)+g(p7)
c     only s34,s56,propagators included.
c     pol=1=LH, pol=2=RH

c     part with propagator s34
      call a7ZZsrnpart(p1,p2,p5,p6,p3,p4,p7,p,n,za,zb,temp)
      do h34=1,2
      do h56=1,2
      A34(:,h34,h56)= temp(:,h56,h34)
      enddo
      enddo

c     part with propagator s56 (i.e. unswapped)
      call a7ZZsrnpart(p1,p2,p3,p4,p5,p6,p7,p,n,za,zb,A56)

      return
      end


      subroutine a7ZZsrnpart(p1,p2,p3,p4,p5,p6,p7,
     & p,n,za,zb,A)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer::p1,p2,p3,p4,p5,p6,p7,i,j
      real(dp):: s356,s456,s3456,s17,s27,s56,s3,n(4),p(mxpart,4)
      complex(dp):: A(2,2,2),zab2,vecm(mxpart,mxpart)
c     result for the ZZ process supplementary diagrams
c     with a factor of 2*rt2*e^4*gs removed
c     contribution 34<-->56 needs to be calculated
c     Amplitudes for the process
c     0 --> qbar(j1)+q(j2)+e^-(j3)+e^+(p4)+mu^-(p5)+mu^+(p6)+g(p7)
c     only s34,s56,propagators included.
c     pol=1=LH, pol=2=RH

      s3(p1,p2,p3)=s(p1,p2)+s(p1,p3)+s(p2,p3)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      s17=s(p1,p7)
      s27=s(p2,p7)
      s56=s(p5,p6)
      s356=s3(p3,p5,p6)
      s456=s3(p4,p5,p6)
      s3456=s(p3,p4)+s(p3,p5)+s(p3,p6)+s(p4,p5)+s(p4,p6)+s(p5,p6)
      do i=1,7
      do j=i,7
      call ndveccur(i,j,n,p,vecm)
      enddo
      enddo

      A(1,1,1)= + s356**(-1)*s17**(-1) * ( za(p3,p5)*zb(p1,p4)*zab2(p2,
     &    p3,p5,p6)*vecm(p1,p1) + za(p3,p5)*zb(p7,p4)*zab2(p2,p3,p5,p6)
     &    *vecm(p7,p1) )
      A(1,1,1) = A(1,1,1) + s356**(-1)*s27**(-1) * (  - za(p3,p5)*zb(p1
     &    ,p4)*zab2(p2,p3,p5,p6)*vecm(p2,p2) - za(p3,p5)*zb(p1,p4)*
     &    zab2(p7,p3,p5,p6)*vecm(p2,p7) )
      A(1,1,1) = A(1,1,1) + s456**(-1)*s17**(-1) * ( za(p3,p2)*zb(p4,p6
     &    )*zab2(p5,p4,p6,p1)*vecm(p1,p1) + za(p3,p2)*zb(p4,p6)*zab2(p5
     &    ,p4,p6,p7)*vecm(p7,p1) )
      A(1,1,1) = A(1,1,1) + s456**(-1)*s27**(-1) * (  - za(p3,p2)*zb(p4
     &    ,p6)*zab2(p5,p4,p6,p1)*vecm(p2,p2) - za(p3,p7)*zb(p4,p6)*
     &    zab2(p5,p4,p6,p1)*vecm(p2,p7) )
      A(2,1,1)= + s356**(-1)*s17**(-1) * ( za(p3,p5)*zb(p2,p4)*zab2(p1,
     &    p3,p5,p6)*vecm(p1,p1) + za(p3,p5)*zb(p2,p4)*zab2(p7,p3,p5,p6)
     &    *vecm(p1,p7) )
      A(2,1,1) = A(2,1,1) + s356**(-1)*s27**(-1) * (  - za(p3,p5)*zb(p2
     &    ,p4)*zab2(p1,p3,p5,p6)*vecm(p2,p2) - za(p3,p5)*zb(p7,p4)*
     &    zab2(p1,p3,p5,p6)*vecm(p7,p2) )
      A(2,1,1) = A(2,1,1) + s456**(-1)*s17**(-1) * ( za(p3,p1)*zb(p4,p6
     &    )*zab2(p5,p4,p6,p2)*vecm(p1,p1) + za(p3,p7)*zb(p4,p6)*zab2(p5
     &    ,p4,p6,p2)*vecm(p1,p7) )
      A(2,1,1) = A(2,1,1) + s456**(-1)*s27**(-1) * (  - za(p3,p1)*zb(p4
     &    ,p6)*zab2(p5,p4,p6,p2)*vecm(p2,p2) - za(p3,p1)*zb(p4,p6)*
     &    zab2(p5,p4,p6,p7)*vecm(p7,p2) )
      A(1,2,1)= + s356**(-1)*s17**(-1) * ( za(p2,p4)*zb(p3,p6)*zab2(p5,
     &    p3,p6,p1)*vecm(p1,p1) + za(p2,p4)*zb(p3,p6)*zab2(p5,p3,p6,p7)
     &    *vecm(p7,p1) )
      A(1,2,1) = A(1,2,1) + s356**(-1)*s27**(-1) * (  - za(p2,p4)*zb(p3
     &    ,p6)*zab2(p5,p3,p6,p1)*vecm(p2,p2) - za(p7,p4)*zb(p3,p6)*
     &    zab2(p5,p3,p6,p1)*vecm(p2,p7) )
      A(1,2,1) = A(1,2,1) + s456**(-1)*s17**(-1) * ( za(p4,p5)*zb(p3,p1
     &    )*zab2(p2,p4,p5,p6)*vecm(p1,p1) + za(p4,p5)*zb(p3,p7)*zab2(p2
     &    ,p4,p5,p6)*vecm(p7,p1) )
      A(1,2,1) = A(1,2,1) + s456**(-1)*s27**(-1) * (  - za(p4,p5)*zb(p3
     &    ,p1)*zab2(p2,p4,p5,p6)*vecm(p2,p2) - za(p4,p5)*zb(p3,p1)*
     &    zab2(p7,p4,p5,p6)*vecm(p2,p7) )
      A(2,2,1)= + s356**(-1)*s17**(-1) * ( za(p1,p4)*zb(p3,p6)*zab2(p5,
     &    p3,p6,p2)*vecm(p1,p1) + za(p7,p4)*zb(p3,p6)*zab2(p5,p3,p6,p2)
     &    *vecm(p1,p7) )
      A(2,2,1) = A(2,2,1) + s356**(-1)*s27**(-1) * (  - za(p1,p4)*zb(p3
     &    ,p6)*zab2(p5,p3,p6,p2)*vecm(p2,p2) - za(p1,p4)*zb(p3,p6)*
     &    zab2(p5,p3,p6,p7)*vecm(p7,p2) )
      A(2,2,1) = A(2,2,1) + s456**(-1)*s17**(-1) * ( za(p4,p5)*zb(p3,p2
     &    )*zab2(p1,p4,p5,p6)*vecm(p1,p1) + za(p4,p5)*zb(p3,p2)*zab2(p7
     &    ,p4,p5,p6)*vecm(p1,p7) )
      A(2,2,1) = A(2,2,1) + s456**(-1)*s27**(-1) * (  - za(p4,p5)*zb(p3
     &    ,p2)*zab2(p1,p4,p5,p6)*vecm(p2,p2) - za(p4,p5)*zb(p3,p7)*
     &    zab2(p1,p4,p5,p6)*vecm(p7,p2) )

      A(1,1,2)= + s356**(-1)*s17**(-1) * ( za(p3,p6)*zb(p1,p4)*zab2(p2,
     &    p3,p6,p5)*vecm(p1,p1) + za(p3,p6)*zb(p7,p4)*zab2(p2,p3,p6,p5)
     &    *vecm(p7,p1) )
      A(1,1,2) = A(1,1,2) + s356**(-1)*s27**(-1) * (  - za(p3,p6)*zb(p1
     &    ,p4)*zab2(p2,p3,p6,p5)*vecm(p2,p2) - za(p3,p6)*zb(p1,p4)*
     &    zab2(p7,p3,p6,p5)*vecm(p2,p7) )
      A(1,1,2) = A(1,1,2) + s456**(-1)*s17**(-1) * ( za(p3,p2)*zb(p4,p5
     &    )*zab2(p6,p4,p5,p1)*vecm(p1,p1) + za(p3,p2)*zb(p4,p5)*zab2(p6
     &    ,p4,p5,p7)*vecm(p7,p1) )
      A(1,1,2) = A(1,1,2) + s456**(-1)*s27**(-1) * (  - za(p3,p2)*zb(p4
     &    ,p5)*zab2(p6,p4,p5,p1)*vecm(p2,p2) - za(p3,p7)*zb(p4,p5)*
     &    zab2(p6,p4,p5,p1)*vecm(p2,p7) )
      A(2,1,2)= + s356**(-1)*s17**(-1) * ( za(p3,p6)*zb(p2,p4)*zab2(p1,
     &    p3,p6,p5)*vecm(p1,p1) + za(p3,p6)*zb(p2,p4)*zab2(p7,p3,p6,p5)
     &    *vecm(p1,p7) )
      A(2,1,2) = A(2,1,2) + s356**(-1)*s27**(-1) * (  - za(p3,p6)*zb(p2
     &    ,p4)*zab2(p1,p3,p6,p5)*vecm(p2,p2) - za(p3,p6)*zb(p7,p4)*
     &    zab2(p1,p3,p6,p5)*vecm(p7,p2) )
      A(2,1,2) = A(2,1,2) + s456**(-1)*s17**(-1) * ( za(p3,p1)*zb(p4,p5
     &    )*zab2(p6,p4,p5,p2)*vecm(p1,p1) + za(p3,p7)*zb(p4,p5)*zab2(p6
     &    ,p4,p5,p2)*vecm(p1,p7) )
      A(2,1,2) = A(2,1,2) + s456**(-1)*s27**(-1) * (  - za(p3,p1)*zb(p4
     &    ,p5)*zab2(p6,p4,p5,p2)*vecm(p2,p2) - za(p3,p1)*zb(p4,p5)*
     &    zab2(p6,p4,p5,p7)*vecm(p7,p2) )
      A(1,2,2)= + s356**(-1)*s17**(-1) * ( za(p2,p4)*zb(p3,p5)*zab2(p6,
     &    p3,p5,p1)*vecm(p1,p1) + za(p2,p4)*zb(p3,p5)*zab2(p6,p3,p5,p7)
     &    *vecm(p7,p1) )
      A(1,2,2) = A(1,2,2) + s356**(-1)*s27**(-1) * (  - za(p2,p4)*zb(p3
     &    ,p5)*zab2(p6,p3,p5,p1)*vecm(p2,p2) - za(p7,p4)*zb(p3,p5)*
     &    zab2(p6,p3,p5,p1)*vecm(p2,p7) )
      A(1,2,2) = A(1,2,2) + s456**(-1)*s17**(-1) * ( za(p4,p6)*zb(p3,p1
     &    )*zab2(p2,p4,p6,p5)*vecm(p1,p1) + za(p4,p6)*zb(p3,p7)*zab2(p2
     &    ,p4,p6,p5)*vecm(p7,p1) )
      A(1,2,2) = A(1,2,2) + s456**(-1)*s27**(-1) * (  - za(p4,p6)*zb(p3
     &    ,p1)*zab2(p2,p4,p6,p5)*vecm(p2,p2) - za(p4,p6)*zb(p3,p1)*
     &    zab2(p7,p4,p6,p5)*vecm(p2,p7) )
      A(2,2,2)= + s356**(-1)*s17**(-1) * ( za(p1,p4)*zb(p3,p5)*zab2(p6,
     &    p3,p5,p2)*vecm(p1,p1) + za(p7,p4)*zb(p3,p5)*zab2(p6,p3,p5,p2)
     &    *vecm(p1,p7) )
      A(2,2,2) = A(2,2,2) + s356**(-1)*s27**(-1) * (  - za(p1,p4)*zb(p3
     &    ,p5)*zab2(p6,p3,p5,p2)*vecm(p2,p2) - za(p1,p4)*zb(p3,p5)*
     &    zab2(p6,p3,p5,p7)*vecm(p7,p2) )
      A(2,2,2) = A(2,2,2) + s456**(-1)*s17**(-1) * ( za(p4,p6)*zb(p3,p2
     &    )*zab2(p1,p4,p6,p5)*vecm(p1,p1) + za(p4,p6)*zb(p3,p2)*zab2(p7
     &    ,p4,p6,p5)*vecm(p1,p7) )
      A(2,2,2) = A(2,2,2) + s456**(-1)*s27**(-1) * (  - za(p4,p6)*zb(p3
     &    ,p2)*zab2(p1,p4,p6,p5)*vecm(p2,p2) - za(p4,p6)*zb(p3,p7)*
     &    zab2(p1,p4,p6,p5)*vecm(p7,p2) )

      A(:,:,:)=A(:,:,:)/(s3456*s56)
      return
      end
