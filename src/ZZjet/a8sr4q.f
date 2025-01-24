!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a8sr4q(p1,p2,p3,p4,p5,p6,p7,p8,
     & u3456,l3456,u5634,l5634)
      implicit none
      include 'types.f'
      integer p1,p2,p3,p4,p5,p6,p7,p8,h1,h2
      complex(dp):: u3456(2,2,2,2),l3456(2,2,2,2),
     & u5634(2,2,2,2),l5634(2,2,2,2),temp1(2,2,2,2),
     & temp2(2,2,2,2)
c     u3456(2,2,2,2) equiv u3456(polg17,pol28,pol34,pol56)
c     result for the ZZ process four quark supplementary diagrams
c     with a factor of i*4*e^4*gs^2 removed
c     setting up exchange contribution 34<-->56
c     Amplitudes for the process
c 0-->qbar(j1)+q(j2)+e^-(j3)+e^+(p4)+mu^-(p5)+mu^+(p6)+g(p7)+g(p8)
c     only s3456,s56,propagators included.
c     pol=1=LH, pol=2=RH
      call a8sr4qpart(p1,p2,p3,p4,p5,p6,p7,p8,u3456,l3456)

      call a8sr4qpart(p1,p2,p5,p6,p3,p4,p7,p8,temp1,temp2)
      do h1=1,2
      do h2=1,2
      u5634(:,:,h1,h2)= temp1(:,:,h2,h1)
      l5634(:,:,h1,h2)= temp2(:,:,h2,h1)
      enddo
      enddo

      return
      end


      subroutine a8sr4qpart(p1,p2,p3,p4,j5,j6,p7,p8,upper,lower)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6,p7,p8,j5,j6,h56
      real(dp)::s17,s28,s56,s3,s356,s456,s178,s278,s128,s127,s3456
      complex(dp)::upper(2,2,2,2),lower(2,2,2,2),zab2
c     result for the ZZ process supplementary diagrams
c     with a factor of i*4*e^4*gs^2 removed
c     sequential contribution 34-->56 only
c     Amplitudes for the process
c 0-->qbar(j1)+q(j7)+e^-(j3)+e^+(p4)+mu^-(p5)+mu^+(p6)+qbar(p2)+q(p8)
c     only s3456,s56,propagators included.
c     pol=1=LH, pol=2=RH
c     A(h17,h28,h34,h56)
c Statement functions
      s3(p1,p2,p3)=s(p1,p2)+s(p1,p3)+s(p2,p3)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
c End statement functions
      s17=s(p1,p7)
      s28=s(p2,p8)
      s127=s3(p1,p2,p7)
      s128=s3(p1,p2,p8)
      s178=s3(p1,p7,p8)
      s278=s3(p2,p7,p8)
      do h56=1,2
      if (h56 == 1) then
      p5=j5
      p6=j6
      elseif (h56 == 2) then
      p5=j6
      p6=j5
      endif
      s56=s(p5,p6)
      s356=s3(p3,p5,p6)
      s456=s3(p4,p5,p6)
      s3456=s(p3,p4)+s(p3,p5)+s(p3,p6)+s(p4,p5)+s(p4,p6)+s(p5,p6)

      upper(1,1,1,h56)=
     & -za(p7,p8)*za(p3,p5)*zb(p1,p4)/(s278*s356)
     & *(zb(p6,p3)*zab2(p3,p7,p8,p2)+zb(p6,p5)*zab2(p5,p7,p8,p2))
     & -za(p7,p8)*zb(p4,p6)*zab2(p3,p7,p8,p2)/(s278*s456)
     & *(za(p5,p4)*zb(p4,p1)+za(p5,p6)*zb(p6,p1))
     & -za(p3,p5)*zb(p1,p2)*zab2(p8,p1,p2,p4)/(s128*s356)
     & *(za(p7,p3)*zb(p3,p6)+za(p7,p5)*zb(p5,p6))
     & -za(p3,p7)*zb(p1,p2)*zb(p4,p6)/(s128*s456)
     & *(zab2(p8,p1,p2,p4)*za(p4,p5)+zab2(p8,p1,p2,p6)*za(p6,p5))

      upper(1,1,2,h56)=
     & +za(p7,p8)*zb(p3,p6)*zab2(p4,p7,p8,p2)/(s278*s356)
     & *zab2(p5,p3,p6,p1)
     & +za(p4,p5)*za(p7,p8)*zb(p1,p3)/(s278*s456)
     & *(zb(p6,p4)*zab2(p4,p7,p8,p2)+zb(p6,p5)*zab2(p5,p7,p8,p2))
     & +za(p4,p7)*zb(p1,p2)*zb(p3,p6)/(s128*s356)
     & *(zab2(p8,p1,p2,p3)*za(p3,p5)+zab2(p8,p1,p2,p6)*za(p6,p5))
     & +za(p4,p5)*zb(p1,p2)*zab2(p8,p1,p2,p3)/(s128*s456)
     & *zab2(p7,p4,p5,p6)

      upper(1,2,1,h56)=
     & +za(p2,p7)*za(p3,p5)*zb(p1,p4)/(s278*s356)
     & *(zb(p6,p3)*zab2(p3,p2,p7,p8)+zb(p6,p5)*zab2(p5,p2,p7,p8))
     & +zab2(p3,p2,p7,p8)*za(p2,p7)*zb(p4,p6)/(s278*s456)
     & *zab2(p5,p4,p6,p1)
     & -zab2(p2,p1,p8,p4)*za(p3,p5)*zb(p1,p8)/(s128*s356)
     & *zab2(p7,p3,p5,p6)
     & -za(p3,p7)*zb(p1,p8)*zb(p4,p6)/(s128*s456)
     & *(zab2(p2,p1,p8,p4)*za(p4,p5)+zab2(p2,p1,p8,p6)*za(p6,p5))

      upper(1,2,2,h56)=
     & -zab2(p4,p2,p7,p8)*za(p2,p7)*zb(p3,p6)/(s278*s356)
     & *zab2(p5,p3,p6,p1)
     & -za(p2,p7)*za(p4,p5)*zb(p1,p3)/(s278*s456)
     & *(zb(p6,p4)*zab2(p4,p2,p7,p8)+zb(p6,p5)*zab2(p5,p2,p7,p8))
     & +za(p4,p7)*zb(p1,p8)*zb(p3,p6)/(s128*s356)
     & *(zab2(p2,p1,p8,p3)*za(p3,p5)+zab2(p2,p1,p8,p6)*za(p6,p5))
     & +zab2(p2,p1,p8,p3)*za(p4,p5)*zb(p1,p8)/(s128*s456)
     & *zab2(p7,p4,p5,p6)

      upper(2,1,1,h56)=
     & -zab2(p8,p2,p7,p4)*za(p3,p5)*zb(p2,p7)/(s278*s356)
     & *zab2(p1,p3,p5,p6)
     & +za(p1,p3)*zb(p2,p7)*zb(p4,p6)/(s278*s456)
     & *(zab2(p8,p2,p7,p4)*za(p4,p5)+zab2(p8,p2,p7,p6)*za(p6,p5))
     & -za(p1,p8)*za(p3,p5)*zb(p4,p7)/(s128*s356)
     & *(zb(p6,p3)*zab2(p3,p1,p8,p2)+zb(p6,p5)*zab2(p5,p1,p8,p2))
     & +zab2(p3,p1,p8,p2)*za(p1,p8)*zb(p4,p6)/(s128*s456)
     & *zab2(p5,p4,p6,p7)

      upper(2,1,2,h56)=
     & -za(p1,p4)*zb(p2,p7)*zb(p3,p6)/(s278*s356)
     & *(zab2(p8,p2,p7,p3)*za(p3,p5)+zab2(p8,p2,p7,p6)*za(p6,p5))
     & +zab2(p8,p2,p7,p3)*za(p4,p5)*zb(p2,p7)/(s278*s456)
     & *zab2(p1,p4,p5,p6)
     & -zab2(p4,p1,p8,p2)*za(p1,p8)*zb(p3,p6)/(s128*s356)
     & *zab2(p5,p3,p6,p7)
     & +za(p1,p8)*za(p4,p5)*zb(p3,p7)/(s128*s456)
     & *(zb(p6,p4)*zab2(p4,p1,p8,p2)+zb(p6,p5)*zab2(p5,p1,p8,p2))

      upper(2,2,1,h56)=
     & +za(p3,p5)*zb(p7,p8)*zab2(p2,p7,p8,p4)/(s278*s356)
     & *zab2(p1,p3,p5,p6)
     & -za(p1,p3)*zb(p4,p6)*zb(p7,p8)/(s278*s456)
     & *(zab2(p2,p7,p8,p4)*za(p4,p5)+zab2(p2,p7,p8,p6)*za(p6,p5))
     & -za(p1,p2)*za(p3,p5)*zb(p4,p7)/(s128*s356)
     & *(zb(p6,p3)*zab2(p3,p1,p2,p8)+zb(p6,p5)*zab2(p5,p1,p2,p8))
     & +za(p1,p2)*zb(p4,p6)*zab2(p3,p1,p2,p8)/(s128*s456)
     & *zab2(p5,p4,p6,p7)

      upper(2,2,2,h56)=
     & +za(p1,p4)*zb(p3,p6)*zb(p7,p8)/(s278*s356)
     & *(zab2(p2,p7,p8,p3)*za(p3,p5)+zab2(p2,p7,p8,p6)*za(p6,p5))
     & -za(p4,p5)*zb(p7,p8)*zab2(p2,p7,p8,p3)/(s278*s456)
     & *zab2(p1,p4,p5,p6)
     & -za(p1,p2)*zb(p3,p6)*zab2(p4,p1,p2,p8)/(s128*s356)
     & *zab2(p5,p3,p6,p7)
     & +za(p1,p2)*za(p4,p5)*zb(p3,p7)/(s128*s456)
     & *(zb(p6,p4)*zab2(p4,p1,p2,p8)+zb(p6,p5)*zab2(p5,p1,p2,p8))

      lower(1,1,1,h56)=
     & -za(p8,p7)*za(p3,p5)*zb(p2,p4)/(s178*s356)
     & *(zb(p6,p3)*zab2(p3,p8,p7,p1)+zb(p6,p5)*zab2(p5,p8,p7,p1))
     & -za(p8,p7)*zb(p4,p6)*zab2(p3,p8,p7,p1)/(s178*s456)
     & *(za(p5,p4)*zb(p4,p2)+za(p5,p6)*zb(p6,p2))
     & -za(p3,p5)*zb(p2,p1)*zab2(p7,p2,p1,p4)/(s127*s356)
     & *(za(p8,p3)*zb(p3,p6)+za(p8,p5)*zb(p5,p6))
     & -za(p3,p8)*zb(p2,p1)*zb(p4,p6)/(s127*s456)
     & *(zab2(p7,p2,p1,p4)*za(p4,p5)+zab2(p7,p2,p1,p6)*za(p6,p5))

      lower(1,1,2,h56)=
     & +za(p8,p7)*zb(p3,p6)*zab2(p4,p8,p7,p1)/(s178*s356)
     & *zab2(p5,p3,p6,p2)
     & +za(p4,p5)*za(p8,p7)*zb(p2,p3)/(s178*s456)
     & *(zb(p6,p4)*zab2(p4,p8,p7,p1)+zb(p6,p5)*zab2(p5,p8,p7,p1))
     & +za(p4,p8)*zb(p2,p1)*zb(p3,p6)/(s127*s356)
     & *(zab2(p7,p2,p1,p3)*za(p3,p5)+zab2(p7,p2,p1,p6)*za(p6,p5))
     & +za(p4,p5)*zb(p2,p1)*zab2(p7,p2,p1,p3)/(s127*s456)
     & *zab2(p8,p4,p5,p6)

      lower(2,1,1,h56)=
     & +za(p1,p8)*za(p3,p5)*zb(p2,p4)/(s178*s356)
     & *(zb(p6,p3)*zab2(p3,p1,p8,p7)+zb(p6,p5)*zab2(p5,p1,p8,p7))
     & +zab2(p3,p1,p8,p7)*za(p1,p8)*zb(p4,p6)/(s178*s456)
     & *zab2(p5,p4,p6,p2)
     & -zab2(p1,p2,p7,p4)*za(p3,p5)*zb(p2,p7)/(s127*s356)
     & *zab2(p8,p3,p5,p6)
     & -za(p3,p8)*zb(p2,p7)*zb(p4,p6)/(s127*s456)
     & *(zab2(p1,p2,p7,p4)*za(p4,p5)+zab2(p1,p2,p7,p6)*za(p6,p5))

      lower(2,1,2,h56)=
     & -zab2(p4,p1,p8,p7)*za(p1,p8)*zb(p3,p6)/(s178*s356)
     & *zab2(p5,p3,p6,p2)
     & -za(p1,p8)*za(p4,p5)*zb(p2,p3)/(s178*s456)
     & *(zb(p6,p4)*zab2(p4,p1,p8,p7)+zb(p6,p5)*zab2(p5,p1,p8,p7))
     & +za(p4,p8)*zb(p2,p7)*zb(p3,p6)/(s127*s356)
     & *(zab2(p1,p2,p7,p3)*za(p3,p5)+zab2(p1,p2,p7,p6)*za(p6,p5))
     & +zab2(p1,p2,p7,p3)*za(p4,p5)*zb(p2,p7)/(s127*s456)
     & *zab2(p8,p4,p5,p6)

      lower(1,2,1,h56)=
     & -zab2(p7,p1,p8,p4)*za(p3,p5)*zb(p1,p8)/(s178*s356)
     & *zab2(p2,p3,p5,p6)
     & +za(p2,p3)*zb(p1,p8)*zb(p4,p6)/(s178*s456)
     & *(zab2(p7,p1,p8,p4)*za(p4,p5)+zab2(p7,p1,p8,p6)*za(p6,p5))
     & -za(p2,p7)*za(p3,p5)*zb(p4,p8)/(s127*s356)
     & *(zb(p6,p3)*zab2(p3,p2,p7,p1)+zb(p6,p5)*zab2(p5,p2,p7,p1))
     & +zab2(p3,p2,p7,p1)*za(p2,p7)*zb(p4,p6)/(s127*s456)
     & *zab2(p5,p4,p6,p8)

      lower(1,2,2,h56)=
     & -za(p2,p4)*zb(p1,p8)*zb(p3,p6)/(s178*s356)
     & *(zab2(p7,p1,p8,p3)*za(p3,p5)+zab2(p7,p1,p8,p6)*za(p6,p5))
     & +zab2(p7,p1,p8,p3)*za(p4,p5)*zb(p1,p8)/(s178*s456)
     & *zab2(p2,p4,p5,p6)
     & -zab2(p4,p2,p7,p1)*za(p2,p7)*zb(p3,p6)/(s127*s356)
     & *zab2(p5,p3,p6,p8)
     & +za(p2,p7)*za(p4,p5)*zb(p3,p8)/(s127*s456)
     & *(zb(p6,p4)*zab2(p4,p2,p7,p1)+zb(p6,p5)*zab2(p5,p2,p7,p1))

      lower(2,2,1,h56)=
     & +za(p3,p5)*zb(p8,p7)*zab2(p1,p8,p7,p4)/(s178*s356)
     & *zab2(p2,p3,p5,p6)
     & -za(p2,p3)*zb(p4,p6)*zb(p8,p7)/(s178*s456)
     & *(zab2(p1,p8,p7,p4)*za(p4,p5)+zab2(p1,p8,p7,p6)*za(p6,p5))
     & -za(p2,p1)*za(p3,p5)*zb(p4,p8)/(s127*s356)
     & *(zb(p6,p3)*zab2(p3,p2,p1,p7)+zb(p6,p5)*zab2(p5,p2,p1,p7))
     & +za(p2,p1)*zb(p4,p6)*zab2(p3,p2,p1,p7)/(s127*s456)
     & *zab2(p5,p4,p6,p8)

      lower(2,2,2,h56)=
     & +za(p2,p4)*zb(p3,p6)*zb(p8,p7)/(s178*s356)
     & *(zab2(p1,p8,p7,p3)*za(p3,p5)+zab2(p1,p8,p7,p6)*za(p6,p5))
     & -za(p4,p5)*zb(p8,p7)*zab2(p1,p8,p7,p3)/(s178*s456)
     & *zab2(p2,p4,p5,p6)
     & -za(p2,p1)*zb(p3,p6)*zab2(p4,p2,p1,p7)/(s127*s356)
     & *zab2(p5,p3,p6,p8)
     & +za(p2,p1)*za(p4,p5)*zb(p3,p8)/(s127*s456)
     & *(zb(p6,p4)*zab2(p4,p2,p1,p7)+zb(p6,p5)*zab2(p5,p2,p1,p7))
      enddo

      upper(:,:,:,:)=upper(:,:,:,:)/(s56*s28*s3456)
      lower(:,:,:,:)=lower(:,:,:,:)/(s56*s17*s3456)

      return
      end
