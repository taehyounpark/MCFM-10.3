!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqbggAxslCoeffs(p1,p2,p3,p4,p5,p6,za,zb,mtsq,
     &   coeffsl)
      implicit none
      include 'types.f'
c Qed like coefficients
      complex(dp)::coeffsl(15,2,2)
      include 'qqbggintnames.f'


c              6---<----- 5
c                    (
c                    )
c                    (
c                    )
c                   / \
c                  /   \
c                 /     \
c       1|       /       \
c        |      /         \
c   3    |cccccc-----<-----cccccccc 4
c   ccccc|
c       2|


cCCC  QED delta function piece
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6,j,k
      complex(dp)::zab2,iza,izb,ans34(2,2),ans43(2,2)
      real(dp)::s3,s56,s123,s124,mtsq

c--- begin statement function
      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p3,p1)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      izb(p1,p2)=cone/zb(p1,p2)
      iza(p1,p2)=cone/za(p1,p2)
c--- end statement functions
      s56=s(p5,p6)
      s123=s3(p1,p2,p3)
      s124=s3(p1,p2,p4)

      coeffsl(:,:,:)=czip

      ans34(2,2)=za(p2,p5)*zb(p4,p6)*zab2(p2,p1,p3,p4)
     & *iza(p1,p3)*iza(p2,p3)/s56
      ans43(2,2)=za(p2,p5)*zb(p3,p6)*zab2(p2,p1,p4,p3)
     & *iza(p1,p4)*iza(p2,p4)/s56

      ans34(2,1)=za(p2,p4)*za(p4,p5)*zab2(p2,p1,p3,p6)
     & *iza(p1,p3)*iza(p2,p3)/s56
      ans43(2,1)=zb(p1,p3)*zb(p3,p6)*zab2(p5,p2,p4,p1)
     & *izb(p1,p4)*izb(p2,p4)/s56

      ans34(1,2)=zb(p1,p4)*zb(p4,p6)*zab2(p5,p2,p3,p1)
     & *izb(p1,p3)*izb(p2,p3)/s56
      ans43(1,2)=za(p2,p3)*za(p3,p5)*zab2(p2,p1,p4,p6)
     & *iza(p1,p4)*iza(p2,p4)/s56

      ans34(1,1)=za(p4,p5)*zb(p1,p6)*zab2(p4,p2,p3,p1)
     & *izb(p1,p3)*izb(p2,p3)/s56
      ans43(1,1)=za(p3,p5)*zb(p1,p6)*zab2(p3,p2,p4,p1)
     & *izb(p1,p4)*izb(p2,p4)/s56


      do j=1,2
      do k=1,2
      coeffsl(c4_123,j,k)=2d0*mtsq/(s56-s123)*ans34(j,k)
      coeffsl(c3_124,j,k)=2d0*mtsq/(s56-s124)*ans43(j,k)
      coeffsl(b123,j,k)=-1d0/(s56-s123)*(1d0+s123/(s56-s123))*ans34(j,k)
      coeffsl(b124,j,k)=-1d0/(s56-s124)*(1d0+s124/(s56-s124))*ans43(j,k)
      coeffsl(b56,j,k)=-coeffsl(b123,j,k)-coeffsl(b124,j,k)
      coeffsl(rat,j,k)=1d0/(s56-s123)*ans34(j,k)+1d0/(s56-s124)*ans43(j,k)
      enddo
      enddo
cId,a34pp=1/[s56-s123]*(2+4*mtsq*C0(p123,p4,mt,mt,mt)+2*(1+s123/[s56-s123])*(B0(p1234,mt,mt)-B0(p123,mt,mt)));
cId,a43pp=1/[s56-s124]*(2+4*mtsq*C0(p3,p124,mt,mt,mt)+2*(1+s124/[s56-s124])*(B0(p1234,mt,mt)-B0(p124,mt,mt)));


      return
      end
