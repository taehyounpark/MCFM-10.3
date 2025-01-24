!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a6Wpgamxndk(p1,p2,p3,p4,p5,p6,za,zb,b)
      implicit none
c---- Matrix element for Wgamma radiation from W decay
c  d(-p1) +ubar(-p2) --> ve(p3)+e^+(p4)+gam(p5)+g(p6)

c               5     3-----<--4
c               \      /
c          gam \    / W
c                 \  /
c                  \/
c                     |W
c     2 ----<-------|-------------1
c                0
c                0
c                0
c             jtype=3

      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'scale.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'ewcharge.f'
      include 'masses.f'
      integer::p1,p2,p3,p4,p5,p6
      real(dp):: s3,s16,s26,s34,s345,Qdiff
      complex(dp):: iza,izb,zab2,b(2,2),
     & Lsm1,L0,L1,lnrat,prop34,prop345,lo(2,2),VLC
c      amplitude b(h6)

      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p1,p3)
      zab2(p1,p2,p3,p4)=+za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      iza(p1,p2)=cone/za(p1,p2)
      izb(p1,p2)=cone/zb(p1,p2)
c Vxncc(s(p2,p6),s(p1,p6))=
c  -(epinv**2+epinv*lnrat(musq,-s(p2,p6))+1/2*lnrat(musq,-s(p2,p6))^2)
c  -(epinv**2+epinv*lnrat(musq,-s(p1,p6))+1/2*lnrat(musq,-s(p1,p6))^2)
c  -2*(epinv+lnrat(musq,-s(p1,p6)))-4
c Vxnsc(s(p1,p6)) =1/2*(epinv+lnrat(musq,-s(p1,p6)))+1
      VLC(s26,s16)=
     & -2._dp*epinv*epinv2-epinv*(lnrat(musq,-s16)+lnrat(musq,-s26))
     & -0.5_dp*(lnrat(musq,-s16)**2+lnrat(musq,-s26)**2)
     & -2._dp*(epinv+lnrat(musq,-s16))-4
     & +0.5_dp*(epinv+lnrat(musq,-s16))+1._dp

      Qdiff=Q(2)-Q(1)
      s16=s(p1,p6)
      s26=s(p2,p6)
      s34=s(p3,p4)
      s345=s3(p3,p4,p5)
      prop34=cmplx(s34-wmass**2,wmass*wwidth,kind=dp)
      prop345=cmplx(s345-wmass**2,wmass*wwidth,kind=dp)

      lo(1,1)=Qdiff*zb(p1,p4)**2/(prop345*zb(p1,p5)*zb(p1,p6)*zb(p2,p6))
     & *(za(p3,p4)*zab2(p5,p3,p4,p1)/prop34+zab2(p3,p4,p5,p1)/zb(p4,p5))

      lo(2,2)=Qdiff*za(p2,p3)**2/(prop345*za(p2,p5)*za(p2,p6)*za(p1,p6))
     & *(zb(p3,p4)*zab2(p2,p3,p4,p5)/prop34+zab2(p2,p4,p5,p3)/za(p4,p5))

      lo(1,2)=-Qdiff*zab2(p2,p1,p6,p4)
     & /(prop345*zb(p1,p5)*za(p1,p6)*za(p2,p6))
     & *((za(p2,p5)*za(p3,p4)*zb(p1,p4)+za(p3,p5)*za(p2,p6)*zb(p1,p6))
     & /prop34+za(p2,p3)*zb(p1,p4)/zb(p4,p5))

      lo(2,1)=+Qdiff*zab2(p3,p2,p6,p1)
     & /(prop345*za(p2,p5)*zb(p2,p6)*zb(p1,p6))
     & *((zb(p1,p5)*zb(p4,p3)*za(p2,p3)+zb(p4,p5)*zb(p1,p6)*za(p2,p6))
     & /prop34-zab2(p2,p4,p5,p1)/za(p4,p5))

c     Radiation from W-line satisfies a simple rule for swapping helicities
      call a6WpgamxnWrad(p1,p2,p3,p4,p5,p6,za,zb,b(2,2),b(2,1))
      call a6WpgamxnWrad(p2,p1,p4,p3,p5,p6,zb,za,b(1,1),b(1,2))

      b(2,2)=b(2,2)+Qdiff/prop345*iza(p2,p5)*iza(p1,p6)*iza(p2,p6)*(
     & (za(p1,p2)*za(p2,p3)*iza(p4,p5)*L0(-s26,-s345)
     & *(0.5_dp*za(p2,p6)*zb(p1,p6)*s(p4,p5)/s345-zab2(p2,p4,p5,p1))
     & +za(p2,p3)**2*zab2(p2,p1,p6,p3)*iza(p4,p5)*Lsm1(-s16,-s345,-s26,-s345)
     & +0.5_dp*za(p1,p2)**2*zab2(p2,p4,p5,p1)*zab2(p3,p2,p6,p1)
     & /s345*iza(p4,p5)*L1(-s26,-s345)))

      b(2,1)=b(2,1)+Qdiff/prop345*iza(p2,p5)*izb(p1,p6)*izb(p2,p6)*(
     &  + za(p2,p3)*zb(p1,p2)*iza(p4,p5)*L0(-s16,-s345)
     &  *(zab2(p2,p4,p5,p1)
     &  -0.5_dp*za(p2,p6)*zb(p1,p6)*s(p4,p5)/s345)
     &  +zab2(p3,p2,p6,p1)*zab2(p2,p4,p5,p1)/za(p4,p5)
     &  *Lsm1(-s26,-s345,-s16,-s345)
     &  -0.5_dp*(za(p2,p3)*zb(p1,p2))**2*zab2(p2,p4,p5,p3)
     &  *L1(-s16,-s345)/(za(p4,p5)*s345))

      b(1,1)=-b(1,1)+Qdiff/prop345*izb(p1,p5)*izb(p1,p6)*izb(p2,p6)
     & *(za(p2,p3)*zb(p1,p2)*zb(p1,p4)**2*izb(p4,p5)*L0(-s16,-s345)
     & +0.5_dp*zab2(p2,p1,p6,p4)*za(p2,p3)*zb(p1,p2)**2*zb(p1,p4)
     & *izb(p4,p5)/s345*L1(-s16,-s345)
     & -0.5_dp*za(p2,p6)*za(p3,p5)*zb(p1,p2)*zb(p1,p4)*zb(p1,p6)
     & *L0(-s16,-s345)/s345
     & +zab2(p3,p2,p6,p1)*zb(p1,p4)**2*izb(p4,p5)
     & *Lsm1(-s26,-s345,-s16,-s345))

      b(1,2)=-b(1,2)-Qdiff/(prop345*zb(p1,p5)*za(p2,p6)*za(p1,p6))
     & *(za(p1,p2)*zb(p1,p4)*L0(-s26,-s345)
     & *(za(p2,p3)*zb(p1,p4)/zb(p4,p5)
     & -0.5_dp*za(p2,p6)*za(p3,p5)*zb(p1,p6)/s345)
     & -za(p2,p3)*zb(p1,p4)*zab2(p2,p1,p6,p4)/zb(p4,p5)
     & *Lsm1(-s16,-s345,-s26,-s345)
     & -0.5_dp*(za(p1,p2)*zb(p1,p4))**2*zab2(p3,p2,p6,p1)
     & *L1(-s26,-s345)/(s345*zb(p4,p5)))

      b(:,1)=b(:,1)+VLC(s26,s16)*lo(:,1)
      b(:,2)=b(:,2)+VLC(s16,s26)*lo(:,2)

      return
      end
