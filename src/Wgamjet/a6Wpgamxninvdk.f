!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a6Wpgamxninvdk(p1,p2,p3,p4,p5,p6,za,zb,a)
      implicit none
c---- Matrix element for Wgamma radiation from W decay
c  u(-p1) +dbar(-p2) --> ve(p3)+e^+(p4)+gam(p5)+g(p6)

c               5     3-----<--4
c               \      /
c           gam \    / W
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
      real(dp):: s3,s12,s16,s26,s34,s345,Qdiff
      complex(dp):: iza,izb,zab2,zaa22,a(2,2),
     & Lsm1,L0,L1,lnrat,prop34,prop345,lo(2,2),VSL
c      amplitude b(h5,h6)

      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p1,p3)
      zab2(p1,p2,p3,p4)=+za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zaa22(p1,p2,p3,p4,p5,p6)=zab2(p1,p2,p3,p4)*za(p4,p6)
     &                        +zab2(p1,p2,p3,p5)*za(p5,p6)
      iza(p1,p2)=cone/za(p1,p2)
      izb(p1,p2)=cone/zb(p1,p2)
c  Vxninvcc(s(p1,p2),s3(p1,p2,p6))=
c -(epinv**2+epinv*lnrat(musq,-s12)+1/2*lnrat(musq,-s12)^2)
c -2*(epinv+lnrat(musq,-s345))-4
c Vxninvsc(s3(p1,p2,p6))=1/2*(epinv+lnrat(musq,-s345))+1/2
      VSL(s12,s345)=-(epinv*epinv2+epinv*lnrat(musq,-s12)
     & +0.5_dp*lnrat(musq,-s12)**2)
     & -2._dp*(epinv+lnrat(musq,-s345))-4._dp
     & +0.5_dp*(epinv+lnrat(musq,-s345))+0.5_dp

      Qdiff=Q(2)-Q(1)
      s12=s(p1,p2)
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
c      call a6WpgamxninvWrad(p1,p2,p3,p4,p5,p6,za,zb,aold(2,2),aold(2,1))
c      call a6WpgamxninvWrad(p2,p1,p4,p3,p5,p6,zb,za,aold(1,1),aold(1,2))

      call a6WpgamxninvWrad_paper(p1,p2,p3,p4,p5,p6,za,zb,a(2,2),a(2,1))
      call a6WpgamxninvWrad_paper(p2,p1,p4,p3,p5,p6,zb,za,a(1,1),a(1,2))

c      write(6,*) 'new/old a(2,2)',a(2,2)/aold(2,2)
c      write(6,*) 'new/old a(2,1)',a(2,1)/aold(2,1)
c      write(6,*) 'new/old a(1,2)',a(1,2)/aold(1,2)
c      write(6,*) 'new/old a(1,1)',a(1,1)/aold(1,1)
c      pause

c--------------------------------
      a(2,2)=a(2,2)+Qdiff/prop345*iza(p1,p6)*iza(p2,p5)*iza(p4,p5)*(
     & +L0(-s26,-s345)/s345*za(p2,p3)*zb(p1,p6)
     & *(2*za(p1,p3)*zab2(p2,p4,p5,p3)+za(p1,p2)*s(p4,p5))

     & -L0(-s345,-s26)/(za(p1,p6)*zb(p2,p6))
     & *za(p1,p3)*zb(p1,p6)*zaa22(p1,p2,p6,p4,p5,p2)

     & +L0(-s345,-s12)/(za(p1,p6)*zb(p1,p2))
     & *zb(p1,p6)*(0.5_dp*za(p1,p6)*za(p2,p3)*s(p4,p5)
     & -za(p3,p6)*zaa22(p1,p2,p6,p4,p5,p2))

     & +Lsm1(-s12,-s345,-s16,-s345)/za(p2,p6)
     & *za(p2,p3)**2*zab2(p2,p4,p5,p3)

     & -Lsm1(-s12,-s345,-s26,-s345)/(za(p1,p6)**2*za(p2,p6))
     & *za(p1,p2)**2*za(p3,p6)*zaa22(p2,p4,p5,p1,p2,p6)

     & -0.5_dp*L1(-s345,-s26)/(s(p2,p6)*zb(p2,p6))
     & *za(p1,p3)*zb(p1,p6)**2*zaa22(p1,p2,p6,p4,p5,p2)

     & +0.5_dp*L1(-s345,-s12)/(s(p1,p2)*za(p2,p6)*zb(p1,p2))
     & *zb(p1,p6)*(2*za(p2,p6)*zab2(p2,p4,p5,p6)*za(p3,p6)*s3(p1,p2,p6)
     & +za(p2,p3)*za(p2,p6)*zb(p4,p5)*za(p4,p5)*(s(p1,p6)+s(p2,p6)))

     & -0.5_dp/(zb(p2,p6)*zb(p1,p2))
     & *(-za(p2,p3)*s(p4,p5)*zb(p2,p6)*zb(p1,p6)
     & -zab2(p3,p1,p2,p6)*
     & (za(p2,p4)*(zb(p2,p6)*zb(p1,p4)+zb(p1,p6)*zb(p2,p4))
     & +za(p2,p5)*(zb(p2,p6)*zb(p1,p5)+zb(p1,p6)*zb(p2,p5)))))

c--------------------------------
      a(2,1)=a(2,1)+Qdiff/prop345*iza(p2,p5)*iza(p4,p5)*izb(p2,p6)*(
     &-L0(-s345,-s16)*iza(p1,p6)*izb(p2,p6)
     &*(s(p2,p5)+s(p2,p4))*za(p2,p6)*zab2(p3,p1,p6,p2)

     &+0.5_dp*L1(-s345,-s16)/s(p1,p6)*iza(p1,p6)
     &*(s(p2,p5)+s(p2,p4))*za(p2,p6)**2*zab2(p3,p4,p5,p2)

     &+L0(-s16,-s345)/s345*za(p2,p6)
     &*(zb(p1,p2)*zb(p4,p5)*za(p2,p3)*za(p4,p5)
     &+2*zab2(p3,p2,p6,p1)*(za(p2,p4)*zb(p2,p4)+zb(p2,p5)*za(p2,p5)))

     &-0.5_dp*L1(-s345,-s12)/s(p1,p2)*iza(p1,p2)
     &*za(p2,p6)*(2*s3(p1,p2,p6)*zab2(p2,p4,p5,p6)*za(p3,p6)
     &-(s(p2,p6)+s(p1,p6))*za(p2,p3)*s(p4,p5))

     &+0.5_dp*L0(-s345,-s12)*iza(p1,p2)*izb(p2,p6)*za(p2,p6)
     &*(2*zab2(p3,p4,p5,p2)*zab2(p2,p4,p5,p6)
     &-za(p2,p3)*s(p4,p5)*zb(p2,p6))

     &-Lsm1(-s12,-s345,-s26,-s345)*izb(p1,p6)
     &*zab2(p2,p4,p5,p1)*zab2(p3,p2,p6,p1)

     &+Lsm1(-s12,-s345,-s16,-s345)*izb(p2,p6)**2*izb(p1,p6)
     &*zb(p1,p2)*(zab2(p3,p2,p6,p1)*zab2(p2,p4,p5,p6)*zb(p2,p6)
     &-za(p1,p3)*zb(p1,p6)**2*s(p4,p5)
     &-zab2(p3,p1,p2,p6)*zb(p1,p6)*s3(p2,p4,p5))

     &+0.5_dp*iza(p1,p2)*iza(p1,p6)
     &*(-(za(p1,p6)*za(p2,p3)+za(p1,p3)*za(p2,p6))
     &*za(p1,p6)*zab2(p2,p4,p5,p1)
     &-za(p2,p6)**2*za(p1,p3)*(s(p2,p4)+s(p2,p5))
     &+za(p2,p3)*za(p2,p6)*za(p1,p6)*(-s3(p2,p4,p5))))

c--------------------------------
      a(1,1)=-a(1,1)
     &+Qdiff/(prop345*zb(p1,p5)*zb(p2,p6)*zb(p4,p5))*(
     & +L0(-s345,-s16)/(s(p1,p6)*zb(p2,p6))
     & *za(p2,p6)*zb(p1,p4)*(
     & +za(p1,p3)*zb(p1,p2)**2*zb(p4,p6)
     & -za(p2,p3)*zb(p2,p4)*zb(p1,p2)*zb(p2,p6)
     & +za(p3,p6)*zb(p1,p4)*zb(p2,p6)**2)

     & +L0(-s345,-s12)/(za(p1,p2)*zb(p2,p6))*za(p2,p6)*zb(p1,p4)
     & *(0.5_dp*za(p3,p5)*zb(p4,p5)*zb(p2,p6)
     & -zab2(p3,p1,p2,p6)*zb(p2,p4))

     & -Lsm1(-s12,-s345,-s26,-s345)/zb(p1,p6)
     & *zb(p1,p4)**2*zab2(p3,p2,p6,p1)

     & +Lsm1(-s12,-s345,-s16,-s345)/(zb(p2,p6)**2*zb(p1,p6))
     & *zb(p1,p2)**2*zb(p1,p4)*zb(p4,p6)*zab2(p3,p1,p2,p6)

     & -0.5_dp*L1(-s345,-s16)/(s(p1,p6)*za(p1,p6))
     & *za(p2,p6)**2*zb(p1,p4)*zb(p2,p4)*zab2(p3,p1,p6,p2)

     & +0.5_dp*L1(-s345,-s12)/(s(p1,p2)*za(p1,p2))*zb(p1,p4)*za(p2,p6)
     & *(+zab2(p3,p1,p2,p4)*(s(p1,p6)+s(p2,p6))
     & +za(p3,p6)*zb(p4,p6)*(s(p1,p6)+s(p2,p6)+2*s(p1,p2)))

     & -0.5_dp/(za(p1,p2)*za(p1,p6))
     & *zb(p1,p4)*(za(p3,p6)*za(p2,p6)*zab2(p1,p2,p6,p4)
     & -za(p1,p6)*za(p2,p3)*zab2(p6,p1,p2,p4)))

c--------------------------------
      a(1,2)=-a(1,2)
     & +Qdiff/(prop345*za(p1,p6)*zb(p1,p5)*zb(p4,p5))*(
     & 0.5_dp*L1(-s345,-s12)/(s(p1,p2)*zb(p1,p2))*zb(p1,p6)*zb(p1,p4)
     & *(-zab2(p3,p1,p2,p4)*(s(p1,p6)+s(p2,p6))
     & -za(p3,p6)*zb(p4,p6)*(s(p1,p6)+s(p2,p6)+2*s(p1,p2)))

     & +L0(-s345,-s26)/s(p2,p6)*za(p1,p2)*zb(p1,p4)*zb(p1,p6)
     & *(+zab2(p3,p2,p6,p4)-za(p3,p1)*zb(p1,p4))

     & -0.5_dp*L1(-s345,-s26)/zb(p2,p6)/s(p2,p6)
     & *zab2(p1,p2,p6,p4)*za(p1,p3)*zb(p1,p4)*zb(p1,p6)**2

     & +L0(-s345,-s26)/(za(p1,p6)*zb(p2,p6))
     & *za(p1,p3)*zb(p1,p4)*zb(p1,p6)
     & *(za(p2,p1)*zb(p2,p4)+za(p1,p6)*zb(p6,p4))

     & +L0(-s345,-s12)/(za(p1,p6)*zb(p1,p2))*zb(p1,p4)*zb(p1,p6)
     & *(za(p3,p6)*(-za(p1,p2)*zb(p2,p4)-0.5_dp*za(p1,p6)*zb(p6,p4))
     & +0.5_dp*za(p1,p6)*zab2(p3,p1,p2,p4))

     & +Lsm1(-s12,-s345,-s26,-s345)/(za(p1,p6)**2*za(p2,p6))
     & *za(p1,p2)**2*za(p3,p6)*zb(p1,p4)*zab2(p6,p1,p2,p4)

     & -Lsm1(-s12,-s345,-s16,-s345)/za(p2,p6)
     & *za(p2,p3)*zb(p1,p4)*zab2(p2,p1,p6,p4)

     & +0.5_dp*zb(p1,p4)/(zb(p1,p2)*zb(p2,p6))
     & *(zab2(p3,p2,p6,p1)*zb(p4,p6)*zb(p2,p6)
     &  -zab2(p3,p1,p2,p6)*zb(p1,p6)*zb(p2,p4)))

c     sign to take account of sign in xninv lowest order
      a(:,:)=a(:,:)-VSL(s12,s345)*lo(:,:)

      return
      end
