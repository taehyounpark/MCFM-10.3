!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a8new(k1,k2,k3,k4,k5,k6,k7,k8,a,b)
      implicit none
c---- Matrix element for W^+ W^- radiation from line 17
c     u(-k1) +c(-k2) --> vm(k3)+mu+(k4)+e^-(k5)+ve~(k6)+u(k7)+c(k8)


c           5-----<-- 6   3-----<--4
c               |W-           |W+
c      7 ----<--|-------------|----1
c                    0
c                    0                       ETC....
c                    0
c      8 -----<-------------------2


      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      include 'masses.f'
      integer k1,k2,k3,k4,k5,k6,k7,k8
      real(dp):: s34,s56,s28,s3456,s278,s134,s567,s128
      complex(dp):: zab2,zba2,zaba22,zbab22,prop34,prop56,a(2,2),b(2,2)

      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
      zba2(k1,k2,k3,k4)=zb(k1,k2)*za(k2,k4)+zb(k1,k3)*za(k3,k4)
      zbab22(k1,k2,k3,k4,k5,k6)=
     & +zba2(k1,k2,k3,k4)*zb(k4,k6)+zba2(k1,k2,k3,k5)*zb(k5,k6)
      zaba22(k1,k2,k3,k4,k5,k6)=
     & +zab2(k1,k2,k3,k4)*za(k4,k6)+zab2(k1,k2,k3,k5)*za(k5,k6)

      s34=s(k3,k4)
      s56=s(k5,k6)
      s28=s(k2,k8)
      s3456=s(k3,k4)+s(k3,k5)+s(k3,k6)+s(k4,k5)+s(k4,k6)+s(k5,k6)
      s278=s(k2,k7)+s(k7,k8)+s(k8,k2)
      s134=s(k1,k3)+s(k3,k4)+s(k4,k1)
      s567=s(k5,k6)+s(k6,k7)+s(k7,k5)
      s128=s(k1,k2)+s(k2,k8)+s(k8,k1)
      prop34=cmplx(s34,kind=dp)/cmplx(s34-wmass**2,wmass*wwidth,kind=dp)
      prop56=cmplx(s56,kind=dp)/cmplx(s56-wmass**2,wmass*wwidth,kind=dp)

c--- the 3 left-handed diagrams of type a
      a(1,1)=-za(k7,k8)*zba2(k2,k7,k8,k5)*zba2(k6,k1,k4,k3)*zb(k4,k1)
     &  /(s28*s278*s134*s34*s56)*prop34*prop56
      a(1,1)=a(1,1)
     & -za(k7,k5)*zba2(k6,k5,k7,k3)*zba2(k4,k1,k2,k8)*zb(k2,k1)
     &  /(s28*s567*s128*s34*s56)*prop34*prop56
      a(1,1)=a(1,1)
     & -za(k7,k5)*zba2(k6,k5,k7,k8)*zba2(k2,k1,k4,k3)*zb(k4,k1)
     &  /(s28*s567*s134*s34*s56)*prop34*prop56

      a(1,2)=-za(k7,k2)*zba2(k8,k7,k2,k5)*zba2(k6,k1,k4,k3)*zb(k4,k1)
     &  /(s28*s278*s134*s34*s56)*prop34*prop56
      a(1,2)=a(1,2)
     & -za(k7,k5)*zba2(k6,k5,k7,k3)*zba2(k4,k1,k8,k2)*zb(k8,k1)
     &  /(s28*s567*s128*s34*s56)*prop34*prop56
      a(1,2)=a(1,2)
     & -za(k7,k5)*zba2(k6,k5,k7,k2)*zba2(k8,k1,k4,k3)*zb(k4,k1)
     &  /(s28*s567*s134*s34*s56)*prop34*prop56
      a(2,1)=czip
      a(2,2)=czip

c---The three boson vertex diagrams
      b(2,2)=-zb(k7,k8)/(s28*s278*s34*s56*s3456)*prop34*prop56
     & *(+zaba22(k2,k7,k8,k5,k6,k1)*za(k5,k3)*zb(k4,k6)
     & +zab2(k2,k7,k8,k4)*za(k3,k1)*zab2(k5,k3,k4,k6)
     & -zab2(k2,k7,k8,k6)*za(k5,k1)*zab2(k3,k5,k6,k4))

      b(2,2)=b(2,2)+za(k2,k1)/(s28*s128*s34*s56*s3456)*prop34*prop56
     & *(+zbab22(k7,k5,k6,k1,k2,k8)*za(k5,k3)*zb(k4,k6)
     & +zb(k7,k4)*zab2(k3,k1,k2,k8)*zab2(k5,k3,k4,k6)
     & -zb(k7,k6)*zab2(k5,k1,k2,k8)*zab2(k3,k5,k6,k4))



      b(1,1)=-za(k7,k8)/(s28*s278*s34*s56*s3456)*prop34*prop56
     & *(+zbab22(k2,k7,k8,k5,k6,k1)*zb(k4,k6)*za(k5,k3)
     & +zba2(k2,k7,k8,k3)*zb(k4,k1)*zab2(k5,k3,k4,k6)
     & -zba2(k2,k7,k8,k5)*zb(k6,k1)*zab2(k3,k5,k6,k4))
      b(1,1)=b(1,1)+zb(k2,k1)/(s28*s128*s34*s56*s3456)*prop34*prop56
     & *(+zaba22(k7,k5,k6,k1,k2,k8)*za(k5,k3)*zb(k4,k6)
     & +za(k7,k3)*zba2(k4,k1,k2,k8)*zab2(k5,k3,k4,k6)
     & -za(k7,k5)*zba2(k6,k1,k2,k8)*zab2(k3,k5,k6,k4))


      b(2,1)=-zb(k7,k2)/(s28*s278*s34*s56*s3456)*prop34*prop56
     & *(+zaba22(k8,k7,k2,k5,k6,k1)*za(k5,k3)*zb(k4,k6)
     & +zab2(k8,k7,k2,k4)*za(k3,k1)*zab2(k5,k3,k4,k6)
     & -zab2(k8,k7,k2,k6)*za(k5,k1)*zab2(k3,k5,k6,k4))

      b(2,1)=b(2,1)+za(k8,k1)/(s28*s128*s34*s56*s3456)*prop34*prop56
     & *(+zbab22(k7,k5,k6,k1,k8,k2)*za(k5,k3)*zb(k4,k6)
     & +zb(k7,k4)*zab2(k3,k1,k8,k2)*zab2(k5,k3,k4,k6)
     & -zb(k7,k6)*zab2(k5,k1,k8,k2)*zab2(k3,k5,k6,k4))

      b(1,2)=-za(k7,k2)/(s28*s278*s34*s56*s3456)*prop34*prop56
     & *(+zbab22(k8,k7,k2,k5,k6,k1)*za(k5,k3)*zb(k4,k6)
     & +zba2(k8,k7,k2,k3)*zb(k4,k1)*zab2(k5,k3,k4,k6)
     & -zba2(k8,k7,k2,k5)*zb(k6,k1)*zab2(k3,k5,k6,k4))
      b(1,2)=b(1,2)+zb(k8,k1)/(s28*s128*s34*s56*s3456)*prop34*prop56
     & *(+zaba22(k7,k5,k6,k1,k8,k2)*za(k5,k3)*zb(k4,k6)
     & +za(k7,k3)*zba2(k4,k1,k8,k2)*zab2(k5,k3,k4,k6)
     & -za(k7,k5)*zba2(k6,k1,k8,k2)*zab2(k3,k5,k6,k4))


      return
      end
