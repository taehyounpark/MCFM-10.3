!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a8ZxZ(k1,k2,k3,k4,k5,k6,k7,k8,a)
      implicit none
c---- Matrix element for Z radiation from line 17 and line 28
c---- Four diagrams, before and after on both lines
c     u(-k1) +c(-k2) --> mu^-(k3)+mu^+(k4)+e^-(k5)+e^+(k6)+u(k7)+c(k8)


c           3-----<--4
c               |Z
c      7 ----<--|-------------------1
c                    0
c                    0                       ETC....
c                    0
c      8 -----<---------------|-----2
c                             |
c                         5----<---6


c     Overall factor of  i_*gs^2*e^4*(T^A)_{i7,i1}*(T^A)_{i8,i2}
c     Tr{T^A T^B} = delta^{AB}
c     including QED type propagators
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer::k1,k2,k3,k4,k5,k6,k7,k8
      real(dp):: s3,s4
      complex(dp):: zab2,zaa22,zbb22,ansmm,ansmp,anspm,anspp
c     amplitude a(h17,h26,h34,h56)
      complex(dp):: a(2,2,2,2)
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
      zaa22(k1,k2,k3,k4,k5,k6)=
     & +zab2(k1,k2,k3,k4)*za(k4,k6)+zab2(k1,k2,k3,k5)*za(k5,k6)
      zbb22(k1,k2,k3,k4,k5,k6)=
     & +zb(k1,k2)*zab2(k2,k4,k5,k6)+zb(k1,k3)*zab2(k3,k4,k5,k6)
      s3(k1,k2,k3)=s(k1,k2)+s(k2,k3)+s(k1,k3)
      s4(k1,k2,k3,k4)=s(k1,k2)+s(k1,k3)+s(k1,k4)
     &               +s(k2,k3)+s(k2,k4)+s(k3,k4)
      ansmm(k1,k2,k3,k4,k5,k6,k7,k8)=(
     & za(k8,k7)*zb(k1,k4)*zb(k2,k6)*zaa22(k5,k2,k6,k1,k4,k3)
     & /(s3(k2,k5,k6)*s3(k1,k3,k4))
     & +zab2(k8,k3,k7,k4)*za(k3,k7)*zb(k2,k6)*zab2(k5,k2,k6,k1)
     & /(s3(k2,k5,k6)*s3(k3,k4,k7))
     & +zab2(k3,k1,k4,k2)*za(k5,k8)*zb(k1,k4)*zab2(k7,k5,k8,k6)
     & /(s3(k5,k6,k8)*s3(k1,k3,k4))
     & +za(k3,k7)*za(k5,k8)*zb(k1,k2)*zbb22(k4,k7,k3,k5,k8,k6)
     & /(s3(k5,k6,k8)*s3(k3,k4,k7)))/(s(k3,k4)*s(k5,k6)*s4(k1,k3,k4,k7))

      ansmp(k1,k2,k3,k4,k5,k6,k7,k8)=(
     & +zab2(k3,k1,k4,k8)*za(k2,k5)*zb(k1,k4)*zab2(k7,k2,k5,k6)
     & /(s3(k2,k5,k6)*s3(k1,k3,k4))
     & +za(k2,k5)*za(k3,k7)*zb(k1,k8)*zbb22(k4,k3,k7,k2,k5,k6)
     & /(s3(k2,k5,k6)*s3(k3,k4,k7))
     & +za(k2,k7)*zb(k1,k4)*zb(k6,k8)*zaa22(k5,k6,k8,k1,k4,k3)
     & /(s3(k5,k6,k8)*s3(k1,k3,k4))
     & +zab2(k2,k3,k7,k4)*za(k3,k7)*zb(k6,k8)*zab2(k5,k6,k8,k1)
     & /(s3(k5,k6,k8)*s3(k3,k4,k7)))/(s(k3,k4)*s(k5,k6)*s4(k1,k3,k4,k7))

      anspm(k1,k2,k3,k4,k5,k6,k7,k8)=(
     & za(k3,k1)*zab2(k8,k1,k3,k4)*zb(k6,k2)*zab2(k5,k2,k6,k7)
     & /(s3(k2,k5,k6)*s3(k1,k3,k4))
     & +za(k8,k1)*zb(k6,k2)*zb(k7,k4)*zaa22(k5,k2,k6,k4,k7,k3)
     & /(s3(k2,k5,k6)*s3(k3,k4,k7))
     & +za(k3,k1)*za(k8,k5)*zb(k7,k2)*zbb22(k4,k1,k3,k5,k8,k6)
     & /(s3(k5,k6,k8)*s3(k1,k3,k4))
     & +zab2(k3,k4,k7,k2)*za(k5,k8)*zb(k4,k7)*zab2(k1,k5,k8,k6)
     & /(s3(k5,k6,k8)*s3(k3,k4,k7)))/(s(k3,k4)*s(k5,k6)*s4(k1,k3,k4,k7))

      anspp(k1,k2,k3,k4,k5,k6,k7,k8)=(
     & za(k1,k3)*za(k5,k2)*zb(k8,k7)*zbb22(k4,k1,k3,k2,k5,k6)
     & /(s3(k2,k5,k6)*s3(k1,k3,k4))
     & +zab2(k3,k4,k7,k8)*za(k5,k2)*zb(k7,k4)*zab2(k1,k2,k5,k6)
     & /(s3(k2,k5,k6)*s3(k3,k4,k7))
     & +zab2(k2,k1,k3,k4)*za(k3,k1)*zb(k8,k6)*zab2(k5,k6,k8,k7)
     & /(s3(k5,k6,k8)*s3(k1,k3,k4))
     & +za(k1,k2)*zb(k4,k7)*zb(k8,k6)*zaa22(k5,k6,k8,k4,k7,k3)
     & /(s3(k5,k6,k8)*s3(k3,k4,k7)))/(s(k3,k4)*s(k5,k6)*s4(k1,k3,k4,k7))

      a(1,1,1,1)=ansmm(k1,k2,k3,k4,k5,k6,k7,k8)
      a(1,1,1,2)=ansmm(k1,k2,k3,k4,k6,k5,k7,k8)
      a(1,1,2,1)=ansmm(k1,k2,k4,k3,k5,k6,k7,k8)
      a(1,1,2,2)=ansmm(k1,k2,k4,k3,k6,k5,k7,k8)

      a(1,2,1,1)=ansmp(k1,k2,k3,k4,k5,k6,k7,k8)
      a(1,2,1,2)=ansmp(k1,k2,k3,k4,k6,k5,k7,k8)
      a(1,2,2,1)=ansmp(k1,k2,k4,k3,k5,k6,k7,k8)
      a(1,2,2,2)=ansmp(k1,k2,k4,k3,k6,k5,k7,k8)

      a(2,1,1,1)=anspm(k1,k2,k3,k4,k5,k6,k7,k8)
      a(2,1,1,2)=anspm(k1,k2,k3,k4,k6,k5,k7,k8)
      a(2,1,2,1)=anspm(k1,k2,k4,k3,k5,k6,k7,k8)
      a(2,1,2,2)=anspm(k1,k2,k4,k3,k6,k5,k7,k8)

      a(2,2,1,1)=anspp(k1,k2,k3,k4,k5,k6,k7,k8)
      a(2,2,1,2)=anspp(k1,k2,k3,k4,k6,k5,k7,k8)
      a(2,2,2,1)=anspp(k1,k2,k4,k3,k5,k6,k7,k8)
      a(2,2,2,2)=anspp(k1,k2,k4,k3,k6,k5,k7,k8)

      return
      end
