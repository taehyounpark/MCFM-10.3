!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a8WZ(k1,k2,k3,k4,k5,k6,k7,k8,a)
      implicit none
c---- Matrix element for WZ radiation from line 17
c---- including W Z interchange
c  d(-k1) +c(-k2) --> e^-(k3)+ve~^+(k4)+mu^-(k5)+mu^+(k6)+u(k7)+c(k8)

c                                         5-----<-- 6 3-----<--4
c                                                \      /
c                                               Z \    / W
c                                                  \  /
c        5-----<-- 6   3-----<--4                   \/
c            |Z            |W                        |W
c   7 ----<--|-------------|----1      7 ----<-------|-------------1
c                 0                               0
c                 0                               0
c                 0                               0
c   8 -----<--------------------2      8 -----<--------------------2
c              jtype=1                         jtype=3

c        3-----<-- 4   5-----<--6
c            |W            |Z
c   7 ----<--|-------------|----1
c                 0
c                 0
c                 0
c   8 -----<--------------------2
c                jtype=2


c Overall factor of  i_*2*gs^2*e^2*gw^2*(T^A)_{i7,i1}*(T^A)_{i8,i2}
c     Tr{T^A T^B} = delta^{AB}
c     including QED type propagators
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer::k1,k2,k3,k4,k5,k6,k7,k8
      real(dp):: s3,s567,s347,s278,s128,s134,s156,s3456,
     & s34,s56,s28
      complex(dp):: zab2,a(3,2,2)
c     amplitude a(jtype,h28,h56)

      s3(k1,k2,k3)=s(k1,k2)+s(k2,k3)+s(k1,k3)
      zab2(k1,k2,k3,k4)=+za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)

      s567=s3(k5,k6,k7)
      s347=s3(k3,k4,k7)
      s278=s3(k2,k7,k8)
      s128=s3(k1,k2,k8)
      s134=s3(k1,k3,k4)
      s156=s3(k1,k5,k6)
      s3456=s(k3,k4)+s(k3,k5)+s(k3,k6)+s(k4,k5)+s(k4,k6)+s(k5,k6)
      s34=s(k3,k4)
      s56=s(k5,k6)
      s28=s(k2,k8)

      a(1,1,1)= + s34**(-1)*s56**(-1)*s28**(-1) * (  - za(k3,k7)*zb(k1,
     &    k2)*zab2(k5,k3,k7,k4)*zab2(k8,k1,k2,k6)*s347**(-1)*s128**(-1)
     &     - za(k3,k7)*zb(k1,k6)*zab2(k5,k1,k6,k2)*zab2(k8,k3,k7,k4)*
     &    s347**(-1)*s156**(-1) + za(k7,k8)*zb(k1,k6)*zab2(k3,k7,k8,k2)
     &    *zab2(k5,k1,k6,k4)*s278**(-1)*s156**(-1) )
      a(1,1,2)= + s34**(-1)*s56**(-1)*s28**(-1) * (  - za(k3,k7)*zb(k1,
     &    k2)*zab2(k6,k3,k7,k4)*zab2(k8,k1,k2,k5)*s347**(-1)*s128**(-1)
     &     - za(k3,k7)*zb(k1,k5)*zab2(k6,k1,k5,k2)*zab2(k8,k3,k7,k4)*
     &    s347**(-1)*s156**(-1) + za(k7,k8)*zb(k1,k5)*zab2(k3,k7,k8,k2)
     &    *zab2(k6,k1,k5,k4)*s278**(-1)*s156**(-1) )
      a(1,2,1)= + s34**(-1)*s56**(-1)*s28**(-1) * (  - za(k2,k7)*zb(k1,
     &    k6)*zab2(k3,k2,k7,k8)*zab2(k5,k1,k6,k4)*s278**(-1)*s156**(-1)
     &     - za(k3,k7)*zb(k1,k6)*zab2(k2,k3,k7,k4)*zab2(k5,k1,k6,k8)*
     &    s347**(-1)*s156**(-1) - za(k3,k7)*zb(k1,k8)*zab2(k2,k1,k8,k6)
     &    *zab2(k5,k3,k7,k4)*s347**(-1)*s128**(-1) )
      a(1,2,2)= + s34**(-1)*s56**(-1)*s28**(-1) * (  - za(k2,k7)*zb(k1,
     &    k5)*zab2(k3,k2,k7,k8)*zab2(k6,k1,k5,k4)*s278**(-1)*s156**(-1)
     &     - za(k3,k7)*zb(k1,k5)*zab2(k2,k3,k7,k4)*zab2(k6,k1,k5,k8)*
     &    s347**(-1)*s156**(-1) - za(k3,k7)*zb(k1,k8)*zab2(k2,k1,k8,k5)
     &    *zab2(k6,k3,k7,k4)*s347**(-1)*s128**(-1) )

      a(2,1,1)= + s34**(-1)*s56**(-1)*s28**(-1) * (  - za(k5,k7)*zb(k1,
     &    k2)*zab2(k3,k5,k7,k6)*zab2(k8,k1,k2,k4)*s128**(-1)*s567**(-1)
     &     - za(k5,k7)*zb(k1,k4)*zab2(k3,k1,k4,k2)*zab2(k8,k5,k7,k6)*
     &    s567**(-1)*s134**(-1) + za(k7,k8)*zb(k1,k4)*zab2(k3,k1,k4,k6)
     &    *zab2(k5,k7,k8,k2)*s278**(-1)*s134**(-1) )
      a(2,1,2)= + s34**(-1)*s56**(-1)*s28**(-1) * (  - za(k6,k7)*zb(k1,
     &    k2)*zab2(k3,k6,k7,k5)*zab2(k8,k1,k2,k4)*s128**(-1)*s567**(-1)
     &     - za(k6,k7)*zb(k1,k4)*zab2(k3,k1,k4,k2)*zab2(k8,k6,k7,k5)*
     &    s567**(-1)*s134**(-1) + za(k7,k8)*zb(k1,k4)*zab2(k3,k1,k4,k5)
     &    *zab2(k6,k7,k8,k2)*s278**(-1)*s134**(-1) )
      a(2,2,1)= + s34**(-1)*s56**(-1)*s28**(-1) * (  - za(k2,k7)*zb(k1,
     &    k4)*zab2(k3,k1,k4,k6)*zab2(k5,k2,k7,k8)*s278**(-1)*s134**(-1)
     &     - za(k5,k7)*zb(k1,k4)*zab2(k2,k5,k7,k6)*zab2(k3,k1,k4,k8)*
     &    s567**(-1)*s134**(-1) - za(k5,k7)*zb(k1,k8)*zab2(k2,k1,k8,k4)
     &    *zab2(k3,k5,k7,k6)*s128**(-1)*s567**(-1) )
      a(2,2,2)= + s34**(-1)*s56**(-1)*s28**(-1) * (  - za(k2,k7)*zb(k1,
     &    k4)*zab2(k3,k1,k4,k5)*zab2(k6,k2,k7,k8)*s278**(-1)*s134**(-1)
     &     - za(k6,k7)*zb(k1,k4)*zab2(k2,k6,k7,k5)*zab2(k3,k1,k4,k8)*
     &    s567**(-1)*s134**(-1) - za(k6,k7)*zb(k1,k8)*zab2(k2,k1,k8,k4)
     &    *zab2(k3,k6,k7,k5)*s128**(-1)*s567**(-1) )

      a(3,1,1)= + s34**(-1)*s56**(-1)*s28**(-1)*s3456**(-1) * ( za(k3,
     &    k5)*za(k3,k7)*zb(k1,k2)*zb(k4,k6)*zab2(k8,k1,k2,k3)*
     &    s128**(-1) + za(k3,k5)*za(k4,k7)*zb(k1,k2)*zb(k4,k6)*zab2(k8,
     &    k1,k2,k4)*s128**(-1) + za(k3,k5)*za(k7,k8)*zb(k1,k3)*zb(k4,k6
     &    )*zab2(k3,k7,k8,k2)*s278**(-1) + za(k3,k5)*za(k7,k8)*zb(k1,k4
     &    )*zb(k4,k6)*zab2(k4,k7,k8,k2)*s278**(-1) + za(k3,k7)*zb(k1,k2
     &    )*zab2(k5,k3,k4,k6)*zab2(k8,k1,k2,k4)*s128**(-1) - za(k5,k7)*
     &    zb(k1,k2)*zab2(k3,k5,k6,k4)*zab2(k8,k1,k2,k6)*s128**(-1) +
     &    za(k7,k8)*zb(k1,k4)*zab2(k3,k7,k8,k2)*zab2(k5,k3,k4,k6)*
     &    s278**(-1) - za(k7,k8)*zb(k1,k6)*zab2(k3,k5,k6,k4)*zab2(k5,k7
     &    ,k8,k2)*s278**(-1) )
      a(3,1,2)= + s34**(-1)*s56**(-1)*s28**(-1)*s3456**(-1) * (  - za(
     &    k3,k6)*za(k3,k7)*zb(k1,k2)*zb(k3,k5)*zab2(k8,k1,k2,k4)*
     &    s128**(-1) + za(k3,k6)*za(k3,k7)*zb(k1,k2)*zb(k4,k5)*zab2(k8,
     &    k1,k2,k3)*s128**(-1) + za(k3,k6)*za(k4,k7)*zb(k1,k2)*zb(k4,k5
     &    )*zab2(k8,k1,k2,k4)*s128**(-1) + za(k3,k6)*za(k7,k8)*zb(k1,k3
     &    )*zb(k4,k5)*zab2(k3,k7,k8,k2)*s278**(-1) - za(k3,k6)*za(k7,k8
     &    )*zb(k1,k4)*zb(k3,k5)*zab2(k3,k7,k8,k2)*s278**(-1) + za(k3,k6
     &    )*za(k7,k8)*zb(k1,k4)*zb(k4,k5)*zab2(k4,k7,k8,k2)*s278**(-1)
     &     - za(k3,k7)*za(k4,k6)*zb(k1,k2)*zb(k4,k5)*zab2(k8,k1,k2,k4)*
     &    s128**(-1) - za(k4,k6)*za(k7,k8)*zb(k1,k4)*zb(k4,k5)*zab2(k3,
     &    k7,k8,k2)*s278**(-1) - za(k6,k7)*zb(k1,k2)*zab2(k3,k5,k6,k4)*
     &    zab2(k8,k1,k2,k5)*s128**(-1) - za(k7,k8)*zb(k1,k5)*zab2(k3,k5
     &    ,k6,k4)*zab2(k6,k7,k8,k2)*s278**(-1) )
      a(3,2,1)= + s34**(-1)*s56**(-1)*s28**(-1)*s3456**(-1) * (  - za(
     &    k2,k7)*za(k3,k5)*zb(k1,k3)*zb(k4,k6)*zab2(k3,k2,k7,k8)*
     &    s278**(-1) - za(k2,k7)*za(k3,k5)*zb(k1,k4)*zb(k4,k6)*zab2(k4,
     &    k2,k7,k8)*s278**(-1) - za(k2,k7)*zb(k1,k4)*zab2(k3,k2,k7,k8)*
     &    zab2(k5,k3,k4,k6)*s278**(-1) + za(k2,k7)*zb(k1,k6)*zab2(k3,k5
     &    ,k6,k4)*zab2(k5,k2,k7,k8)*s278**(-1) + za(k3,k5)*za(k3,k7)*
     &    zb(k1,k8)*zb(k4,k6)*zab2(k2,k1,k8,k3)*s128**(-1) + za(k3,k5)*
     &    za(k4,k7)*zb(k1,k8)*zb(k4,k6)*zab2(k2,k1,k8,k4)*s128**(-1) +
     &    za(k3,k7)*zb(k1,k8)*zab2(k2,k1,k8,k4)*zab2(k5,k3,k4,k6)*
     &    s128**(-1) - za(k5,k7)*zb(k1,k8)*zab2(k2,k1,k8,k6)*zab2(k3,k5
     &    ,k6,k4)*s128**(-1) )
      a(3,2,2)= + s34**(-1)*s56**(-1)*s28**(-1)*s3456**(-1) * (  - za(
     &    k2,k7)*za(k3,k6)*zb(k1,k3)*zb(k4,k5)*zab2(k3,k2,k7,k8)*
     &    s278**(-1) + za(k2,k7)*za(k3,k6)*zb(k1,k4)*zb(k3,k5)*zab2(k3,
     &    k2,k7,k8)*s278**(-1) - za(k2,k7)*za(k3,k6)*zb(k1,k4)*zb(k4,k5
     &    )*zab2(k4,k2,k7,k8)*s278**(-1) + za(k2,k7)*za(k4,k6)*zb(k1,k4
     &    )*zb(k4,k5)*zab2(k3,k2,k7,k8)*s278**(-1) + za(k2,k7)*zb(k1,k5
     &    )*zab2(k3,k5,k6,k4)*zab2(k6,k2,k7,k8)*s278**(-1) - za(k3,k6)*
     &    za(k3,k7)*zb(k1,k8)*zb(k3,k5)*zab2(k2,k1,k8,k4)*s128**(-1) +
     &    za(k3,k6)*za(k3,k7)*zb(k1,k8)*zb(k4,k5)*zab2(k2,k1,k8,k3)*
     &    s128**(-1) + za(k3,k6)*za(k4,k7)*zb(k1,k8)*zb(k4,k5)*zab2(k2,
     &    k1,k8,k4)*s128**(-1) - za(k3,k7)*za(k4,k6)*zb(k1,k8)*zb(k4,k5
     &    )*zab2(k2,k1,k8,k4)*s128**(-1) - za(k6,k7)*zb(k1,k8)*zab2(k2,
     &    k1,k8,k5)*zab2(k3,k5,k6,k4)*s128**(-1) )
      return
      end
