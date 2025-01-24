!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a8ZZ(k1,k2,k3,k4,k5,k6,k7,k8,a)
      implicit none
c---- Matrix element for ZZ radiation from line 17
c---- Six diagrams, before,aft and straddle + (34<-->56) interchange
c     u(-k1) +c(-k2) --> mu^-(k3)+mu^+(k4)+e^-(k5)+e^+(k6)+u(k7)+c(k8)


c           5-----<-- 6   3-----<--4
c               |Z            |Z
c      7 ----<--|-------------|----1
c                    0
c                    0                       ETC....
c                    0
c      8 -----<--------------------2


c     Overall factor of  i_*gs^2*e^4*(T^A)_{i7,i1}*(T^A)_{i8,i2}
c     Tr{T^A T^B} = delta^{AB}
c     including QED type propagators
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer::k1,k2,k3,k4,k5,k6,k7,k8
      real(dp):: s3
      complex(dp):: zab2,ansm,ansp
c     amplitude a(h17,h28,h34,h56)
      complex(dp):: a(2,2,2,2)
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
      s3(k1,k2,k3)=s(k1,k2)+s(k2,k3)+s(k1,k3)
      ansm(k1,k2,k3,k4,k5,k6,k7,k8)=
     &(-za(k5,k7)*zb(k1,k2)*zab2(k3,k5,k7,k6)*zab2(k8,k1,k2,k4)
     & /(s3(k5,k6,k7)*s3(k1,k2,k8))
     & -za(k5,k7)*zb(k1,k4)*zab2(k3,k1,k4,k2)*zab2(k8,k5,k7,k6)
     & /(s3(k5,k6,k7)*s3(k1,k3,k4))
     & +za(k7,k8)*zb(k1,k4)*zab2(k3,k1,k4,k6)*zab2(k5,k7,k8,k2)
     & /(s3(k2,k7,k8)*s3(k1,k3,k4)))/(s(k3,k4)*s(k5,k6)*s(k2,k8))
      ansp(k1,k2,k3,k4,k5,k6,k7,k8)=
     &(-za(k1,k3)*zb(k2,k7)*zab2(k5,k1,k3,k4)*zab2(k8,k2,k7,k6)
     & /(s3(k2,k7,k8)*s3(k1,k3,k4))
     & -za(k1,k3)*zb(k6,k7)*zab2(k5,k6,k7,k2)*zab2(k8,k1,k3,k4)
     & /(s3(k5,k6,k7)*s3(k1,k3,k4))
     & -za(k1,k8)*zb(k6,k7)*zab2(k3,k1,k8,k2)*zab2(k5,k6,k7,k4)
     & /(s3(k5,k6,k7)*s3(k1,k2,k8)))/(s(k3,k4)*s(k5,k6)*s(k2,k8))

c    ans(h17,h28,h34,h56) including (34 <--> 56) interchange
      a(1,1,1,1)=
     & ansm(k1,k2,k3,k4,k5,k6,k7,k8)+ansm(k1,k2,k5,k6,k3,k4,k7,k8)
      a(1,1,1,2)=
     & ansm(k1,k2,k3,k4,k6,k5,k7,k8)+ansm(k1,k2,k6,k5,k3,k4,k7,k8)
      a(1,1,2,1)=
     & ansm(k1,k2,k4,k3,k5,k6,k7,k8)+ansm(k1,k2,k5,k6,k4,k3,k7,k8)
      a(1,1,2,2)=
     & ansm(k1,k2,k4,k3,k6,k5,k7,k8)+ansm(k1,k2,k6,k5,k4,k3,k7,k8)
      a(1,2,1,1)=
     & ansm(k1,k8,k3,k4,k5,k6,k7,k2)+ansm(k1,k8,k5,k6,k3,k4,k7,k2)
      a(1,2,1,2)=
     & ansm(k1,k8,k3,k4,k6,k5,k7,k2)+ansm(k1,k8,k6,k5,k3,k4,k7,k2)
      a(1,2,2,1)=
     & ansm(k1,k8,k4,k3,k5,k6,k7,k2)+ansm(k1,k8,k5,k6,k4,k3,k7,k2)
      a(1,2,2,2)=
     & ansm(k1,k8,k4,k3,k6,k5,k7,k2)+ansm(k1,k8,k6,k5,k4,k3,k7,k2)

      a(2,1,1,1)=
     & ansp(k1,k2,k3,k4,k5,k6,k7,k8)+ansp(k1,k2,k5,k6,k3,k4,k7,k8)
      a(2,1,1,2)=
     & ansp(k1,k2,k3,k4,k6,k5,k7,k8)+ansp(k1,k2,k6,k5,k3,k4,k7,k8)
      a(2,1,2,1)=
     & ansp(k1,k2,k4,k3,k5,k6,k7,k8)+ansp(k1,k2,k5,k6,k4,k3,k7,k8)
      a(2,1,2,2)=
     & ansp(k1,k2,k4,k3,k6,k5,k7,k8)+ansp(k1,k2,k6,k5,k4,k3,k7,k8)
      a(2,2,1,1)=
     & ansp(k1,k8,k3,k4,k5,k6,k7,k2)+ansp(k1,k8,k5,k6,k3,k4,k7,k2)
      a(2,2,1,2)=
     & ansp(k1,k8,k3,k4,k6,k5,k7,k2)+ansp(k1,k8,k6,k5,k3,k4,k7,k2)
      a(2,2,2,1)=
     & ansp(k1,k8,k4,k3,k5,k6,k7,k2)+ansp(k1,k8,k5,k6,k4,k3,k7,k2)
      a(2,2,2,2)=
     & ansp(k1,k8,k4,k3,k6,k5,k7,k2)+ansp(k1,k8,k6,k5,k4,k3,k7,k2)

      return
      end

