!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a6WpgamxninvWrad_paper(p1,p2,p3,p4,p5,p6,za,zb,b22,b21)
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
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      integer::p1,p2,p3,p4,p5,p6
      real(dp):: sx,t
      complex(dp):: zab2,zab3,b22,b21,Lsm1,L0,L1,P
c      amplitude b(h5,h6)

      zab2(p1,p2,p3,p4)=+za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zab3(p1,p2,p3,p4,p5)=+za(p1,p2)*zb(p2,p5)+za(p1,p3)*zb(p3,p5)+za(p1,p4)*zb(p4,p5)
      t(p1,p2,p3)=s(p1,p2)+s(p1,p3)+s(p2,p3)
      P(sx)=cone/cmplx(sx-wmass**2,wmass*wwidth,kind=dp)

      b22=
     & P(s(p3,p4))*P(t(p3,p4,p5))*(1._dp/(2*za(p1,p6)*za(p2,p5)*zb(p2,p6)
     & *zb(p1,p2))*(za(p1,p2)*za(p3,p5)*zb(p1,p6)*(zb(p1,p5)*zb(p2,p6)+zb(p1,p6)*zb(p2,p5))*zb(p4,p5)
     & -(zb(p1,p4)*zb(p2,p6)+zb(p1,p6)*zb(p2,p4))*zab2(p2,p3,p4,p5)*zab2(p3,p1,p2,p6)
     & -za(p2,p3)*zb(p4,p5)*zb(p5,p6)*zb(p2,p6)*zab2(p5,p3,p4,p1)
     & +za(p2,p3)*zb(p4,p5)*zb(p1,p6)*zb(p2,p5)*zab2(p5,p3,p4,p6)
     & )
     & +(zb(p1,p6))/(za(p1,p6)**2*za(p2,p5))*(za(p1,p2)**2*za(p3,p6)*zb(p4,p5)*zab2(p5,p3,p4,p5)
     & +(za(p1,p6)*za(p2,p3)-za(p1,p2)*za(p3,p6))*za(p1,p3)*zb(p3,p4)*zab2(p2,p3,p4,p5)
     & )*(L0(-s(p2,p6),-t(p3,p4,p5)))/(t(p3,p4,p5))
     & +(zb(p1,p6))/(za(p1,p6)**2*za(p2,p5)*zb(p1,p2))*L0(-t(p3,p4,p5),-s(p1,p2))
     & *(-1._dp/2._dp*za(p2,p3)*zb(p4,p5)*za(p1,p6)*zab2(p5,p3,p4,p5)
     & +za(p2,p5)*za(p3,p6)*zb(p4,p5)*zab2(p1,p3,p4,p5)
     & -za(p3,p6)*zab2(p1,p2,p6,p4)*zab2(p2,p3,p4,p5))
     & -(za(p1,p2)**2*za(p3,p6))/(za(p1,p6)**3*za(p2,p6)*za(p2,p5))*(za(p2,p5)*zab2(p6,p3,p4,p5)*zb(p4,p5)
     & +zab2(p2,p3,p4,p5)*zab2(p6,p3,p5,p4)
     & )*Lsm1(-s(p1,p2),-t(p3,p4,p5),-s(p2,p6),-t(p3,p4,p5))
     & +(za(p2,p3)**2*zb(p3,p4)*zab2(p2,p3,p4,p5))/(za(p1,p6)*za(p2,p6)*za(p2,p5))
     & *Lsm1(-s(p1,p2),-t(p3,p4,p5),-s(p1,p6),-t(p3,p4,p5))
     & +(za(p1,p2)*zb(p4,p5)*zab2(p5,p3,p4,p5)
     & -za(p1,p3)*zb(p3,p4)*zab2(p2,p3,p4,p5))/(2*za(p1,p6)*za(p2,p6)*za(p2,p5)*zb(p2,p6)**2)*za(p1,p3)*zb(p1,p6)**2
     & *L1(-t(p3,p4,p5),-s(p2,p6))
     & +(zb(p1,p6))/(za(p1,p6)*za(p2,p5)*za(p1,p2)*zb(p1,p2)**2)*L1(-t(p3,p4,p5),-s(p1,p2))*((zab2(p2,p3,p4,p5)*za(p3,p6)*zb(p4,p6)
     & +za(p2,p6)*za(p3,p5)*zb(p4,p5)*zb(p5,p6)
     & +1._dp/2._dp*za(p2,p3)*zb(p4,p5)*(t(p3,p4,p6)-s(p5,p6)))*t(p3,p4,p5)
     & -1._dp/2._dp*za(p2,p3)*zb(p4,p5)*s(p1,p2)*s(p3,p4)
     & )
     & )

      b21=
     & P(s(p3,p4))*P(t(p3,p4,p5))*(L0(-t(p3,p4,p5),-s(p1,p6))*(za(p2,p6)*zab2(p3,p4,p5,p2))
     & /(za(p1,p6)*za(p2,p5)*zb(p2,p6)*zb(p2,p6))*(zab2(p2,p1,p6,p2)*zb(p4,p5)
     & +za(p2,p3)*zb(p3,p4)*zb(p2,p5))
     & -L0(-t(p3,p4,p5),-s(p1,p6))*(za(p2,p6))
     & /(za(p1,p6)*za(p2,p5)*zb(p2,p6)*zb(p1,p6))*(-zb(p1,p2)*za(p2,p3)*zb(p4,p5)*(s(p3,p5)+s(p4,p5))
     & -2*za(p2,p6)*zb(p1,p6)*zb(p4,p5)*zab2(p3,p4,p5,p2)
     & +2*za(p2,p3)*zb(p3,p4)*zb(p1,p5)*zab2(p3,p4,p5,p2)
     & )
     & -L1(-t(p3,p4,p5),-s(p1,p6))*(za(p2,p6)**2*zab2(p3,p4,p5,p2))
     & /(2*za(p1,p6)**2*zb(p1,p6)*za(p2,p5)*zb(p2,p6))*(za(p2,p3)*zb(p3,p4)*zb(p2,p5)-zb(p4,p5)*zab3(p2,p3,p4,p5,p2))
     & +(L1(-t(p3,p4,p5),-s(p1,p2)))
     & /(s(p1,p2))*(za(p2,p6)*za(p3,p6)*(za(p2,p5)*zb(p4,p5)*zb(p5,p6)+zab2(p2,p3,p4,p5)*zb(p4,p6))*t(p3,p4,p5))
     & /(za(p1,p2)*za(p2,p5)*zb(p2,p6))
     & +L0(-t(p3,p4,p5),-s(p1,p2))*(za(p2,p6))
     & /(za(p1,p2)*za(p2,p5)*zb(p2,p6)**2)*(zab2(p3,p1,p6,p2)*(za(p1,p2)*zb(p1,p6)*zb(p4,p5)+za(p2,p3)*zb(p3,p4)*zb(p5,p6))
     & +za(p2,p3)*zb(p2,p6)*zb(p4,p5)*(t(p3,p4,p5)-s(p3,p4))
     & )
     & -(zab2(p3,p2,p6,p1)*(za(p2,p3)*zb(p1,p5)*zb(p3,p4)-za(p2,p6)*zb(p1,p6)*zb(p4,p5)))/(za(p2,p5)*zb(p2,p6)*zb(p1,p6))
     & *Lsm1(-s(p1,p2),-t(p3,p4,p5),-s(p2,p6),-t(p3,p4,p5))
     & -
     & (zb(p1,p2)**2)/(za(p2,p5)*zb(p2,p6)**2)*
     & ((za(p1,p2)*zb(p4,p5)*zab2(p3,p4,p5,p6))/(zb(p2,p6))
     & +(za(p2,p3)**2*zb(p3,p4)*zb(p5,p6))/(zb(p1,p6))
     & +(za(p2,p3)*za(p1,p3)*zb(p3,p4)*zb(p5,p6))/(zb(p2,p6))
     & )*Lsm1(-s(p1,p2),-t(p3,p4,p5),-s(p1,p6),-t(p3,p4,p5))
     & +((za(p1,p6)*za(p2,p3)+za(p1,p3)*za(p2,p6))*zab2(p6,p1,p2,p5)*zb(p3,p4)*za(p2,p3)
     & -za(p2,p6)*zb(p4,p5)*(2*za(p1,p6)*za(p2,p3)*s(p3,p4)
     & +za(p1,p2)*za(p3,p6)*t(p1,p2,p6)
     & ))/(2*zb(p2,p6)*za(p2,p5)*za(p1,p6)*za(p1,p2))
     & )


      return
      end
