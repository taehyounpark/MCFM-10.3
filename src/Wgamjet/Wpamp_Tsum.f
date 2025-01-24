!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function Wpamp_Tsum(p1,p2,p3,p4,p5,p6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6
      complex(dp):: Wpamp_Tsum,zab,zab2
      real(dp):: t,Delta3,DeltaTM1526

c Statement functions
      zab(p1,p2,p3)=za(p1,p2)*zb(p2,p3)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      Delta3(p1,p2,p3,p4,p5,p6)=
     & (s(p1,p2)**2+s(p3,p4)**2+s(p5,p6)**2
     & -2._dp*(s(p1,p2)*s(p3,p4)+s(p1,p2)*s(p5,p6)+s(p3,p4)*s(p5,p6)))

      DeltaTM1526=Delta3(p1,p5,p2,p6,p3,p4)

      Wpamp_Tsum=
     & s(p3,p4)*(
     & 2*(za(p1,p5)*zb(p1,p6)*za(p2,p3)*zab2(p5,p1,p6,p4))
     & /(za(p1,p6)*t(p2,p3,p4)**2)
     & -(za(p1,p2)*zb(p1,p5)*za(p2,p3)*zab2(p2,p1,p6,p4)*zab2(p5,p1,p6,p2))
     & /(za(p1,p6)*zab2(p2,p1,p6,p5)*zab2(p6,p1,p2,p5)*t(p2,p3,p4))
     & -2*(za(p1,p5)*zb(p1,p6)*za(p2,p3)*zb(p2,p4)*za(p5,p6))
     & /(za(p1,p6)*zab2(p6,p1,p5,p2)*t(p2,p3,p4))
     & +(za(p1,p2)*zab2(p3,p2,p5,p1)*zab2(p5,p1,p6,p4))
     & /(za(p1,p6)*zab2(p6,p1,p2,p5)*t(p2,p3,p4))
     & -(za(p1,p2)*zb(p2,p6)*za(p3,p6)*zb(p4,p5)*zab2(p5,p2,p6,p1))
     & /(zab2(p1,p2,p6,p5)*zab2(p6,p1,p5,p2)*zab2(p6,p1,p2,p5))
     & +(2*za(p2,p6)*zb(p2,p6)*zab(p1,p5,p6)*zab2(p3,p1,p5,p4)*zab2(p5,p2,p6,p1))
     & /(zab2(p1,p2,p6,p5)*zab2(p6,p1,p5,p2)*DeltaTM1526)
     & +(za(p1,p2)*zb(p2,p6)*zab2(p3,p1,p5,p4)*zab2(p5,p2,p6,p1)
     & *(s(p1,p2)+s(p2,p5)+s(p1,p6)+s(p5,p6)))
     & /(zab2(p1,p2,p6,p5)*zab2(p6,p1,p5,p2)*DeltaTM1526)
     & -(za(p1,p2)*zb(p1,p6)*zb(p2,p4)*za(p3,p6)*za(p5,p6))
     & /(za(p1,p6)*zab2(p6,p1,p5,p2)*zab2(p6,p1,p2,p5))
     & +(za(p1,p2)*zb(p1,p4)*za(p3,p5))
     & /(za(p1,p6)*zab2(p6,p1,p2,p5))
     & )

      return
      end
