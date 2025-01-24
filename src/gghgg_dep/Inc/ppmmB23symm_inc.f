!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^- 4_g^- +H
      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      real(dp):: s14,s23,s123,s124,s134,s234,s1234,Delta
      complex(dp):: ppmmB23symm_res,zab2

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      s14=s(p1,p4)
      s23=s(p2,p3)
      s123=s(p1,p2)+s(p1,p3)+s(p2,p3)
      s124=s(p1,p2)+s(p1,p4)+s(p2,p4)
      s134=s(p1,p3)+s(p1,p4)+s(p3,p4)
      s234=s(p2,p3)+s(p2,p4)+s(p3,p4)
      s1234=s(p1,p2)+s(p1,p3)+s(p1,p4)+s(p2,p3)+s(p2,p4)+s(p3,p4)

      Delta=(s1234-s23-s14)**2-4._dp*s23*s14

      ppmmB23symm_res =
     & -4._dp*zb(p2,p4)**2*za(p3,p4)*s234
     & /(zb(p3,p4)*zab2(p4,p2,p3,p4)*zab2(p1,p2,p3,p4)**2)
     & -4._dp*zb(p1,p3)**2*za(p2,p3)*zab2(p1,p2,p3,p1)*zab2(p3,p4,p1,p2)
     & *(2._dp*s124+s134-3._dp*s14)
     &  /(zb(p1,p4)*zab2(p1,p2,p3,p4)*zab2(p2,p1,p4,p3)**2*Delta)
     &  +12._dp*za(p1,p4)*zb(p1,p3)**2*za(p2,p3)**2
     & *zb(p2,p3)*zab2(p3,p4,p1,p2)
     &  /(zab2(p1,p2,p3,p4)*zab2(p2,p1,p4,p3)**2*Delta)
     & +4._dp*za(p2,p4)*zb(p2,p4)*zab2(p3,p4,p1,p2)
     & /(zab2(p1,p2,p3,p4)**2*zab2(p2,p1,p4,p3))
     & -8._dp*zb(p2,p4)*za(p2,p3)*zab2(p3,p4,p1,p2)*(
     &     za(p1,p4)*zb(p3,p4)*zab2(p4,p2,p3,p1)
     &    +za(p1,p4)*za(p2,p3)*zb(p1,p3)*zb(p2,p3)
     &    +za(p2,p4)*zb(p2,p3)*zab2(p4,p2,p3,p4))
     &  /(zab2(p1,p2,p3,p4)**2*zab2(p2,p1,p4,p3)*Delta)
     & -4._dp*za(p2,p4)*zb(p1,p3)*zab2(p3,p4,p1,p2)
     &  /(zab2(p1,p2,p3,p4)*zab2(p2,p1,p4,p3)**2)
     & +8._dp*za(p1,p4)*zb(p3,p4)*zb(p1,p3)*za(p2,p3)
     & *zab2(p4,p2,p3,p1)*zab2(p3,p4,p1,p2)
     &  /(zab2(p1,p2,p3,p4)*zab2(p2,p1,p4,p3)**2*Delta)
     &     +3._dp*zab2(p4,p2,p3,p1)*zab2(p3,p4,p1,p2)
     &  *(s124-s134)*(s123-s234)*(s234+s123)
     &  /(zab2(p1,p2,p3,p4)*zab2(p2,p1,p4,p3)*Delta**2)
     & +zab2(p4,p2,p3,p1)*zab2(p3,p4,p1,p2)
     & *(2._dp*s23+5._dp*s(p2,p4)+3._dp*s(p3,p4)
     &  +3._dp*s(p1,p2)+5._dp*s(p1,p3))
     &  /(zab2(p1,p2,p3,p4)*zab2(p2,p1,p4,p3)*Delta)

      return

