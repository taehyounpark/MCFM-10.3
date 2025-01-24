!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^- 4_g^- +H
      include 'Inc/zprods_decl.f'

      integer p1,p2,p3,p4
      real(dp):: s12,s13,s14,s23,s24,s34,p14Dp23,Delta
      complex(dp):: ppmmC23x41m2_unsym_res,zab2

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      s12=s(p1,p2)
      s13=s(p1,p3)
      s14=s(p1,p4)
      s23=s(p2,p3)
      s24=s(p2,p4)
      s34=s(p3,p4)
      p14Dp23=(s12+s13+s24+s34)/2._dp
      Delta=4._dp*(p14Dp23**2-s14*s23)

      ppmmC23x41m2_unsym_res =
     & 4._dp*za(p2,p4)*zab2(p3,p1,p4,p2)
     & /(za(p2,p3)*za(p1,p4)*zab2(p2,p1,p4,p3))
     & *((za(p2,p4)*(s13-s24)/zab2(p2,p1,p4,p3)
     & -2._dp*za(p3,p4))/zab2(p1,p2,p3,p4)
     & -zab2(p4,p2,p3,p1)
     & *(za(p1,p2)*zab2(p3,p1,p4,p2)/zab2(p1,p2,p3,p4)-za(p3,p4))/Delta)

      return


