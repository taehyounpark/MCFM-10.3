!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_qb^+ 2_q^- 3_g^- 4_g^+ +H
c     Implementation of Eq.~(9.8) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'

      integer p1,p2,p3,p4
      real(dp):: s14p24p34
      complex(dp):: aqmpC4x123m2_res,zab2

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      s14p24p34=s(p1,p4)+s(p2,p4)+s(p3,p4)
      aqmpC4x123m2_res=
     & -8*zb(p1,p4)**2*zb(p2,p4)
     & /(zb(p1,p2)*zb(p2,p3)*zb(p3,p4)*s14p24p34)
     & +4*za(p2,p3)*za(p2,p4)*s14p24p34
     & /(za(p1,p2)*za(p3,p4)*zab2(p4,p1,p2,p3)**2)
     & +4*zb(p1,p3)*zb(p1,p4)*s14p24p34
     & /(zb(p1,p2)*zb(p3,p4)*zab2(p4,p1,p2,p3)**2)
      return
