!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_qb^+ 2_q^- 3_g^- 4_g^+ +H
c     Implementation of Eq.~(9.7) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'

      integer p1,p2,p3,p4
      real(dp):: s12,s13,s14,s23,s24,s34,s123,s14p24p34
      complex(dp):: aqmpC4x123m0_res,zab2

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      s12=s(p1,p2)
      s13=s(p1,p3)
      s14=s(p1,p4)
      s23=s(p2,p3)
      s24=s(p2,p4)
      s34=s(p3,p4)

      s123=s12+s13+s23
      s14p24p34=s14+s24+s34

      aqmpC4x123m0_res =
     & -2*s14p24p34*zb(p1,p3)*zab2(p4,p2,p3,p1)*s123
     & /(zb(p1,p2)*zab2(p4,p1,p2,p3)**3)
     & -s14p24p34*za(p2,p3)**2
     & /(za(p1,p2)*za(p3,p4)*zab2(p4,p1,p2,p3))
     & -s14p24p34*zb(p1,p4)**2
     & /(zb(p1,p2)*zb(p3,p4)*zab2(p4,p1,p2,p3))
     & +2*zb(p1,p4)**2*zb(p2,p4)/(zb(p1,p2)*zb(p2,p3)*zb(p3,p4))

      return

