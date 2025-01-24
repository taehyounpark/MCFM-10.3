!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_qb^+ 2_q^- 3_g^+ 4_g^+ +H
c     Implementation of Eq.~(8.10) from arXiv:2002.04018 v2

      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      real(dp):: s12,s13,s14,s23,s24,s34,s124,s13p23p34,Tr3x12x4x12
      complex(dp):: aqppC3x412m2_res,zab2

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      s12=s(p1,p2)
      s13=s(p1,p3)
      s14=s(p1,p4)
      s23=s(p2,p3)
      s24=s(p2,p4)
      s34=s(p3,p4)
      s124=s12+s14+s24
      s13p23p34=s13+s23+s34
      Tr3x12x4x12=(s13+s23)*(s14+s24)-s12*s34

      aqppC3x412m2_res=
     & -8*zab2(p2,p1,p4,p3)**2*zab2(p1,p2,p4,p3)
     & /(za(p2,p1)*za(p1,p4)*zab2(p4,p1,p2,p3)*s13p23p34*s124)
     & +4*za(p2,p4)*za(p2,p3)*zb(p4,p3)*s13p23p34
     & /(za(p2,p1)*za(p4,p3)**2*Tr3x12x4x12)
     & +4*zb(p1,p4)*zb(p1,p3)*s13p23p34
     & /(zb(p2,p1)*za(p4,p3)*Tr3x12x4x12)
      return
