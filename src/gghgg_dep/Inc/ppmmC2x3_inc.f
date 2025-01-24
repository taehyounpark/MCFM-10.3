!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^- 4_g^- +H
c     Implementation of Eq.~(7.10) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'

      integer p1,p2,p3,p4
      real(dp):: s12,s13,s14,s23,s24,s34,s124,s134
      complex(dp)::ppmmC2x3_res,zab2

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      s12=s(p1,p2)
      s13=s(p1,p3)
      s14=s(p1,p4)
      s23=s(p2,p3)
      s24=s(p2,p4)
      s34=s(p3,p4)

      s134=s13+s14+s34
      s124=s12+s14+s24

      ppmmC2x3_res=
     &   za(p1,p2)*zb(p1,p3)**2*za(p1,p4)*s134*zab2(p2,p3,p4,p1)
     &  +zb(p3,p4)*za(p2,p4)**2*zb(p1,p4)*s124*zab2(p4,p1,p2,p3)
      ppmmC2x3_res=ppmmC2x3_res
     & *2._dp*s23/(s14*zab2(p2,p1,p4,p3)**3*za(p1,p2)*zb(p3,p4))
      return

