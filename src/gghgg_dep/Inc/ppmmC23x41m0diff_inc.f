!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^- 4_g^- +H
c     Implementation of Eq.~(7.12) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      complex(dp):: ppmmC23x41m0diff_res,zab2

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      ppmmC23x41m0diff_res=
     & (za(p3,p4)**2*zab2(p2,p1,p4,p2)*zab2(p4,p2,p3,p4)
     & -(za(p2,p3)*zb(p2,p1)*za(p1,p4))**2)
     & /(za(p2,p3)*za(p1,p4)*zab2(p2,p1,p4,p3)*zab2(p1,p2,p3,p4))
      return

