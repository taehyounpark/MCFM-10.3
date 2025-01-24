!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_qb^+ 2_q^- 3_g^- 4_g^+ +H
c     Implementation of Eq.~(9.4) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      real(dp):: s123
      complex(dp):: aqmpC3x12_res,zab2

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      s123=s(p1,p2)+s(p1,p3)+s(p2,p3)

      aqmpC3x12_res=2*zb(p1,p3)*zab2(p3,p1,p2,p3)*zab2(p4,p2,p3,p1)*s123
     & /(zb(p1,p2)*zab2(p4,p1,p2,p3)**3)

      return

