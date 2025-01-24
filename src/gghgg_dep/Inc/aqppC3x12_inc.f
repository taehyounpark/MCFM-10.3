!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_qb^+ 2_q^- 3_g^+ 4_g^+ +H
c     Implementation of Eq.~(8.4) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      complex(dp):: aqppC3x12_res
      aqppC3x12_res=2*za(p2,p3)*za(p2,p4)*(s(p1,p3)+s(p2,p3))
     & /(za(p1,p2)*za(p3,p4)**3)
      return
