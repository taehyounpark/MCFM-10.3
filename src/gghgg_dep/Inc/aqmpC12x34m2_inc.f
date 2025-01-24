!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_qb^+ 2_q^- 3_g^- 4_g^+ +H
c     Implementation of Eq.~(9.6) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      complex(dp):: aqmpC12x34m2_res
      aqmpC12x34m2_res=aqmpC12x34m2unsym(p1,p2,p3,p4,za,zb)
     &                -aqmpC12x34m2unsym(p2,p1,p4,p3,zb,za)
      return
