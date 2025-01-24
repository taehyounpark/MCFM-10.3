!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^- 4_g^- +H
c     Implementation of Eq.~(5.32) from arXiv:2002.04018 v2

      include 'Inc/zprods_decl.f'
      complex(dp)::pppmB1234_res
      integer p1,p2,p3,p4

      pppmB1234_res=
     & -pppmB34(p1,p2,p3,p4,za,zb)
     & -pppmB34(p3,p2,p1,p4,za,zb)
     & -pppmB234(p1,p2,p3,p4,za,zb)
     & -pppmB341(p1,p2,p3,p4,za,zb)
     & -pppmB234(p3,p2,p1,p4,za,zb)
      return
