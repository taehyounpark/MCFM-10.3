!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      use ppmmC23x41m2_unsym_generic
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^- 4_g^- +H
c     Implementation of Eq.~(7.14) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      complex(dp):: ppmmC23x41m2_res

      ppmmC23x41m2_res=
     &        ppmmC23x41m2_unsym(p1,p2,p3,p4,za,zb)
     &       +ppmmC23x41m2_unsym(p3,p4,p1,p2,zb,za)
     &       +ppmmC23x41m2_unsym(p4,p3,p2,p1,zb,za)
     &       +ppmmC23x41m2_unsym(p2,p1,p4,p3,za,zb)

      return

