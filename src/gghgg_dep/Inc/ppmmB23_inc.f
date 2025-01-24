!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^- 4_g^- +H
c     Implementation of Eq.~(7.17) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      complex(dp):: ppmmB23_res

      ppmmB23_res=ppmmB23symm(p1,p2,p3,p4,za,zb)
     &           +ppmmB23symm(p4,p3,p2,p1,zb,za)

      return


