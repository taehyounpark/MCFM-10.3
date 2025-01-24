!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^+ 4_g^- +H
c     Implementation of Eq.~(5.20) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'

      integer p1,p2,p3,p4
      real(dp):: s23,s24
      complex(dp):: pppmC2x34_res

      s23=s(p2,p3)
      s24=s(p2,p4)

      pppmC2x34_res=-2._dp*(s23+s24)*za(p1,p4)*za(p2,p4)**2
     & /(za(p1,p2)**3*za(p2,p3)*za(p3,p4))

      return

