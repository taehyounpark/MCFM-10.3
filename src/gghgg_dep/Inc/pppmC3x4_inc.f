!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^+ 4_g^- +H
c     Implementation of Eq.~(5.19) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      real(dp):: s34
      complex(dp):: pppmC3x4_res
      s34=s(p3,p4)
      pppmC3x4_res=2._dp*za(p1,p4)*za(p4,p3)*s34
     & /(za(p1,p2)*za(p2,p3)*za(p1,p3)**2)
      return


