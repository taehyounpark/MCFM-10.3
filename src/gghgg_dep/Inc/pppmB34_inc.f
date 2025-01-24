!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^+ 4_g^- +H
c     Implementation of Eq.~(5.30) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'

      complex(dp)::pppmB34_res
      integer p1,p2,p3,p4
      real(dp)::s13,s14,s23,s24

      s13=s(p1,p3)
      s14=s(p1,p4)
      s23=s(p2,p3)
      s24=s(p2,p4)

      pppmB34_res=4._dp
     & *(za(p2,p4)**2*za(p1,p3)*zb(p2,p3)/(s23+s24)
     &  -za(p1,p4)**2*za(p2,p3)*zb(p1,p3)/(s13+s14))
     & /(za(p1,p2)**2*za(p1,p3)*za(p2,p3))
      return
