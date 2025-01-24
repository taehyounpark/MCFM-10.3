!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_qb^+ 2_q^- 3_g^- 4_g^+ +H
c     Implementation of Eq.~(9.2) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      real(dp):: s34,s123
      complex(dp):: aqmpD4x3x21_res,zab2

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      s34=s(p3,p4)
      s123=s(p1,p2)+s(p1,p3)+s(p2,p3)

      aqmpD4x3x21_res =
     & -2*zb(p1,p3)*zab2(p4,p2,p3,p1)*s34*s123**2
     & /(zb(p1,p2)*zab2(p4,p1,p2,p3)**3)
     & +0.5_dp*za(p2,p3)**2*zb(p3,p4)*s123
     & /(za(p1,p2)*zab2(p4,p1,p2,p3))
     & +0.5_dp*zb(p1,p4)**2*za(p3,p4)*s123
     & /(zb(p1,p2)*zab2(p4,p1,p2,p3))
     & +mtsq*(-6*zb(p1,p3)*zab2(p3,p1,p2,p4)*zab2(p4,p2,p3,p1)
     & /(zb(p1,p2)*zab2(p4,p1,p2,p3)**2)
     & -2*za(p2,p3)*zab2(p2,p1,p3,p4)/(za(p1,p2)*zab2(p4,p1,p2,p3)))

      return
