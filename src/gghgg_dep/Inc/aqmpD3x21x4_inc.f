!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_qb^+ 2_q^- 3_g^- 4_g^+ +H
c     Implementation of Eq.~(9.1) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      complex(dp):: aqmpD3x21x4_res,zab2

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      aqmpD3x21x4_res = zab2(p3,p1,p2,p4)/2._dp*(
     &  +zb(p1,p4)**2/(zb(p1,p2)*zb(p3,p4))
     &  -za(p2,p3)**2/(za(p1,p2)*za(p3,p4))
     &  +4._dp*mtsq/zab2(p4,p1,p2,p3)
     & *(za(p2,p3)*za(p2,p4)/(za(p1,p2)*za(p3,p4))
     &  -zb(p1,p3)*zb(p1,p4)/(zb(p1,p2)*zb(p3,p4))))

      return
