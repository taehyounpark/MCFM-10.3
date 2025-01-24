!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^- 4_g^- +H
c     Implementation of Eq.~(7.7) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      include 'Inc/Cred.f'
      include 'Inc/I5to4.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      complex(dp):: ppmmD2x34x1_res,e2x3x4x1,zab2
cc
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      e2x3x4x1=mtsq*za(p3,p4)/zb(p3,p4)
     & *((s(p3,p4)-4*mtsq)*zb(p2,p3)*za(p3,p4)*zb(p4,p1)
     & /(za(p2,p3)*zb(p3,p4)*za(p4,p1))-zb(p1,p2)**2)

      ppmmD2x34x1_res=Cred(3,I5to4(p2,p3,p4))*e2x3x4x1
     & +0.5_dp*za(p3,p4)*zab2(p1,p3,p4,p2)*zab2(p2,p3,p4,p1)
     & /(za(p1,p2)*za(p1,p4)*za(p2,p3)*zb(p3,p4)**2)
     & *(s(p3,p4)-4._dp*mtsq)

      return

