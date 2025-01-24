!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^- 4_g^- +H
c     Implementation of Eq.~(7.5) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      include 'Inc/Cred.f'
      include 'Inc/I5to4.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      complex(dp):: ppmmD1x2x34_res,e1x2x3x4

      e1x2x3x4=mtsq*(s(p1,p2)+s(p3,p4)-4._dp*mtsq)*zb(p1,p2)*za(p3,p4)
     & /(za(p1,p2)*zb(p3,p4))
      ppmmD1x2x34_res=e1x2x3x4*Cred(4,I5to4(p1,p2,p3))
      return
