!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^+ 4_g^- +H
c     Implementation of Eq.~(5.15) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      include 'Inc/Cred.f'
      include 'Inc/I5to4.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      real(dp):: s123
      complex(dp):: pppmD1x23x4_res,pppme1x2x3x4,zab2

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      s123=s(p1,p2)+s(p1,p3)+s(p2,p3)

      pppme1x2x3x4=(s123-4._dp*mtsq)*mtsq
     & *zb(p2,p3)*zab2(p4,p2,p3,p1)/(za(p2,p3)*zab2(p1,p2,p3,p4))

      pppmD1x23x4_res=pppme1x2x3x4*Cred(3,I5to4(p1,p2,p3))
      return
