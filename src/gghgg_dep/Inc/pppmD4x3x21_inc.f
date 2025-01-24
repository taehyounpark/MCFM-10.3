!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^+ 4_g^- +H
c     Implementation of Eq.~(5.14) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      include 'Inc/Cred.f'
      include 'Inc/I5to4.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      real(dp):: s34,s123
      complex(dp):: pppmD4x3x21_res,zab2,e1x2x3x4

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      s34=s(p3,p4)
      s123=s(p1,p2)+s(p1,p3)+s(p2,p3)
      e1x2x3x4=(s123-4._dp*mtsq)*mtsq
     & *zab2(p4,p2,p3,p1)/zab2(p1,p2,p3,p4)*zb(p2,p3)/za(p2,p3)

      pppmD4x3x21_res=(4._dp*mtsq-s123)*(
     & +0.5_dp*s34*s123**2
     & /(za(p1,p2)*za(p2,p3)*zab2(p1,p2,p3,p4)*zab2(p3,p1,p2,p4)))
     & +e1x2x3x4*Cred(2,I5to4(p1,p2,p3))

      return

