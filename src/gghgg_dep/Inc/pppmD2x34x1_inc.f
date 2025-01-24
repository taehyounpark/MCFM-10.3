!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^+ 4_g^- +H
c     Implementation of Eq.~(5.13) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      include 'Inc/Cred.f'
      include 'Inc/I5to4.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      real(dp):: s123,s234,s1234
      complex(dp):: pppmD2x34x1_res,pppme2x3x4x1,zab2

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      s123=s(p1,p2)+s(p1,p3)+s(p2,p3)
      s234=s(p2,p3)+s(p2,p4)+s(p3,p4)
      s1234=s123+s234-s(p2,p3)+s(p1,p4)

      pppme2x3x4x1=-mtsq*zb(p2,p3)/(za(p2,p3)*zb(p3,p4))*(
     & zb(p2,p3)*zab2(p2,p3,p4,p1)*zab2(p4,p1,p3,p2)/zab2(p1,p3,p4,p2)
     &+zb(p1,p3)*zab2(p4,p2,p3,p1)
     &+4._dp*mtsq*zb(p2,p3)*za(p3,p4)*zab2(p2,p3,p4,p1)
     & /(za(p2,p3)*zab2(p1,p3,p4,p2)))

      pppmD2x34x1_res=pppme2x3x4x1*Cred(3,I5to4(p2,p3,p4))
     & +za(p2,p4)/(za(p1,p2)*za(p2,p3))*(
     &  2._dp*za(p1,p4)*za(p2,p4)*zab2(p1,p3,p4,p2)*zab2(p2,p3,p4,p1)
     &   /(za(p1,p2)**2*za(p3,p4))
     & +zab2(p4,p1,p3,p2)*zab2(p4,p2,p3,p1)/(2._dp*za(p3,p4))
     & +s1234*zb(p1,p3)*zb(p2,p3)/(2._dp*zb(p3,p4))
     & -2._dp*mtsq*(3._dp*za(p1,p4)*za(p2,p4)*zb(p1,p2)
     & /(za(p1,p2)*za(p3,p4))
     & +2._dp*zb(p1,p3)*zb(p2,p3)/zb(p3,p4)
     & +za(p2,p4)*zb(p1,p4)*zb(p2,p3)/(za(p2,p3)*zb(p3,p4)) ) )
      return
