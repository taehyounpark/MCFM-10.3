!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^- 4_g^- +H
c     Implementation of Eq.~(7.8) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      include 'Inc/Cred.f'
      include 'Inc/I5to4.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      real(dp):: s12,s34,C3
      complex(dp):: ppmmD1x23x4_res,zab2,e1x2x3x4,zab1234,zab4231

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      s12=s(p1,p2)
      s34=s(p3,p4)

      C3=Cred(3,I5to4(p1,p2,p3))

      e1x2x3x4=mtsq*(s12+s34-4._dp*mtsq)*zb(p1,p2)*za(p3,p4)
     & /(za(p1,p2)*zb(p3,p4))

      zab1234=zab2(p1,p2,p3,p4)
      zab4231=zab2(p4,p2,p3,p1)

      ppmmD1x23x4_res = C3*e1x2x3x4
     & -0.5_dp*zab4231
     & *(zb(p2,p1)*zb(p2,p4)/(zb(p1,p4)*zb(p2,p3)*zb(p3,p4))
     & *(zb(p2,p1)-4*mtsq*zb(p2,p4)/zab1234)
     &  +za(p4,p3)*za(p1,p3)/(za(p2,p3)*za(p1,p4)*za(p1,p2))
     & *(za(p4,p3)-4*mtsq*za(p1,p3)/zab1234))

      return
