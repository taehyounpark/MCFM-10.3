!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^- 4_g^- +H
c     Implementation of Eq.~(7.9) from arXiv:2002.04018 v2
      implicit none
      complex(dp):: ppmmD1x2x3_res
      include 'Inc/zprods_decl.f'
      include 'Inc/Cred.f'
      include 'Inc/I5to4.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      real(dp):: s12,s34,C1,C5
      complex(dp)::e1x2x3x4,e4x1x2x3

      s12=s(p1,p2)
      s34=s(p3,p4)

      e1x2x3x4=mtsq*(s12+s34-4._dp*mtsq)*zb(p1,p2)*za(p3,p4)
     & /(za(p1,p2)*zb(p3,p4))
      C5=Cred(5,I5to4(p1,p2,p3))

      e4x1x2x3=mtsq*zb(p1,p2)/za(p1,p2)
     & *((s12-4*mtsq)*za(p4,p1)*zb(p1,p2)*za(p2,p3)
     &               /(zb(p4,p1)*za(p1,p2)*zb(p2,p3))-za(p3,p4)**2)
      C1=Cred(1,I5to4(p4,p1,p2))

      ppmmD1x2x3_res=C1*e4x1x2x3+C5*e1x2x3x4
     & +zb(p1,p2)**2*za(p2,p3)*(s12-4*mtsq)
     & /(2._dp*zb(p3,p4)*za(p1,p2)*zb(p1,p4))
