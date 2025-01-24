!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^+ 4_g^- +H
c     Implementation of Eq.~(5.16) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      include 'Inc/Cred.f'
      include 'Inc/I5to4.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      real(dp):: s123,s234,s1234,fac
      complex(dp):: pppmD2x3x4_res,pppme1x2x3x4,mpppe2x3x4x1,zab2

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      s123=s(p1,p2)+s(p1,p3)+s(p2,p3)
      s234=s(p2,p3)+s(p2,p4)+s(p3,p4)
      s1234=s123+s234-s(p2,p3)+s(p1,p4)
      fac=(s123-4._dp*mtsq)*mtsq

      pppme1x2x3x4=fac
     & *zb(p2,p3)*zab2(p4,p2,p3,p1)/(za(p2,p3)*zab2(p1,p2,p3,p4))
      mpppe2x3x4x1=-mtsq*zb(p2,p3)/(za(p2,p3)*zb(p3,p4))*(
     & zb(p2,p3)*zab2(p2,p3,p4,p1)*zab2(p4,p1,p3,p2)/zab2(p1,p3,p4,p2)
     & +zb(p1,p3)*zab2(p4,p2,p3,p1)
     & +4._dp*mtsq*zb(p2,p3)*za(p3,p4)*zab2(p2,p3,p4,p1)
     & /(za(p2,p3)*zab2(p1,p3,p4,p2)))

      pppmD2x3x4_res=(4._dp*mtsq*s234-s(p2,p3)*s1234)*za(p3,p4)*zb(p2,p3)**2
     & /(2._dp*za(p2,p3)*zab2(p1,p3,p4,p2)*zab2(p1,p2,p3,p4))
     &     +pppme1x2x3x4*Cred(1,I5to4(p1,p2,p3))
     &     +mpppe2x3x4x1*Cred(5,I5to4(p2,p3,p4))

      return

