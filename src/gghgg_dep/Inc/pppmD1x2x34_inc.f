!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^+ 4_g^- +H
c     Implementation of Eq.~(5.10) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      include 'Inc/Cred.f'
      include 'Inc/I5to4.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      real(dp):: s12,s13,s14,s23,s24,s34,s123,s1234
      complex(dp):: pppmD1x2x34_res,zab2,e1x2x3x4

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      s12=s(p1,p2)
      s13=s(p1,p3)
      s14=s(p1,p4)
      s23=s(p2,p3)
      s24=s(p2,p4)
      s34=s(p3,p4)

      s123=s12+s13+s23
      s1234=s12+s13+s14+s23+s24+s34

      e1x2x3x4=(s123-4._dp*mtsq)*mtsq
     & *zab2(p4,p2,p3,p1)/zab2(p1,p2,p3,p4)*zb(p2,p3)/za(p2,p3)

      pppmD1x2x34_res=Cred(4,I5to4(p1,p2,p3))*e1x2x3x4
     & +0.5_dp*zb(p1,p2)*zb(p2,p3)*zab2(p1,p2,p4,p3)
     & *(zb(p2,p3)*s1234-4._dp*zab2(p1,p2,p4,p3)*mtsq/za(p1,p2))
     & /(zb(p3,p4)*zab2(p1,p3,p4,p2)*zab2(p1,p2,p3,p4))

     & -0.5_dp*zb(p1,p2)*za(p2,p4)*zab2(p4,p2,p3,p1)
     & *(zab2(p4,p2,p3,p1)+4._dp*mtsq*za(p2,p4)/za(p1,p2))
     & /(za(p2,p3)*za(p3,p4)*zab2(p2,p3,p4,p1))

      return
