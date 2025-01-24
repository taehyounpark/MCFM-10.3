!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^+ 4_g^- +H
c     Implementation of Eq.~(5.17) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      include 'Inc/Cred.f'
      include 'Inc/I5to4.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      real(dp):: s123,fac
      complex(dp):: pppmD1x2x3_res,pppme1x2x3x4,mpppe4x1x2x3,zab2

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      s123=s(p1,p2)+s(p1,p3)+s(p2,p3)
      fac=(s123-4._dp*mtsq)*mtsq

      pppme1x2x3x4=fac
     & *zb(p2,p3)*zab2(p4,p2,p3,p1)/(za(p2,p3)*zab2(p1,p2,p3,p4))
      mpppe4x1x2x3=fac
     & *zb(p2,p1)*zab2(p4,p2,p1,p3)/(za(p2,p1)*zab2(p3,p2,p1,p4))

      pppmD1x2x3_res=(4._dp*mtsq-s123)*(s123*zb(p1,p2)*zb(p2,p3))
     & /(2._dp*zab2(p3,p1,p2,p4)*zab2(p1,p2,p3,p4))
     &     +pppme1x2x3x4*Cred(5,I5to4(p1,p2,p3))
     &     +mpppe4x1x2x3*Cred(1,I5to4(p4,p1,p2))

      return

