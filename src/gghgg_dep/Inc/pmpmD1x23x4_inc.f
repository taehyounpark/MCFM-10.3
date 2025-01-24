!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^- 3_g^+ 4_g^- +H
c     Implementation of Eq.~(6.7) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      include 'Inc/Cred.f'
      include 'Inc/I5to4.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      complex(dp):: pmpmD1x23x4_res,pmpme1x2x3x4,zab2

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      pmpme1x2x3x4=
     & -mtsq*zb(p3,p4)*zab2(p4,p2,p3,p1)/(za(p3,p4)*zab2(p1,p2,p3,p4))
     & *za(p2,p4)**2
     & *(1._dp+4._dp*mtsq*za(p1,p2)/(za(p2,p4)*zab2(p1,p2,p3,p4)))

     & -mtsq*za(p2,p1)*zab2(p4,p2,p3,p1)/(zb(p2,p1)*zab2(p1,p2,p3,p4))
     & *zb(p3,p1)**2
     & *(1._dp+4._dp*mtsq*zb(p4,p3)/(zb(p3,p1)*zab2(p1,p2,p3,p4)))

      pmpmD1x23x4_res=pmpme1x2x3x4*Cred(3,I5to4(p1,p2,p3))
     & +0.5_dp*za(p2,p4)**3*zab2(p4,p2,p3,p1)
     & /(za(p2,p3)*za(p3,p4)*za(p4,p1))
     & *(1._dp+4._dp*mtsq*za(p1,p2)/(za(p2,p4)*zab2(p1,p2,p3,p4)))
     & +0.5_dp*zb(p3,p1)**3*zab2(p4,p2,p3,p1)
     & /(zb(p3,p2)*zb(p2,p1)*zb(p1,p4))
     & *(1._dp+4._dp*mtsq*zb(p4,p3)/(zb(p3,p1)*zab2(p1,p2,p3,p4)))

      return

