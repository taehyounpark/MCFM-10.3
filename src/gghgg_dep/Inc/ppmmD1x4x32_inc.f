!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^- 4_g^- +H
c     Implementation of Eq.~(7.6) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      include 'Inc/Cred.f'
      include 'Inc/I5to4.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      real(dp):: s14,s23,s24,s34,s234,C2
      complex(dp):: ppmmD1x4x32_res,zab2,
     & e2x3x4x1,zab1342,zab4231,zab1234,zab3241

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      s14=s(p1,p4)
      s23=s(p2,p3)
      s24=s(p2,p4)
      s34=s(p3,p4)
      s234=s23+s24+s34

      C2=Cred(2,I5to4(p2,p3,p4))
      e2x3x4x1=mtsq*za(p3,p4)/zb(p3,p4)
     & *((s34-4._dp*mtsq)*zb(p2,p3)*za(p3,p4)*zb(p4,p1)
     & /(za(p2,p3)*zb(p3,p4)*za(p4,p1))-zb(p1,p2)**2)

      zab1342=zab2(p1,p3,p4,p2)
      zab4231=zab2(p4,p2,p3,p1)
      zab1234=zab2(p1,p2,p3,p4)
      zab3241=zab2(p3,p2,p4,p1)

      ppmmD1x4x32_res=C2*e2x3x4x1
     & -zb(p2,p4)/(zb(p3,p4)*zab1234)*(
     & 2*s234**2*s14*zb(p2,p4)*zab1342/(zb(p2,p3)*zab1234**2)
     & +0.5_dp*s234
     & *(zb(p1,p4)*za(p3,p4)**2/za(p2,p3)
     &  +za(p1,p4)*zb(p1,p2)**2/zb(p2,p3))
     & +2._dp*mtsq
     & *(3._dp*zb(p2,p4)*zab1342*zab4231/(zb(p2,p3)*zab1234)
     & +(zb(p1,p4)*za(p3,p4)*s234-s34*zab3241)/(za(p2,p3)*zb(p3,p4))))

      return

