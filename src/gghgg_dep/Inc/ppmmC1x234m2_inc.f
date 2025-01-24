!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^- 4_g^- +H
c     Implementation of Eq.~(7.16) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'

      integer p1,p2,p3,p4
      real(dp):: s12,s13,s14,s23,s12ps13ps14
      complex(dp):: ppmmC1x234m2_res,zab,zab2,rat

      zab(p1,p2,p3)=za(p1,p2)*zb(p2,p3)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      s12=s(p1,p2)
      s13=s(p1,p3)
      s14=s(p1,p4)
      s23=s(p2,p3)

      s12ps13ps14=s12+s13+s14

      rat=zab2(p1,p3,p4,p2)/zab2(p1,p2,p3,p4)
      ppmmC1x234m2_res=1._dp/(s23*zab(p1,p4,p3))
     & *((zab(p3,p4,p1)-zab(p3,p2,p1))*zab(p4,p1,p2)**2
     & /(s14*s12ps13ps14)
     & +zb(p1,p2)*za(p3,p4)/zab2(p1,p2,p3,p4)
     & *(zab2(p1,p3,p4,p2)+s12ps13ps14/s12*zab(p1,p3,p2))
     & -2._dp*rat*zab2(p4,p2,p3,p1)*zab(p3,p1,p2)/s12ps13ps14
     & +rat**2*zab2(p4,p2,p3,p1)*zab(p3,p1,p4)/s12ps13ps14)
      ppmmC1x234m2_res=-4._dp*ppmmC1x234m2_res
      return
