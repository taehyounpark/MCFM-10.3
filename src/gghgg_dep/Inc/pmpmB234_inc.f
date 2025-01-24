!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^- 3_g^+ 4_g^- +H
c     Implementation of Eq.~(6.18) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      real(dp):: s12,s13,s14,s23,s24,s34,s234,s12ps13ps14
      complex(dp):: pmpmB234_res,zab2,zab
      zab(p1,p2,p3)=za(p1,p2)*zb(p2,p3)
      zab2(p1,p2,p3,p4)=(za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4))

      s12=s(p1,p2)
      s13=s(p1,p3)
      s14=s(p1,p4)
      s23=s(p2,p3)
      s24=s(p2,p4)
      s34=s(p3,p4)
      s234=s23+s24+s34
      s12ps13ps14=s12+s13+s14

      pmpmB234_res=4._dp*s234*zb(p3,p4)/zb(p2,p4)
     & *(zab(p2,p4,p3)
     & /((s24+s34)*zab2(p1,p2,p3,p4)**2)

     & +zb(p1,p3)*zab2(p1,p2,p4,p3)
     & /(s12ps13ps14*zab2(p1,p2,p3,p4)*zb(p2,p3)*zb(p3,p4))
     & *(zb(p1,p3)/s12ps13ps14-zb(p3,p4)/zab2(p1,p2,p3,p4)))

     &  +4._dp*s234*zb(p3,p2)/zb(p4,p2)
     & *(zab(p4,p2,p3)
     & /((s24+s23)*zab2(p1,p3,p4,p2)**2)

     & +zb(p1,p3)*zab2(p1,p4,p2,p3)
     & /(s12ps13ps14*zab2(p1,p4,p3,p2)*zb(p4,p3)*zb(p3,p2))
     & *(zb(p1,p3)/s12ps13ps14-zb(p3,p2)/zab2(p1,p4,p3,p2)))

      return
