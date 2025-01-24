!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_qb^+ 2_q^- 3_g^- 4_g^+ +H
c     Implementation of Eq.~(9.9) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      real(dp):: s12,s13,s14,s23,s24,s34,s123,s124,s1234,lam,lambda
      complex(dp):: zab2,aqpmmpB34_res

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      lambda(s12,s13,s14)=(s12-s13-s14)**2-4*s13*s14

      s12=s(p1,p2)
      s13=s(p1,p3)
      s14=s(p1,p4)
      s23=s(p2,p3)
      s24=s(p2,p4)
      s34=s(p3,p4)

      s123=s12+s13+s23
      s124=s12+s14+s24
      s1234=s12+s13+s14+s23+s24+s34
      lam=lambda(s1234,s12,s34)
      aqpmmpB34_res=
     & 2._dp*zab2(p3,p1,p2,p4)/(s12*zab2(p4,p1,p2,p3)*lam)
     & *(3._dp*zab2(p2,p3,p4,p1)
     & *(s124**2-s123**2)*(s13+s14+s23+s24)/lam
     & -3._dp*(s123+s124)*(za(p2,p3)*zb(p1,p3)-za(p2,p4)*zb(p1,p4))
     &  +2._dp*(s124-s123)/zab2(p4,p1,p2,p3)
     & *(za(p1,p2)*zb(p1,p3)*zab2(p4,p2,p3,p1)
     &  -za(p2,p4)*zb(p1,p2)*zab2(p2,p1,p4,p3)))

      return
