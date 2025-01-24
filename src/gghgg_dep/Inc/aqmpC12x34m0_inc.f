!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_qb^+ 2_q^- 3_g^- 4_g^+ +H
c     Implementation of Eq.~(9.5) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'

      integer p1,p2,p3,p4
      real(dp)::s12,s13,s14,s23,s24,s34,lambda,s123,s124,s1234,lam
      complex(dp):: aqmpC12x34m0_res,zab2
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

      aqmpC12x34m0_res=
     & 8*(s124-s123)*(s12+s34+2*s13+2*s23)
     & *za(p2,p4)*zb(p1,p3)*zab2(p3,p1,p2,p4)/(zab2(p4,p1,p2,p3)**2*lam)

     & +((9*s13-7*s23-s14-s24+4*s34)*za(p2,p4)*zb(p1,p4)
     &  -(9*s14-7*s24-s13-s23+4*s34)*za(p2,p3)*zb(p1,p3))
     & /(zab2(p4,p1,p2,p3)**2)

     & +12*s1234*((s13+s23)**2-(s14+s24)**2)
     & *zab2(p2,p3,p4,p1)*zab2(p3,p1,p2,p4)
     & /(zab2(p4,p1,p2,p3)*lam**2)

     & +4*zab2(p3,p1,p2,p4)/(zab2(p4,p1,p2,p3)*lam)
     & *((3*s12+3*s34+4*s13+4*s23+4*s14)*za(p2,p3)*zb(p1,p3)
     &  -(3*s12+3*s34+4*s13+4*s14+4*s24)*za(p2,p4)*zb(p1,p4))

     & -24*zb(p1,p3)*za(p2,p4)*zab2(p3,p1,p2,p4)**2
     & /(zab2(p4,p1,p2,p3)*lam)
     & -8*zb(p1,p4)*za(p2,p3)*zab2(p3,p1,p2,p4)/lam
     & +8*zb(p1,p4)*za(p2,p3)/zab2(p4,p1,p2,p3)

      aqmpC12x34m0_res=aqmpC12x34m0_res
     & +aqmpC12x34m0unsym(p1,p2,p3,p4,za,zb)
     & -aqmpC12x34m0unsym(p2,p1,p4,p3,zb,za)
      return

