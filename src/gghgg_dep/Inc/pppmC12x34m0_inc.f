!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^+ 4_g^- +H
c     Implementation of Eq.~(5.28) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'

      integer p1,p2,p3,p4
      real(dp):: s123,s124,s1234,s12,s13,s14,s23,s24,s34
      complex(dp):: pppmC12x34m0_res,zab2,
     & zab1342,zab4231,zab4132,zab2341,zab3124

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      s12=s(p1,p2)
      s13=s(p1,p3)
      s14=s(p1,p4)
      s23=s(p2,p3)
      s24=s(p2,p4)
      s34=s(p3,p4)

      s123=s12+s13+s23
      s124=s12+s14+s24
      s1234=s12+s13+s14+s23+s24+s34

      zab1342 = zab2(p1,p3,p4,p2)
      zab4231 = zab2(p4,p2,p3,p1)
      zab4132 = zab2(p4,p1,p3,p2)
      zab2341 = zab2(p2,p3,p4,p1)
      zab3124 = zab2(p3,p1,p2,p4)

      pppmC12x34m0_res=
     & +(zb(p1,p3)*(s12+s23)*zab4231)
     & /(za(p1,p2)*za(p2,p3)*zb(p2,p4)*zab2341)
     & +(zb(p1,p2)*za(p2,p4)**2*zab4231)
     & /(za(p1,p2)*za(p2,p3)*za(p3,p4)*zab2341)
     & -(zb(p2,p3)*zab4132**2)/(za(p1,p2)*zb(p2,p4)*za(p3,p4)*zab1342)
     & -(zb(p2,p3)**3*s1234)/(za(p1,p2)*zb(p2,p4)*zb(p3,p4)*zab1342)
     & -(zb(p1,p3)**2*zb(p2,p3)*s1234)
     & /(za(p1,p2)*zb(p2,p4)*zb(p3,p4)*zab2341)
     & +(zb(p1,p2)*s123*(s123-s124))
     & /(za(p1,p2)*za(p2,p3)*zb(p2,p4)*zab3124)
     & -(zb(p1,p2)*za(p1,p4)*zab4231)
     & /(za(p1,p2)*za(p2,p3)*zb(p2,p4)*za(p3,p4))

      return


