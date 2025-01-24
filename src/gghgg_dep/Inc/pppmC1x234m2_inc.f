!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^+ 4_g^- +H
c     Implementation of Eq.~(5.25) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      complex(dp)::pppmC1x234m2_res
      integer p1,p2,p3,p4
      real(dp)::s12,s13,s14,s23,s24,s34,s234,s12p13p14
      complex(dp):: zab2,zab1234,zab4231,zab2341,zab1342
c---- statement function
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
c---- end statement function
      s12=s(p1,p2)
      s13=s(p1,p3)
      s14=s(p1,p4)
      s23=s(p2,p3)
      s24=s(p2,p4)
      s34=s(p3,p4)
      s234=s23+s24+s34
      s12p13p14=s12+s13+s14
      zab1234=zab2(p1,p2,p3,p4)
      zab4231=zab2(p4,p2,p3,p1)
      zab2341=zab2(p2,p3,p4,p1)
      zab1342=zab2(p1,p3,p4,p2)

      pppmC1x234m2_res=
     & +8._dp*zab4231**3/(s234*s12p13p14
     & *za(p2,p3)*za(p3,p4)*zab2341)

     & -4._dp*s12p13p14/(za(p1,p2)*za(p3,p4))
     & *(za(p1,p4)**2*zb(p2,p3)/(za(p1,p2)*zab1342*zab1234)

     & -zb(p2,p3)/(zab1234*zab1342*zab2341)
     & *(za(p1,p2)*zb(p1,p3)*zb(p2,p3)*za(p3,p4)/zb(p3,p4)
     &  +zb(p1,p2)*za(p1,p4)*za(p2,p4))

     & -(zb(p1,p3)*za(p1,p4)*za(p2,p4)/za(p1,p2)
     &  +za(p2,p4)*zab4231/za(p2,p3))
     & /(zab1234*zab2341))

      return

