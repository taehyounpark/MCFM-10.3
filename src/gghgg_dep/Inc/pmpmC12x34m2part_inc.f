!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^- 3_g^+ 4_g^- +H
c     Implementation of Eq.~(6.12) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'

      integer p1,p2,p3,p4
      real(dp):: s12,s13,s14,s23,s24,s34,s1234,Delta
      complex(dp)::pmpmC12x34m2part_res,zab2,zab2341,zab1342,zab3124

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      s12=s(p1,p2)
      s13=s(p1,p3)
      s14=s(p1,p4)
      s23=s(p2,p3)
      s24=s(p2,p4)
      s34=s(p3,p4)
      s1234=s12+s13+s14+s23+s24+s34
      Delta=(s1234-s12-s34)**2-4._dp*s12*s34
      zab2341=zab2(p2,p3,p4,p1)
      zab1342=zab2(p1,p3,p4,p2)
      zab3124=zab2(p3,p1,p2,p4)
      pmpmC12x34m2part_res=
     & 4._dp*zab2341/(zab1342*zab3124)
     & *(za(p2,p4)*(s13+s23-s14-s24)/Delta
     & *(zb(p2,p3)-za(p1,p4)*(s1234-s12-s34)
     & /(2._dp*za(p1,p2)*za(p3,p4)))
     & +zb(p2,p3)**2*(s23-s14)/(zb(p1,p2)*zb(p3,p4)*zab1342)
     & +1.5_dp*zb(p1,p3)*zb(p2,p3)/(zb(p1,p2)*zb(p3,p4)))
      return

