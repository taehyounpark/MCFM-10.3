!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^- 3_g^+ 4_g^- +H
c     Implementation of Eq.~(6.9) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'

      integer p1,p2,p3,p4
      real(dp):: s12,s13,s14,s23,s24,s34,s234,s134,
     & s13ps23,s23ps34,s14ps24,s14ps34
      complex(dp):: pmpmC3x4_res,zab2,zab

      zab(p1,p2,p3)=za(p1,p2)*zb(p2,p3)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      s12=s(p1,p2)
      s13=s(p1,p3)
      s14=s(p1,p4)
      s23=s(p2,p3)
      s24=s(p2,p4)
      s34=s(p3,p4)

      s234=s23+s24+s34
      s134=s13+s14+s34
      s13ps23=s13+s23
      s23ps34=s23+s34
      s14ps24=s14+s24
      s14ps34=s14+s34

      pmpmC3x4_res=-2._dp*s34
     & /(s12*zab2(p3,p1,p2,p4)**3*(za(p1,p3)*zb(p2,p4))**2)
     & *(zab(p3,p1,p4)**3*(s24*s14ps24+s12*s23ps34)
     &  +zab(p3,p1,p4)**2*zab(p3,p2,p4)
     & *(s12**2-s14*s14ps24+s12*(3._dp*s234-5._dp*s24))
     & +zab(p3,p2,p4)**3*(s13*s13ps23+s12*s14ps34)
     & +zab(p3,p1,p4)*zab(p3,p2,p4)**2
     & *(s12**2-s23*s13ps23+s12*(3._dp*s134-5._dp*s13)))
      return

