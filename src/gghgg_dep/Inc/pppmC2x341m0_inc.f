!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^- 4_g^- +H
c     Implementation of Eq.~(5.26) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'

      complex(dp)::pppmC2x341m0_res,scpppmC2x341m0,zab2
      integer p1,p2,p3,p4
      real(dp):: s12,s13,s14,s23,s24,s34,s134

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      s12=s(p1,p2)
      s13=s(p1,p3)
      s14=s(p1,p4)
      s23=s(p2,p3)
      s24=s(p2,p4)
      s34=s(p3,p4)

      s134=s13+s14+s34

      scpppmC2x341m0=(2*za(p2,p4)**2
     & *(za(p1,p4)**2*za(p2,p3)**2+za(p1,p2)**2*za(p3,p4)**2)
     & *(s12+s23+s24))
     & /((za(p1,p2)*za(p2,p3))**3*za(p1,p4)*za(p3,p4))

      pppmC2x341m0_res=
     & +zab2(p4,p1,p3,p2)**4
     & /(za(p1,p4)*za(p3,p4)*zab2(p1,p3,p4,p2)*zab2(p3,p1,p4,p2)*s134)
     & +zb(p1,p3)**2*(s12+s23)
     & /(za(p1,p2)*zb(p1,p4)*za(p2,p3)*zb(p3,p4))

     & +zb(p1,p2)*zb(p1,p3)**3*(s12+s23+s24)
     & /(zb(p1,p4)*zb(p3,p4)*zab2(p2,p3,p4,p1)*s134)
     & -zb(p1,p2)*za(p2,p4)**2*zab2(p4,p3,p2,p1)
     & /(za(p1,p2)*za(p2,p3)*za(p3,p4)*zab2(p2,p3,p4,p1))
     & -zb(p1,p2)*zb(p1,p3)**2*za(p2,p4)
     & /(za(p2,p3)*zb(p3,p4)*zab2(p2,p3,p4,p1))

     & +zb(p3,p2)*zb(p3,p1)**3*(s12+s23+s24)
     & /(zb(p3,p4)*zb(p1,p4)*zab2(p2,p1,p4,p3)*s134)
     & -zb(p3,p2)*za(p2,p4)**2*zab2(p4,p1,p2,p3)
     & /(za(p3,p2)*za(p2,p1)*za(p1,p4)*zab2(p2,p1,p4,p3))
     & -zb(p3,p2)*zb(p3,p1)**2*za(p2,p4)
     & /(za(p2,p1)*zb(p1,p4)*zab2(p2,p1,p4,p3))

      pppmC2x341m0_res=pppmC2x341m0_res+scpppmC2x341m0
      return


