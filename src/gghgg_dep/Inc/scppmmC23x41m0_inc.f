!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      use scpmpmC12x34m0_generic
      implicit none
c     Scalar loop
c     This is a calculation for the process 1_g^+ 2_g^- 3_g^- 4_g^- +H
c     Implementation of Eq.~(7.13) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      complex(dp):: scppmmC23x41m0_res,zab2,den
      real(dp):: s13,s14,s23,s24,s1234,Delta

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      s13=s(p1,p3)
      s14=s(p1,p4)
      s23=s(p2,p3)
      s24=s(p2,p4)
      s1234=s(p1,p2)+s13+s14+s23+s24+s(p3,p4)
      Delta=(s1234-s14-s23)**2-4._dp*s14*s23
      den=zab2(p2,p1,p4,p3)*zab2(p1,p2,p3,p4)

      scppmmC23x41m0_res=-scpmpmC12x34m0(p2,p3,p1,p4,za,zb)
     &               -2._dp*Delta*((s13-s24)/den)**2
     &               -4._dp*zab2(p3,p1,p4,p2)*zab2(p4,p2,p3,p1)/den

      return

