!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^- 3_g^+ 4_g^- +H
c     Implementation of Eq.~(6.10) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'

      integer p1,p2,p3,p4
      real(dp):: s234
      complex(dp):: pmpmC2x34_res,zab2

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      s234=s(p2,p3)+s(p2,p4)+s(p3,p4)
      pmpmC2x34_res=2._dp*s234*zab2(p2,p3,p4,p2)*zb(p2,p3)**2
     &   *(zab2(p1,p2,p4,p3)*zb(p2,p4)+zab2(p1,p3,p4,p2)*zb(p3,p4))
     &   /(zab2(p1,p3,p4,p2)**3*zb(p2,p4)**2*zb(p3,p4))
      return

