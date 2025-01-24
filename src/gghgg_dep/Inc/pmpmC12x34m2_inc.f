!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      use pmpmC12x34m2part_generic
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^- 3_g^+ 4_g^- +H

      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      complex(dp):: pmpmC12x34m2_res
      pmpmC12x34m2_res=
     & +pmpmC12x34m2part(p1,p2,p3,p4,za,zb)
     & +pmpmC12x34m2part(p3,p4,p1,p2,za,zb)
     & +pmpmC12x34m2part(p2,p1,p4,p3,zb,za)
     & +pmpmC12x34m2part(p4,p3,p2,p1,zb,za)
      return



