!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      use pppmC12x34m0_generic
      use pppmC12x34m2_generic
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^+ 4_g^- +H
      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      complex(dp):: pppmC12x34_res,msqCoeff
      msqCoeff=pppmC12x34m2(p1,p2,p3,p4,za,zb)
      pppmC12x34_res=pppmC12x34m0(p1,p2,p3,p4,za,zb)+mtsq*msqCoeff
      return



