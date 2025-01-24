!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      use ppmmC1x234m2_generic
      use ppmmC1x234m0_generic
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^- 4_g^- +H
      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      complex(dp):: ppmmC1x234_res,msqCoeff
      msqCoeff=ppmmC1x234m2(p1,p2,p3,p4,za,zb)
      ppmmC1x234_res=ppmmC1x234m0(p1,p2,p3,p4,za,zb)+mtsq*msqCoeff
      return


