!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      use pppmC2x341m0_generic
      use pppmC2x341m2_generic
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^- 4_g^- +H
      include 'Inc/zprods_decl.f'
      complex(dp)::pppmC2x341_res,msqCoeff
      integer p1,p2,p3,p4
      real(dp)::mtsq
      msqCoeff=pppmC2x341m2(p1,p2,p3,p4,za,zb)
      pppmC2x341_res=pppmC2x341m0(p1,p2,p3,p4,za,zb)+mtsq*msqCoeff
      return





