!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      use ppppC1x234m2_generic
      use ppppC1x234m0_generic
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^+ 4_g^+ +H
      include 'Inc/zprods_decl.f'
      complex(dp)::ppppC1x234_res,msqCoeff
      integer p1,p2,p3,p4
      real(dp)::mtsq
      msqCoeff=ppppC1x234m2(p1,p2,p3,p4,za,zb)
      ppppC1x234_res=ppppC1x234m0(p1,p2,p3,p4,za,zb)+mtsq*msqCoeff
      return



