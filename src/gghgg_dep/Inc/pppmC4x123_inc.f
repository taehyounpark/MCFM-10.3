!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c---- Coefficient for the fermionic case
      use pppmC4x123m2_generic
      use pppmC4x123m0_generic
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^+ 4_g^- +H
      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      complex(dp):: pppmC4x123_res,msqCoeff
      msqCoeff=pppmC4x123m2(p1,p2,p3,p4,za,zb)
      pppmC4x123_res=pppmC4x123m0(p1,p2,p3,p4,za,zb)+mtsq*msqCoeff
      return



