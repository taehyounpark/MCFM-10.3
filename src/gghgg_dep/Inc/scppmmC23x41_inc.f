!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      use scppmmC23x41m0_generic
      use ppmmC23x41m2_generic
      implicit none
c     Scalar loop
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^- 4_g^- +H
      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      complex(dp):: scppmmC23x41_res,msqCoeff
cc
      msqCoeff=ppmmC23x41m2(p1,p2,p3,p4,za,zb)
      scppmmC23x41_res=scppmmC23x41m0(p1,p2,p3,p4,za,zb)
     & +mtsq*msqCoeff

      return




