!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      use scpmpmC12x34m0_generic
      use pmpmC12x34m2_generic
      implicit none
c     Scalar loop
c     This is a calculation for the process 1_g^+ 2_g^- 3_g^+ 4_g^- +H
c     Implementation of Eq.~(6.13) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      complex(dp):: scpmpmC12x34_res,msqCoeff
cc
      msqCoeff=pmpmC12x34m2(p1,p2,p3,p4,za,zb)
      scpmpmC12x34_res=scpmpmC12x34m0(p1,p2,p3,p4,za,zb)
     & +mtsq*msqCoeff
      return
