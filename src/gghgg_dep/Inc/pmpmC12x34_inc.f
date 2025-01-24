!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      use pmpmC12x34m2_generic
      use pmpmC12x34m0diff_generic
      use scpmpmC12x34m0_generic
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^- 3_g^+ 4_g^- +H
c     Implementation of Eq.~(6.11) from arXiv:2002.04018 v2

      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      complex(dp)::pmpmC12x34_res,msqCoeff

      msqCoeff=pmpmC12x34m2(p1,p2,p3,p4,za,zb)
      pmpmC12x34_res= +scpmpmC12x34m0(p1,p2,p3,p4,za,zb)
     &            +pmpmC12x34m0diff(p1,p2,p3,p4,za,zb)
     &            +pmpmC12x34m0diff(p3,p4,p1,p2,za,zb)
     &            +pmpmC12x34m0diff(p2,p1,p4,p3,zb,za)
     &            +pmpmC12x34m0diff(p4,p3,p2,p1,zb,za)
     &            +mtsq*msqCoeff
      return

