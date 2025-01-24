!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      use aqppC3x412m0_generic
      use aqppC3x412m2_generic
      implicit none
c     This is a calculation for the process 1_qb^+ 2_q^- 3_g^+ 4_g^+ +H
c     Implementation of Eq.~(8.9) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      complex(dp):: aqppC3x412_res,msqCoeff
      msqCoeff=aqppC3x412m2(p1,p2,p3,p4,za,zb)
      aqppC3x412_res=aqppC3x412m0(p1,p2,p3,p4,za,zb)+mtsq*msqCoeff
      return

