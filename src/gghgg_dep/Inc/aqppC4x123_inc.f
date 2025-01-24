!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      use aqppC4x123m0_generic
      use aqppC4x123m2_generic
      implicit none
c     This is a calculation for the process 1_qb^+ 2_q^- 3_g^+ 4_g^+ +H
c     Implementation of Eq.~(8.7) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      complex(dp):: aqppC4x123_res,msqCoeff
      msqCoeff=aqppC4x123m2(p1,p2,p3,p4,za,zb)
      aqppC4x123_res=aqppC4x123m0(p1,p2,p3,p4,za,zb)
     &     +mtsq*msqCoeff
      return



