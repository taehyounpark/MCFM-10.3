!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      use ppmmC23x41m2_generic
      use scppmmC23x41m0_generic
      use ppmmC23x41m0diff_generic
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^+ 3_g^- 4_g^- +H
c     Implementation of Eq.~(7.12) from arXiv:2002.04018 v2

      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      complex(dp):: ppmmC23x41_res,msqCoeff
c,scppmmC23x41m0
cscppmmC23x41m0,ppmmC23x41m2,ppmmC23x41m0diff,

      msqCoeff=ppmmC23x41m2(p1,p2,p3,p4,za,zb)
      ppmmC23x41_res=+scppmmC23x41m0(p1,p2,p3,p4,za,zb)
     &           +ppmmC23x41m0diff(p1,p2,p3,p4,za,zb)
     &           +ppmmC23x41m0diff(p2,p1,p4,p3,za,zb)
     &           +ppmmC23x41m0diff(p3,p4,p1,p2,zb,za)
     &           +ppmmC23x41m0diff(p4,p3,p2,p1,zb,za)
     &           +mtsq*msqCoeff

      return


