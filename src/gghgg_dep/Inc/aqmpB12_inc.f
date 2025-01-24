!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      use aqmpB12_unsym_generic
      implicit none
      complex(dp):: aqmpB12_res
c     This is a calculation for the process 1_qb^+ 2_q^- 3_g^- 4_g^+ +H
c     Implementation of Eq.~(9.10) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4

      aqmpB12_res = (+aqmpB12_unsym(p1,p2,p3,p4,za,zb)
     &               -aqmpB12_unsym(p2,p1,p4,p3,zb,za))
      return


