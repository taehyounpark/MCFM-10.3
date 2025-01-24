!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function qloop_d56x12x34m2(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
c     This is an implementation of deLaurentis, Campbell and Ellis
c     arXiv:2203.17170, Eq(4.16)
c     this function returns simplified qloop_d56x12x34m2
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp):: qloop_d56x12x34m2,qloop_d12x34x56m2_sym,
     & qloop_diff
      qloop_d56x12x34m2=
     & +qloop_diff(p1,p2,p3,p4,p5,p6,p7,za,zb)
     & +qloop_d12x34x56m2_sym(p5,p6,p1,p2,p3,p4,p7,za,zb)
      return
      end
