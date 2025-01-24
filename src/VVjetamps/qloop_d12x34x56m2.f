!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function qloop_d12x34x56m2(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
c     this function returns simplified qloop_d12x34x56m2
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp):: qloop_d12x34x56m2,
     & qloop_d12x34x56m2_sym,qloop_d12x34x56m2_asy
      qloop_d12x34x56m2=
     & +qloop_d12x34x56m2_sym(p1,p2,p3,p4,p5,p6,p7,za,zb)
     & +qloop_d12x34x56m2_asy(p1,p2,p3,p4,p5,p6,p7,za,zb)
      return
      end
