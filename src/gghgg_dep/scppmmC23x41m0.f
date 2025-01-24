!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module scppmmC23x41m0_generic
      implicit none
      public scppmmC23x41m0,scppmmC23x41m0_qp

      interface scppmmC23x41m0
      module procedure scppmmC23x41m0,scppmmC23x41m0_qp
      end interface

      contains

      function scppmmC23x41m0(p1,p2,p3,p4,za,zb) result(scppmmC23x41m0_res)
      use double_precision
      use sprod_dp
      include 'Inc/scppmmC23x41m0_inc.f'
c      scppmmC23x41m0_res=0
      end function scppmmC23x41m0

      function scppmmC23x41m0_qp(p1,p2,p3,p4,za,zb) result(scppmmC23x41m0_res)
      use quad_precision
      use sprod_qp
      include 'Inc/scppmmC23x41m0_inc.f'
c      scppmmC23x41m0_res=0
      end function scppmmC23x41m0_qp

      end module scppmmC23x41m0_generic

