!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module ppmmC2x3_generic
      implicit none
      public ppmmC2x3,ppmmC2x3_qp

      interface ppmmC2x3
      module procedure ppmmC2x3,ppmmC2x3_qp
      end interface

      contains

      function ppmmC2x3(p1,p2,p3,p4,za,zb) result(ppmmC2x3_res)
      use double_precision
      use sprod_dp
      include 'Inc/ppmmC2x3_inc.f'
      end function ppmmC2x3

      function ppmmC2x3_qp(p1,p2,p3,p4,za,zb) result(ppmmC2x3_res)
      use quad_precision
      use sprod_qp
      include 'Inc/ppmmC2x3_inc.f'
      end function ppmmC2x3_qp

      end module ppmmC2x3_generic

