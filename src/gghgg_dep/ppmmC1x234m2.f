!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module ppmmC1x234m2_generic
      implicit none
      public ppmmC1x234m2,ppmmC1x234m2_qp

      interface ppmmC1x234m2
      module procedure ppmmC1x234m2,ppmmC1x234m2_qp
      end interface

      contains

      function ppmmC1x234m2(p1,p2,p3,p4,za,zb) result(ppmmC1x234m2_res)
      use double_precision
      use sprod_dp
      include 'Inc/ppmmC1x234m2_inc.f'
      end function ppmmC1x234m2

      function ppmmC1x234m2_qp(p1,p2,p3,p4,za,zb) result(ppmmC1x234m2_res)
      use quad_precision
      use sprod_qp
      include 'Inc/ppmmC1x234m2_inc.f'
      end function ppmmC1x234m2_qp

      end module ppmmC1x234m2_generic

