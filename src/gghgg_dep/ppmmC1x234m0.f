!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module ppmmC1x234m0_generic
      implicit none
      public ppmmC1x234m0,ppmmC1x234m0_qp

      interface ppmmC1x234m0
      module procedure ppmmC1x234m0,ppmmC1x234m0_qp
      end interface

      contains

      function ppmmC1x234m0(p1,p2,p3,p4,za,zb) result(ppmmC1x234m0_res)
      use double_precision
      use sprod_dp
      include 'Inc/ppmmC1x234m0_inc.f'
      end function ppmmC1x234m0

      function ppmmC1x234m0_qp(p1,p2,p3,p4,za,zb) result(ppmmC1x234m0_res)
      use quad_precision
      use sprod_qp
      include 'Inc/ppmmC1x234m0_inc.f'
      end function ppmmC1x234m0_qp

      end module ppmmC1x234m0_generic

