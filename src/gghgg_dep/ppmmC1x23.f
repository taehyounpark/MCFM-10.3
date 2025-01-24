!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module ppmmC1x23_generic
      implicit none
      public ppmmC1x23,ppmmC1x23_qp

      interface ppmmC1x23
      module procedure ppmmC1x23,ppmmC1x23_qp
      end interface

      contains

      function ppmmC1x23(p1,p2,p3,p4,za,zb) result(ppmmC1x23_res)
      use double_precision
      use sprod_dp
      include 'Inc/ppmmC1x23_inc.f'
      end function ppmmC1x23

      function ppmmC1x23_qp(p1,p2,p3,p4,za,zb) result(ppmmC1x23_res)
      use quad_precision
      use sprod_qp
      include 'Inc/ppmmC1x23_inc.f'
      end function ppmmC1x23_qp

      end module ppmmC1x23_generic

