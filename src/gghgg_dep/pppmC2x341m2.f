!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pppmC2x341m2_generic
      implicit none
      public pppmC2x341m2,pppmC2x341m2_qp

      interface pppmC2x341m2
      module procedure pppmC2x341m2,pppmC2x341m2_qp
      end interface

      contains

      function pppmC2x341m2(p1,p2,p3,p4,za,zb) result(pppmC2x341m2_res)
      use double_precision
      use sprod_dp
      include 'Inc/pppmC2x341m2_inc.f'
      end function pppmC2x341m2

      function pppmC2x341m2_qp(p1,p2,p3,p4,za,zb) result(pppmC2x341m2_res)
      use quad_precision
      use sprod_qp
      include 'Inc/pppmC2x341m2_inc.f'
      end function pppmC2x341m2_qp

      end module pppmC2x341m2_generic

