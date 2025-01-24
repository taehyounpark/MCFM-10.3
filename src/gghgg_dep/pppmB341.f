!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pppmB341_generic
      implicit none
      public pppmB341,pppmB341_qp

      interface pppmB341
      module procedure pppmB341,pppmB341_qp
      end interface

      contains

      function pppmB341(p1,p2,p3,p4,za,zb) result(pppmB341_res)
      use double_precision
      use sprod_dp
      use pppmB234_generic
      include 'Inc/pppmB341_inc.f'
      end function pppmB341

      function pppmB341_qp(p1,p2,p3,p4,za,zb) result(pppmB341_res)
      use quad_precision
      use sprod_qp
      use pppmB234_generic
      include 'Inc/pppmB341_inc.f'
      end function pppmB341_qp

      end module pppmB341_generic

