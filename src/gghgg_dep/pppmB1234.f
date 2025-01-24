!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pppmB1234_generic
      implicit none
      public pppmB1234,pppmB1234_qp

      interface pppmB1234
      module procedure pppmB1234,pppmB1234_qp
      end interface

      contains

      function pppmB1234(p1,p2,p3,p4,za,zb) result(pppmB1234_res)
      use double_precision
      use sprod_dp
      use pppmB34_generic
      use pppmB341_generic
      use pppmB234_generic
      include 'Inc/pppmB1234_inc.f'
      end function pppmB1234

      function pppmB1234_qp(p1,p2,p3,p4,za,zb) result(pppmB1234_res)
      use quad_precision
      use sprod_qp
      use pppmB34_generic
      use pppmB341_generic
      use pppmB234_generic
      include 'Inc/pppmB1234_inc.f'
      end function pppmB1234_qp

      end module pppmB1234_generic

