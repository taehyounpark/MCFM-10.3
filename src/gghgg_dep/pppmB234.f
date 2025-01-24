!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pppmB234_generic
      implicit none
      public pppmB234,pppmB234_qp

      interface pppmB234
      module procedure pppmB234,pppmB234_qp
      end interface

      contains

      function pppmB234(p1,p2,p3,p4,za,zb) result(pppmB234_res)
      use double_precision
      use sprod_dp
      include 'Inc/pppmB234_inc.f'
      end function pppmB234

      function pppmB234_qp(p1,p2,p3,p4,za,zb) result(pppmB234_res)
      use quad_precision
      use sprod_qp
      include 'Inc/pppmB234_inc.f'
      end function pppmB234_qp

      end module pppmB234_generic

