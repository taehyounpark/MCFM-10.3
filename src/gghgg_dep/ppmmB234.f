!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module ppmmB234_generic
      implicit none
      public ppmmB234,ppmmB234_qp

      interface ppmmB234
      module procedure ppmmB234,ppmmB234_qp
      end interface

      contains

      function ppmmB234(p1,p2,p3,p4,za,zb) result(ppmmB234_res)
      use double_precision
      use sprod_dp
      include 'Inc/ppmmB234_inc.f'
      end function ppmmB234

      function ppmmB234_qp(p1,p2,p3,p4,za,zb) result(ppmmB234_res)
      use quad_precision
      use sprod_qp
      include 'Inc/ppmmB234_inc.f'
      end function ppmmB234_qp

      end module ppmmB234_generic

