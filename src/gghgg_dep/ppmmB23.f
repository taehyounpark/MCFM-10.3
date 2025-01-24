!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module ppmmB23_generic
      implicit none
      public ppmmB23,ppmmB23_qp

      interface ppmmB23
      module procedure ppmmB23,ppmmB23_qp
      end interface

      contains

      function ppmmB23(p1,p2,p3,p4,za,zb) result(ppmmB23_res)
      use double_precision
      use sprod_dp
      use ppmmB23symm_generic
      include 'Inc/ppmmB23_inc.f'
      end function ppmmB23

      function ppmmB23_qp(p1,p2,p3,p4,za,zb) result(ppmmB23_res)
      use quad_precision
      use sprod_qp
      use ppmmB23symm_generic
      include 'Inc/ppmmB23_inc.f'
      end function ppmmB23_qp

      end module ppmmB23_generic

