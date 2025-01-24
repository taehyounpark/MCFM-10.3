!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module ppmmB23symm_generic
      implicit none
      public ppmmB23symm,ppmmB23symm_qp

      interface ppmmB23symm
      module procedure ppmmB23symm,ppmmB23symm_qp
      end interface

      contains

      function ppmmB23symm(p1,p2,p3,p4,za,zb) result(ppmmB23symm_res)
      use double_precision
      use sprod_dp
      include 'Inc/ppmmB23symm_inc.f'
      end function ppmmB23symm

      function ppmmB23symm_qp(p1,p2,p3,p4,za,zb) result(ppmmB23symm_res)
      use quad_precision
      use sprod_qp
      include 'Inc/ppmmB23symm_inc.f'
      end function ppmmB23symm_qp

      end module ppmmB23symm_generic

