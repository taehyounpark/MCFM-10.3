!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pppmC3x4_generic
      implicit none
      public pppmC3x4,pppmC3x4_qp

      interface pppmC3x4
      module procedure pppmC3x4,pppmC3x4_qp
      end interface

      contains

      function pppmC3x4(p1,p2,p3,p4,za,zb) result(pppmC3x4_res)
      use double_precision
      use sprod_dp
      include 'Inc/pppmC3x4_inc.f'
      end function pppmC3x4

      function pppmC3x4_qp(p1,p2,p3,p4,za,zb) result(pppmC3x4_res)
      use quad_precision
      use sprod_qp
      include 'Inc/pppmC3x4_inc.f'
      end function pppmC3x4_qp

      end module pppmC3x4_generic

