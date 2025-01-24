!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pppmC2x34_generic
      implicit none
      public pppmC2x34,pppmC2x34_qp

      interface pppmC2x34
      module procedure pppmC2x34,pppmC2x34_qp
      end interface

      contains

      function pppmC2x34(p1,p2,p3,p4,za,zb) result(pppmC2x34_res)
      use double_precision
      use sprod_dp
      include 'Inc/pppmC2x34_inc.f'
      end function pppmC2x34

      function pppmC2x34_qp(p1,p2,p3,p4,za,zb) result(pppmC2x34_res)
      use quad_precision
      use sprod_qp
      include 'Inc/pppmC2x34_inc.f'
      end function pppmC2x34_qp

      end module pppmC2x34_generic

