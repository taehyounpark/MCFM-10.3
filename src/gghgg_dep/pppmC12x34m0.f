!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pppmC12x34m0_generic
      implicit none
      public pppmC12x34m0,pppmC12x34m0_qp

      interface pppmC12x34m0
      module procedure pppmC12x34m0,pppmC12x34m0_qp
      end interface

      contains

      function pppmC12x34m0(p1,p2,p3,p4,za,zb) result(pppmC12x34m0_res)
      use double_precision
      use sprod_dp
      include 'Inc/pppmC12x34m0_inc.f'
      end function pppmC12x34m0

      function pppmC12x34m0_qp(p1,p2,p3,p4,za,zb) result(pppmC12x34m0_res)
      use quad_precision
      use sprod_qp
      include 'Inc/pppmC12x34m0_inc.f'
      end function pppmC12x34m0_qp

      end module pppmC12x34m0_generic

