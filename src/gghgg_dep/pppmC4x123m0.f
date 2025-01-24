!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pppmC4x123m0_generic
      implicit none
      public pppmC4x123m0,pppmC4x123m0_qp

      interface pppmC4x123m0
      module procedure pppmC4x123m0,pppmC4x123m0_qp
      end interface

      contains

      function pppmC4x123m0(p1,p2,p3,p4,za,zb) result(pppmC4x123m0_res)
      use double_precision
      use sprod_dp
      include 'Inc/pppmC4x123m0_inc.f'
      end function pppmC4x123m0

      function pppmC4x123m0_qp(p1,p2,p3,p4,za,zb) result(pppmC4x123m0_res)
      use quad_precision
      use sprod_qp
      include 'Inc/pppmC4x123m0_inc.f'
      end function pppmC4x123m0_qp

      end module pppmC4x123m0_generic

