!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pppmC1x43_generic
      implicit none
      public pppmC1x43,pppmC1x43_qp

      interface pppmC1x43
      module procedure pppmC1x43,pppmC1x43_qp
      end interface

      contains

      function pppmC1x43(p1,p2,p3,p4,za,zb) result(pppmC1x43_res)
      use double_precision
      use sprod_dp
      include 'Inc/pppmC1x43_inc.f'
      end function pppmC1x43

      function pppmC1x43_qp(p1,p2,p3,p4,za,zb) result(pppmC1x43_res)
      use quad_precision
      use sprod_qp
      include 'Inc/pppmC1x43_inc.f'
      end function pppmC1x43_qp

      end module pppmC1x43_generic

