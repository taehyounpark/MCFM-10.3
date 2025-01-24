!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pppmB34_generic
      implicit none
      public pppmB34,pppmB34_qp

      interface pppmB34
      module procedure pppmB34,pppmB34_qp
      end interface

      contains

      function pppmB34(p1,p2,p3,p4,za,zb) result(pppmB34_res)
      use double_precision
      use sprod_dp
      include 'Inc/pppmB34_inc.f'
      end function pppmB34

      function pppmB34_qp(p1,p2,p3,p4,za,zb) result(pppmB34_res)
      use quad_precision
      use sprod_qp
      include 'Inc/pppmB34_inc.f'
      end function pppmB34_qp

      end module pppmB34_generic

