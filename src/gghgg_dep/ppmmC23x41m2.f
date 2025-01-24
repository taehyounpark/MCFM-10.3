!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module ppmmC23x41m2_generic
      implicit none
      public ppmmC23x41m2,ppmmC23x41m2_qp

      interface ppmmC23x41m2
      module procedure ppmmC23x41m2,ppmmC23x41m2_qp
      end interface

      contains

      function ppmmC23x41m2(p1,p2,p3,p4,za,zb) result(ppmmC23x41m2_res)
      use double_precision
      use sprod_dp
      include 'Inc/ppmmC23x41m2_inc.f'
      end function ppmmC23x41m2

      function ppmmC23x41m2_qp(p1,p2,p3,p4,za,zb) result(ppmmC23x41m2_res)
      use quad_precision
      use sprod_qp
      include 'Inc/ppmmC23x41m2_inc.f'
      end function ppmmC23x41m2_qp

      end module ppmmC23x41m2_generic

