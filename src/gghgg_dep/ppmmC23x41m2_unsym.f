!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module ppmmC23x41m2_unsym_generic
      implicit none
      public ppmmC23x41m2_unsym,ppmmC23x41m2_unsym_qp

      interface ppmmC23x41m2_unsym
      module procedure ppmmC23x41m2_unsym,ppmmC23x41m2_unsym_qp
      end interface

      contains

      function ppmmC23x41m2_unsym(p1,p2,p3,p4,za,zb) result(ppmmC23x41m2_unsym_res)
      use double_precision
      use sprod_dp
      include 'Inc/ppmmC23x41m2_unsym_inc.f'
      end function ppmmC23x41m2_unsym

      function ppmmC23x41m2_unsym_qp(p1,p2,p3,p4,za,zb) result(ppmmC23x41m2_unsym_res)
      use quad_precision
      use sprod_qp
      include 'Inc/ppmmC23x41m2_unsym_inc.f'
      end function ppmmC23x41m2_unsym_qp

      end module ppmmC23x41m2_unsym_generic

