!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module ppmmC23x41_generic
      implicit none
      public ppmmC23x41,ppmmC23x41_qp

      interface ppmmC23x41
      module procedure ppmmC23x41,ppmmC23x41_qp
      end interface

      contains

      function ppmmC23x41(p1,p2,p3,p4,mtsq,za,zb,msqCoeff) result(ppmmC23x41_res)
      use double_precision
      use sprod_dp
      include 'Inc/ppmmC23x41_inc.f'
      end function ppmmC23x41

      function ppmmC23x41_qp(p1,p2,p3,p4,mtsq,za,zb,msqCoeff) result(ppmmC23x41_res)
      use quad_precision
      use sprod_qp
      include 'Inc/ppmmC23x41_inc.f'
      end function ppmmC23x41_qp

      end module ppmmC23x41_generic

