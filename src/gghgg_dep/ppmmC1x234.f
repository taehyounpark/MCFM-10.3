!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module ppmmC1x234_generic
      implicit none
      public ppmmC1x234,ppmmC1x234_qp

      interface ppmmC1x234
      module procedure ppmmC1x234,ppmmC1x234_qp
      end interface

      contains

      function ppmmC1x234(p1,p2,p3,p4,mtsq,za,zb,msqCoeff) result(ppmmC1x234_res)
      use double_precision
      use sprod_dp
      include 'Inc/ppmmC1x234_inc.f'
      end function ppmmC1x234

      function ppmmC1x234_qp(p1,p2,p3,p4,mtsq,za,zb,msqCoeff) result(ppmmC1x234_res)
      use quad_precision
      use sprod_qp
      include 'Inc/ppmmC1x234_inc.f'
      end function ppmmC1x234_qp

      end module ppmmC1x234_generic

