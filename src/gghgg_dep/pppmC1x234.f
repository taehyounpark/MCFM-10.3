!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pppmC1x234_generic
      implicit none
      public pppmC1x234,pppmC1x234_qp

      interface pppmC1x234
      module procedure pppmC1x234,pppmC1x234_qp
      end interface

      contains

      function pppmC1x234(p1,p2,p3,p4,mtsq,za,zb,msqCoeff) result(pppmC1x234_res)
      use double_precision
      use sprod_dp
      include 'Inc/pppmC1x234_inc.f'
      end function pppmC1x234

      function pppmC1x234_qp(p1,p2,p3,p4,mtsq,za,zb,msqCoeff) result(pppmC1x234_res)
      use quad_precision
      use sprod_qp
      include 'Inc/pppmC1x234_inc.f'
      end function pppmC1x234_qp

      end module pppmC1x234_generic

