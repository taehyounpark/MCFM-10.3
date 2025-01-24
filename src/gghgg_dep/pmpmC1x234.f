!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pmpmC1x234_generic
      implicit none
      public pmpmC1x234,pmpmC1x234_qp

      interface pmpmC1x234
      module procedure pmpmC1x234,pmpmC1x234_qp
      end interface

      contains

      function pmpmC1x234(p1,p2,p3,p4,mtsq,za,zb,msqCoeff) result(pmpmC1x234_res)
      use double_precision
      use sprod_dp
      include 'Inc/pmpmC1x234_inc.f'
      end function pmpmC1x234

      function pmpmC1x234_qp(p1,p2,p3,p4,mtsq,za,zb,msqCoeff) result(pmpmC1x234_res)
      use quad_precision
      use sprod_qp
      include 'Inc/pmpmC1x234_inc.f'
      end function pmpmC1x234_qp

      end module pmpmC1x234_generic

