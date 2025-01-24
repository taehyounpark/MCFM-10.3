!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pmpmC1x234m2_generic
      implicit none
      public pmpmC1x234m2,pmpmC1x234m2_qp

      interface pmpmC1x234m2
      module procedure pmpmC1x234m2,pmpmC1x234m2_qp
      end interface

      contains

      function pmpmC1x234m2(p1,p2,p3,p4,za,zb) result(pmpmC1x234m2_res)
      use double_precision
      use sprod_dp
      include 'Inc/pmpmC1x234m2_inc.f'
      end function pmpmC1x234m2

      function pmpmC1x234m2_qp(p1,p2,p3,p4,za,zb) result(pmpmC1x234m2_res)
      use quad_precision
      use sprod_qp
      include 'Inc/pmpmC1x234m2_inc.f'
      end function pmpmC1x234m2_qp

      end module pmpmC1x234m2_generic

