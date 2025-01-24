!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pmpmC1x234m0_generic
      implicit none
      public pmpmC1x234m0,pmpmC1x234m0_qp

      interface pmpmC1x234m0
      module procedure pmpmC1x234m0,pmpmC1x234m0_qp
      end interface

      contains

      function pmpmC1x234m0(p1,p2,p3,p4,za,zb) result(pmpmC1x234m0_res)
      use double_precision
      use sprod_dp
      include 'Inc/pmpmC1x234m0_inc.f'
      end function pmpmC1x234m0

      function pmpmC1x234m0_qp(p1,p2,p3,p4,za,zb) result(pmpmC1x234m0_res)
      use quad_precision
      use sprod_qp
      include 'Inc/pmpmC1x234m0_inc.f'
      end function pmpmC1x234m0_qp

      end module pmpmC1x234m0_generic

