!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pmpmC3x4_generic
      implicit none
      public pmpmC3x4,pmpmC3x4_qp

      interface pmpmC3x4
      module procedure pmpmC3x4,pmpmC3x4_qp
      end interface

      contains

      function pmpmC3x4(p1,p2,p3,p4,za,zb) result(pmpmC3x4_res)
      use double_precision
      use sprod_dp
      include 'Inc/pmpmC3x4_inc.f'
      end function pmpmC3x4

      function pmpmC3x4_qp(p1,p2,p3,p4,za,zb) result(pmpmC3x4_res)
      use quad_precision
      use sprod_qp
      include 'Inc/pmpmC3x4_inc.f'
      end function pmpmC3x4_qp

      end module pmpmC3x4_generic

