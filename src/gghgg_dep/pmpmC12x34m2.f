!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pmpmC12x34m2_generic
      implicit none
      public pmpmC12x34m2,pmpmC12x34m2_qp

      interface pmpmC12x34m2
      module procedure pmpmC12x34m2,pmpmC12x34m2_qp
      end interface

      contains

      function pmpmC12x34m2(p1,p2,p3,p4,za,zb) result(pmpmC12x34m2_res)
      use double_precision
      use sprod_dp
      include 'Inc/pmpmC12x34m2_inc.f'
      end function pmpmC12x34m2

      function pmpmC12x34m2_qp(p1,p2,p3,p4,za,zb) result(pmpmC12x34m2_res)
      use quad_precision
      use sprod_qp
      include 'Inc/pmpmC12x34m2_inc.f'
      end function pmpmC12x34m2_qp

      end module pmpmC12x34m2_generic

