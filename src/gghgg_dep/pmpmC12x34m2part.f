!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pmpmC12x34m2part_generic
      implicit none
      public pmpmC12x34m2part,pmpmC12x34m2part_qp

      interface pmpmC12x34m2part
      module procedure pmpmC12x34m2part,pmpmC12x34m2part_qp
      end interface

      contains

      function pmpmC12x34m2part(p1,p2,p3,p4,za,zb) result(pmpmC12x34m2part_res)
      use double_precision
      use sprod_dp
      include 'Inc/pmpmC12x34m2part_inc.f'
      end function pmpmC12x34m2part

      function pmpmC12x34m2part_qp(p1,p2,p3,p4,za,zb) result(pmpmC12x34m2part_res)
      use quad_precision
      use sprod_qp
      include 'Inc/pmpmC12x34m2part_inc.f'
      end function pmpmC12x34m2part_qp

      end module pmpmC12x34m2part_generic

