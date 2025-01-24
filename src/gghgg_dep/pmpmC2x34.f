!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pmpmC2x34_generic
      implicit none
      public pmpmC2x34,pmpmC2x34_qp

      interface pmpmC2x34
      module procedure pmpmC2x34,pmpmC2x34_qp
      end interface

      contains

      function pmpmC2x34(p1,p2,p3,p4,za,zb) result(pmpmC2x34_res)
      use double_precision
      use sprod_dp
      include 'Inc/pmpmC2x34_inc.f'
      end function pmpmC2x34

      function pmpmC2x34_qp(p1,p2,p3,p4,za,zb) result(pmpmC2x34_res)
      use quad_precision
      use sprod_qp
      include 'Inc/pmpmC2x34_inc.f'
      end function pmpmC2x34_qp

      end module pmpmC2x34_generic

