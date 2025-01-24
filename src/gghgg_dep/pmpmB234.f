!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pmpmB234_generic
      implicit none
      public pmpmB234,pmpmB234_qp

      interface pmpmB234
      module procedure pmpmB234,pmpmB234_qp
      end interface

      contains

      function pmpmB234(p1,p2,p3,p4,za,zb) result(pmpmB234_res)
      use double_precision
      use sprod_dp
      include 'Inc/pmpmB234_inc.f'
      end function pmpmB234

      function pmpmB234_qp(p1,p2,p3,p4,za,zb) result(pmpmB234_res)
      use quad_precision
      use sprod_qp
      include 'Inc/pmpmB234_inc.f'
      end function pmpmB234_qp

      end module pmpmB234_generic

