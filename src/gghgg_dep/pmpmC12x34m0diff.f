!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pmpmC12x34m0diff_generic
      implicit none
      public pmpmC12x34m0diff,pmpmC12x34m0diff_qp

      interface pmpmC12x34m0diff
      module procedure pmpmC12x34m0diff,pmpmC12x34m0diff_qp
      end interface

      contains

      function pmpmC12x34m0diff(p1,p2,p3,p4,za,zb) result(pmpmC12x34m0diff_res)
      use double_precision
      use sprod_dp
      include 'Inc/pmpmC12x34m0diff_inc.f'
      end function pmpmC12x34m0diff

      function pmpmC12x34m0diff_qp(p1,p2,p3,p4,za,zb) result(pmpmC12x34m0diff_res)
      use quad_precision
      use sprod_qp
      include 'Inc/pmpmC12x34m0diff_inc.f'
      end function pmpmC12x34m0diff_qp

      end module pmpmC12x34m0diff_generic

