!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pmpmB34_generic
      implicit none
      public pmpmB34,pmpmB34_qp

      interface pmpmB34
      module procedure pmpmB34,pmpmB34_qp
      end interface

      contains

      function pmpmB34(p1,p2,p3,p4,za,zb) result(pmpmB34_res)
      use double_precision
      use sprod_dp
      use pmpmB34symm_generic
      include 'Inc/pmpmB34_inc.f'
      end function pmpmB34

      function pmpmB34_qp(p1,p2,p3,p4,za,zb) result(pmpmB34_res)
      use quad_precision
      use sprod_qp
      use pmpmB34symm_generic
      include 'Inc/pmpmB34_inc.f'
      end function pmpmB34_qp

      end module pmpmB34_generic

