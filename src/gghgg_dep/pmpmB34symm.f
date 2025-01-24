!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pmpmB34symm_generic
      implicit none
      public pmpmB34symm,pmpmB34symm_qp

      interface pmpmB34symm
      module procedure pmpmB34symm,pmpmB34symm_qp
      end interface

      contains

      function pmpmB34symm(p1,p2,p3,p4,za,zb) result(pmpmB34symm_res)
      use double_precision
      use sprod_dp
      include 'Inc/pmpmB34symm_inc.f'
      end function pmpmB34symm

      function pmpmB34symm_qp(p1,p2,p3,p4,za,zb) result(pmpmB34symm_res)
      use quad_precision
      use sprod_qp
      include 'Inc/pmpmB34symm_inc.f'
      end function pmpmB34symm_qp

      end module pmpmB34symm_generic

