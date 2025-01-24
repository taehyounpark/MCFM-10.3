!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pmpmC12x34_generic
      implicit none
      public pmpmC12x34,pmpmC12x34_qp

      interface pmpmC12x34
      module procedure pmpmC12x34,pmpmC12x34_qp
      end interface

      contains

      function pmpmC12x34(p1,p2,p3,p4,mtsq,za,zb,msqCoeff) result(pmpmC12x34_res)
      use double_precision
      use sprod_dp
      include 'Inc/pmpmC12x34_inc.f'
      end function pmpmC12x34

      function pmpmC12x34_qp(p1,p2,p3,p4,mtsq,za,zb,msqCoeff) result(pmpmC12x34_res)
      use quad_precision
      use sprod_qp
      include 'Inc/pmpmC12x34_inc.f'
      end function pmpmC12x34_qp

      end module pmpmC12x34_generic

