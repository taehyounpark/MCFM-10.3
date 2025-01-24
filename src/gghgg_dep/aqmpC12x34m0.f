!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqmpC12x34m0_generic
      implicit none
      public aqmpC12x34m0,aqmpC12x34m0_qp

      interface aqmpC12x34m0
      module procedure aqmpC12x34m0,aqmpC12x34m0_qp
      end interface

      contains

      function aqmpC12x34m0(p1,p2,p3,p4,za,zb) result(aqmpC12x34m0_res)
      use double_precision
      use sprod_dp
      use aqmpC12x34m0unsym_generic
      include 'Inc/aqmpC12x34m0_inc.f'
      end function aqmpC12x34m0

      function aqmpC12x34m0_qp(p1,p2,p3,p4,za,zb) result(aqmpC12x34m0_res)
      use quad_precision
      use sprod_qp
      use aqmpC12x34m0unsym_generic
      include 'Inc/aqmpC12x34m0_inc.f'
      end function aqmpC12x34m0_qp

      end module aqmpC12x34m0_generic

