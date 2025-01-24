!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqmpC12x34m2_generic
      implicit none
      public aqmpC12x34m2,aqmpC12x34m2_qp

      interface aqmpC12x34m2
      module procedure aqmpC12x34m2,aqmpC12x34m2_qp
      end interface

      contains

      function aqmpC12x34m2(p1,p2,p3,p4,za,zb) result(aqmpC12x34m2_res)
      use double_precision
      use sprod_dp
      use aqmpC12x34m2unsym_generic
      include 'Inc/aqmpC12x34m2_inc.f'
      end function aqmpC12x34m2

      function aqmpC12x34m2_qp(p1,p2,p3,p4,za,zb) result(aqmpC12x34m2_res)
      use quad_precision
      use sprod_qp
      use aqmpC12x34m2unsym_generic
      include 'Inc/aqmpC12x34m2_inc.f'
      end function aqmpC12x34m2_qp

      end module aqmpC12x34m2_generic

