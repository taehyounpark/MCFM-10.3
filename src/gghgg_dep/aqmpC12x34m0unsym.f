!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqmpC12x34m0unsym_generic
      implicit none
      public aqmpC12x34m0unsym,aqmpC12x34m0unsym_qp

      interface aqmpC12x34m0unsym
      module procedure aqmpC12x34m0unsym,aqmpC12x34m0unsym_qp
      end interface

      contains

      function aqmpC12x34m0unsym(p1,p2,p3,p4,za,zb) result(aqmpC12x34m0unsym_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqmpC12x34m0unsym_inc.f'
      end function aqmpC12x34m0unsym

      function aqmpC12x34m0unsym_qp(p1,p2,p3,p4,za,zb) result(aqmpC12x34m0unsym_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqmpC12x34m0unsym_inc.f'
      end function aqmpC12x34m0unsym_qp

      end module aqmpC12x34m0unsym_generic

