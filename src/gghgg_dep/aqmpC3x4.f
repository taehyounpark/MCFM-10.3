!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqmpC3x4_generic
      implicit none
      public aqmpC3x4,aqmpC3x4_qp

      interface aqmpC3x4
      module procedure aqmpC3x4,aqmpC3x4_qp
      end interface

      contains

      function aqmpC3x4(p1,p2,p3,p4,za,zb) result(aqmpC3x4_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqmpC3x4_inc.f'
      end function aqmpC3x4

      function aqmpC3x4_qp(p1,p2,p3,p4,za,zb) result(aqmpC3x4_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqmpC3x4_inc.f'
      end function aqmpC3x4_qp

      end module aqmpC3x4_generic

