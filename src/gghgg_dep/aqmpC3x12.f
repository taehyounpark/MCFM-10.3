!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqmpC3x12_generic
      implicit none
      public aqmpC3x12,aqmpC3x12_qp

      interface aqmpC3x12
      module procedure aqmpC3x12,aqmpC3x12_qp
      end interface

      contains

      function aqmpC3x12(p1,p2,p3,p4,za,zb) result(aqmpC3x12_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqmpC3x12_inc.f'
      end function aqmpC3x12

      function aqmpC3x12_qp(p1,p2,p3,p4,za,zb) result(aqmpC3x12_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqmpC3x12_inc.f'
      end function aqmpC3x12_qp

      end module aqmpC3x12_generic

