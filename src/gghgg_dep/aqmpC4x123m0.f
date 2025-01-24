!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqmpC4x123m0_generic
      implicit none
      public aqmpC4x123m0,aqmpC4x123m0_qp

      interface aqmpC4x123m0
      module procedure aqmpC4x123m0,aqmpC4x123m0_qp
      end interface

      contains

      function aqmpC4x123m0(p1,p2,p3,p4,za,zb) result(aqmpC4x123m0_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqmpC4x123m0_inc.f'
      end function aqmpC4x123m0

      function aqmpC4x123m0_qp(p1,p2,p3,p4,za,zb) result(aqmpC4x123m0_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqmpC4x123m0_inc.f'
      end function aqmpC4x123m0_qp

      end module aqmpC4x123m0_generic

