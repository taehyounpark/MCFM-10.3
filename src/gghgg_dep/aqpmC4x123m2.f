!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqpmC4x123m2_generic
      implicit none
      public aqpmC4x123m2,aqpmC4x123m2_qp

      interface aqpmC4x123m2
      module procedure aqpmC4x123m2,aqpmC4x123m2_qp
      end interface

      contains

      function aqpmC4x123m2(p1,p2,p3,p4,za,zb) result(aqpmC4x123m2_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqpmC4x123m2_inc.f'
      end function aqpmC4x123m2

      function aqpmC4x123m2_qp(p1,p2,p3,p4,za,zb) result(aqpmC4x123m2_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqpmC4x123m2_inc.f'
      end function aqpmC4x123m2_qp

      end module aqpmC4x123m2_generic

