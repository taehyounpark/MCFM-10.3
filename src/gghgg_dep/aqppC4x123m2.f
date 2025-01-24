!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqppC4x123m2_generic
      implicit none
      public aqppC4x123m2,aqppC4x123m2_qp

      interface aqppC4x123m2
      module procedure aqppC4x123m2,aqppC4x123m2_qp
      end interface

      contains

      function aqppC4x123m2(p1,p2,p3,p4,za,zb) result(aqppC4x123m2_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqppC4x123m2_inc.f'
      end function aqppC4x123m2

      function aqppC4x123m2_qp(p1,p2,p3,p4,za,zb) result(aqppC4x123m2_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqppC4x123m2_inc.f'
      end function aqppC4x123m2_qp

      end module aqppC4x123m2_generic

