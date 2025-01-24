!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqppC4x123m0_generic
      implicit none
      public aqppC4x123m0,aqppC4x123m0_qp

      interface aqppC4x123m0
      module procedure aqppC4x123m0,aqppC4x123m0_qp
      end interface

      contains

      function aqppC4x123m0(p1,p2,p3,p4,za,zb) result(aqppC4x123m0_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqppC4x123m0_inc.f'
      end function aqppC4x123m0

      function aqppC4x123m0_qp(p1,p2,p3,p4,za,zb) result(aqppC4x123m0_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqppC4x123m0_inc.f'
      end function aqppC4x123m0_qp

      end module aqppC4x123m0_generic

