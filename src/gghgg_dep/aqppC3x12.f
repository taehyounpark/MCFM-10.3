!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqppC3x12_generic
      implicit none
      public aqppC3x12,aqppC3x12_qp

      interface aqppC3x12
      module procedure aqppC3x12,aqppC3x12_qp
      end interface

      contains

      function aqppC3x12(p1,p2,p3,p4,za,zb) result(aqppC3x12_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqppC3x12_inc.f'
      end function aqppC3x12

      function aqppC3x12_qp(p1,p2,p3,p4,za,zb) result(aqppC3x12_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqppC3x12_inc.f'
      end function aqppC3x12_qp

      end module aqppC3x12_generic

