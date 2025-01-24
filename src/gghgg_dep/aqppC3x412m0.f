!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqppC3x412m0_generic
      implicit none
      public aqppC3x412m0,aqppC3x412m0_qp

      interface aqppC3x412m0
      module procedure aqppC3x412m0,aqppC3x412m0_qp
      end interface

      contains

      function aqppC3x412m0(p1,p2,p3,p4,za,zb) result(aqppC3x412m0_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqppC3x412m0_inc.f'
      end function aqppC3x412m0

      function aqppC3x412m0_qp(p1,p2,p3,p4,za,zb) result(aqppC3x412m0_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqppC3x412m0_inc.f'
      end function aqppC3x412m0_qp

      end module aqppC3x412m0_generic

