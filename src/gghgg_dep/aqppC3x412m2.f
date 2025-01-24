!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqppC3x412m2_generic
      implicit none
      public aqppC3x412m2,aqppC3x412m2_qp

      interface aqppC3x412m2
      module procedure aqppC3x412m2,aqppC3x412m2_qp
      end interface

      contains

      function aqppC3x412m2(p1,p2,p3,p4,za,zb) result(aqppC3x412m2_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqppC3x412m2_inc.f'
      end function aqppC3x412m2

      function aqppC3x412m2_qp(p1,p2,p3,p4,za,zb) result(aqppC3x412m2_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqppC3x412m2_inc.f'
      end function aqppC3x412m2_qp

      end module aqppC3x412m2_generic

