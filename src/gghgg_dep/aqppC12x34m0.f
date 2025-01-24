!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqppC12x34m0_generic
      implicit none
      public aqppC12x34m0,aqppC12x34m0_qp

      interface aqppC12x34m0
      module procedure aqppC12x34m0,aqppC12x34m0_qp
      end interface

      contains

      function aqppC12x34m0(p1,p2,p3,p4,za,zb) result(aqppC12x34m0_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqppC12x34m0_inc.f'
      end function aqppC12x34m0

      function aqppC12x34m0_qp(p1,p2,p3,p4,za,zb) result(aqppC12x34m0_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqppC12x34m0_inc.f'
      end function aqppC12x34m0_qp

      end module aqppC12x34m0_generic

