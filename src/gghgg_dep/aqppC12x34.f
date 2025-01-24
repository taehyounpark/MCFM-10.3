!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqppC12x34_generic
      implicit none
      public aqppC12x34,aqppC12x34_qp

      interface aqppC12x34
      module procedure aqppC12x34,aqppC12x34_qp
      end interface

      contains

      function aqppC12x34(p1,p2,p3,p4,mtsq,za,zb,msqCoeff) result(aqppC12x34_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqppC12x34_inc.f'
      end function aqppC12x34

      function aqppC12x34_qp(p1,p2,p3,p4,mtsq,za,zb,msqCoeff) result(aqppC12x34_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqppC12x34_inc.f'
      end function aqppC12x34_qp

      end module aqppC12x34_generic

