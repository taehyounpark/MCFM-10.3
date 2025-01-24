!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqppC3x412_generic
      implicit none
      public aqppC3x412,aqppC3x412_qp

      interface aqppC3x412
      module procedure aqppC3x412,aqppC3x412_qp
      end interface

      contains

      function aqppC3x412(p1,p2,p3,p4,mtsq,za,zb,msqCoeff) result(aqppC3x412_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqppC3x412_inc.f'
      end function aqppC3x412

      function aqppC3x412_qp(p1,p2,p3,p4,mtsq,za,zb,msqCoeff) result(aqppC3x412_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqppC3x412_inc.f'
      end function aqppC3x412_qp

      end module aqppC3x412_generic

