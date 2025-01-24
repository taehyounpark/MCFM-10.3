!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqppC4x123_generic
      implicit none
      public aqppC4x123,aqppC4x123_qp

      interface aqppC4x123
      module procedure aqppC4x123,aqppC4x123_qp
      end interface

      contains

      function aqppC4x123(p1,p2,p3,p4,mtsq,za,zb,msqCoeff) result(aqppC4x123_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqppC4x123_inc.f'
      end function aqppC4x123

      function aqppC4x123_qp(p1,p2,p3,p4,mtsq,za,zb,msqCoeff) result(aqppC4x123_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqppC4x123_inc.f'
      end function aqppC4x123_qp

      end module aqppC4x123_generic

