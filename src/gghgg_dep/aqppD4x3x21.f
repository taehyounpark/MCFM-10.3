!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqppD4x3x21_generic
      implicit none
      public aqppD4x3x21,aqppD4x3x21_qp

      interface aqppD4x3x21
      module procedure aqppD4x3x21,aqppD4x3x21_qp
      end interface

      contains

      function aqppD4x3x21(p1,p2,p3,p4,mtsq,za,zb) result(aqppD4x3x21_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqppD4x3x21_inc.f'
      end function aqppD4x3x21

      function aqppD4x3x21_qp(p1,p2,p3,p4,mtsq,za,zb) result(aqppD4x3x21_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqppD4x3x21_inc.f'
      end function aqppD4x3x21_qp

      end module aqppD4x3x21_generic

