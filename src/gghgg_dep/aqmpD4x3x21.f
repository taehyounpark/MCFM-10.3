!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqmpD4x3x21_generic
      implicit none
      public aqmpD4x3x21,aqmpD4x3x21_qp

      interface aqmpD4x3x21
      module procedure aqmpD4x3x21,aqmpD4x3x21_qp
      end interface

      contains

      function aqmpD4x3x21(p1,p2,p3,p4,mtsq,za,zb) result(aqmpD4x3x21_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqmpD4x3x21_inc.f'
      end function aqmpD4x3x21

      function aqmpD4x3x21_qp(p1,p2,p3,p4,mtsq,za,zb) result(aqmpD4x3x21_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqmpD4x3x21_inc.f'
      end function aqmpD4x3x21_qp

      end module aqmpD4x3x21_generic

