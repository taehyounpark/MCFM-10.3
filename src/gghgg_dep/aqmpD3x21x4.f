!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqmpD3x21x4_generic
      implicit none
      public aqmpD3x21x4,aqmpD3x21x4_qp

      interface aqmpD3x21x4
      module procedure aqmpD3x21x4,aqmpD3x21x4_qp
      end interface

      contains

      function aqmpD3x21x4(p1,p2,p3,p4,mtsq,za,zb) result(aqmpD3x21x4_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqmpD3x21x4_inc.f'
      end function aqmpD3x21x4

      function aqmpD3x21x4_qp(p1,p2,p3,p4,mtsq,za,zb) result(aqmpD3x21x4_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqmpD3x21x4_inc.f'
      end function aqmpD3x21x4_qp

      end module aqmpD3x21x4_generic

