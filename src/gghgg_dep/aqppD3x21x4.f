!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqppD3x21x4_generic
      implicit none
      public aqppD3x21x4,aqppD3x21x4_qp

      interface aqppD3x21x4
      module procedure aqppD3x21x4,aqppD3x21x4_qp
      end interface

      contains

      function aqppD3x21x4(p1,p2,p3,p4,mtsq,za,zb) result(aqppD3x21x4_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqppD3x21x4_inc.f'
      end function aqppD3x21x4

      function aqppD3x21x4_qp(p1,p2,p3,p4,mtsq,za,zb) result(aqppD3x21x4_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqppD3x21x4_inc.f'
      end function aqppD3x21x4_qp

      end module aqppD3x21x4_generic

