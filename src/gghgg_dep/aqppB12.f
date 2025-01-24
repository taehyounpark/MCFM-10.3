!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqppB12_generic
      implicit none
      public aqppB12,aqppB12_qp

      interface aqppB12
      module procedure aqppB12,aqppB12_qp
      end interface

      contains

      function aqppB12(p1,p2,p3,p4,za,zb) result(aqppB12_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqppB12_inc.f'
      end function aqppB12

      function aqppB12_qp(p1,p2,p3,p4,za,zb) result(aqppB12_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqppB12_inc.f'
      end function aqppB12_qp

      end module aqppB12_generic

