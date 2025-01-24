!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqpmmmB123_generic
      implicit none
      public aqpmmmB123,aqpmmmB123_qp

      interface aqpmmmB123
      module procedure aqpmmmB123,aqpmmmB123_qp
      end interface

      contains

      function aqpmmmB123(p1,p2,p3,p4,za,zb) result(aqpmmmB123_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqpmmmB123_inc.f'
      end function aqpmmmB123

      function aqpmmmB123_qp(p1,p2,p3,p4,za,zb) result(aqpmmmB123_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqpmmmB123_inc.f'
      end function aqpmmmB123_qp

      end module aqpmmmB123_generic

