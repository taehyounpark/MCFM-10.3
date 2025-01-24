!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqpmppB123_generic
      implicit none
      public aqpmppB123,aqpmppB123_qp

      interface aqpmppB123
      module procedure aqpmppB123,aqpmppB123_qp
      end interface

      contains

      function aqpmppB123(p1,p2,p3,p4,za,zb) result(aqpmppB123_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqpmppB123_inc.f'
      end function aqpmppB123

      function aqpmppB123_qp(p1,p2,p3,p4,za,zb) result(aqpmppB123_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqpmppB123_inc.f'
      end function aqpmppB123_qp

      end module aqpmppB123_generic

