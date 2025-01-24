!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqpmmpB123_generic
      implicit none
      public aqpmmpB123,aqpmmpB123_qp

      interface aqpmmpB123
      module procedure aqpmmpB123,aqpmmpB123_qp
      end interface

      contains

      function aqpmmpB123(p1,p2,p3,p4,za,zb) result(aqpmmpB123_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqpmmpB123_inc.f'
      end function aqpmmpB123

      function aqpmmpB123_qp(p1,p2,p3,p4,za,zb) result(aqpmmpB123_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqpmmpB123_inc.f'
      end function aqpmmpB123_qp

      end module aqpmmpB123_generic

