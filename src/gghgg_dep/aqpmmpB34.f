!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqpmmpB34_generic
      implicit none
      public aqpmmpB34,aqpmmpB34_qp

      interface aqpmmpB34
      module procedure aqpmmpB34,aqpmmpB34_qp
      end interface

      contains

      function aqpmmpB34(p1,p2,p3,p4,za,zb) result(aqpmmpB34_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqpmmpB34_inc.f'
      end function aqpmmpB34

      function aqpmmpB34_qp(p1,p2,p3,p4,za,zb) result(aqpmmpB34_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqpmmpB34_inc.f'
      end function aqpmmpB34_qp

      end module aqpmmpB34_generic

