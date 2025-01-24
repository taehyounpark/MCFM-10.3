!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqpmmpB412_generic
      implicit none
      public aqpmmpB412,aqpmmpB412_qp

      interface aqpmmpB412
      module procedure aqpmmpB412,aqpmmpB412_qp
      end interface

      contains

      function aqpmmpB412(p1,p2,p3,p4,za,zb) result(aqpmmpB412_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqpmmpB412_inc.f'
      end function aqpmmpB412

      function aqpmmpB412_qp(p1,p2,p3,p4,za,zb) result(aqpmmpB412_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqpmmpB412_inc.f'
      end function aqpmmpB412_qp

      end module aqpmmpB412_generic

