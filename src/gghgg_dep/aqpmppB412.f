!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqpmppB412_generic
      implicit none
      public aqpmppB412,aqpmppB412_qp

      interface aqpmppB412
      module procedure aqpmppB412,aqpmppB412_qp
      end interface

      contains

      function aqpmppB412(p1,p2,p3,p4,za,zb) result(aqpmppB412_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqpmppB412_inc.f'
      end function aqpmppB412

      function aqpmppB412_qp(p1,p2,p3,p4,za,zb) result(aqpmppB412_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqpmppB412_inc.f'
      end function aqpmppB412_qp

      end module aqpmppB412_generic

