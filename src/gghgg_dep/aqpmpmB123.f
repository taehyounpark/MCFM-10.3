!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqpmpmB123_generic
      implicit none
      public aqpmpmB123,aqpmpmB123_qp

      interface aqpmpmB123
      module procedure aqpmpmB123,aqpmpmB123_qp
      end interface

      contains

      function aqpmpmB123(p1,p2,p3,p4,za,zb) result(aqpmpmB123_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqpmpmB123_inc.f'
      end function aqpmpmB123

      function aqpmpmB123_qp(p1,p2,p3,p4,za,zb) result(aqpmpmB123_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqpmpmB123_inc.f'
      end function aqpmpmB123_qp

      end module aqpmpmB123_generic

