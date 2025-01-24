!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqmpC4x123_generic
      implicit none
      public aqmpC4x123,aqmpC4x123_qp

      interface aqmpC4x123
      module procedure aqmpC4x123,aqmpC4x123_qp
      end interface

      contains

      function aqmpC4x123(p1,p2,p3,p4,mtsq,za,zb,msqCoeff) result(aqmpC4x123_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqmpC4x123_inc.f'
      end function aqmpC4x123

      function aqmpC4x123_qp(p1,p2,p3,p4,mtsq,za,zb,msqCoeff) result(aqmpC4x123_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqmpC4x123_inc.f'
      end function aqmpC4x123_qp

      end module aqmpC4x123_generic

