!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqmpB12_unsym_generic
      implicit none
      public aqmpB12_unsym,aqmpB12_unsym_qp

      interface aqmpB12_unsym
      module procedure aqmpB12_unsym,aqmpB12_unsym_qp
      end interface

      contains

      function aqmpB12_unsym(p1,p2,p3,p4,za,zb) result(aqmpB12_unsym_res)
      use double_precision
      use sprod_dp
      include 'Inc/aqmpB12_unsym_inc.f'
      end function aqmpB12_unsym

      function aqmpB12_unsym_qp(p1,p2,p3,p4,za,zb) result(aqmpB12_unsym_res)
      use quad_precision
      use sprod_qp
      include 'Inc/aqmpB12_unsym_inc.f'
      end function aqmpB12_unsym_qp
      end module aqmpB12_unsym_generic

