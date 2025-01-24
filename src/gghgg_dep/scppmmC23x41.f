!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module scppmmC23x41_generic
      implicit none
      public scppmmC23x41,scppmmC23x41_qp

      interface scppmmC23x41
      module procedure scppmmC23x41,scppmmC23x41_qp
      end interface

      contains

      function scppmmC23x41(p1,p2,p3,p4,mtsq,za,zb) result(scppmmC23x41_res)
      use double_precision
      use sprod_dp
      include 'Inc/scppmmC23x41_inc.f'
      end function scppmmC23x41

      function scppmmC23x41_qp(p1,p2,p3,p4,mtsq,za,zb) result(scppmmC23x41_res)
      use quad_precision
      use sprod_qp
      include 'Inc/scppmmC23x41_inc.f'
      end function scppmmC23x41_qp

      end module scppmmC23x41_generic

