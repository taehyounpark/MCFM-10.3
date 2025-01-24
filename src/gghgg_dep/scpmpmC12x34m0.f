!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module scpmpmC12x34m0_generic
      implicit none
      public scpmpmC12x34m0,scpmpmC12x34m0_qp

      interface scpmpmC12x34m0
      module procedure scpmpmC12x34m0,scpmpmC12x34m0_qp
      end interface

      contains

      function scpmpmC12x34m0(p1,p2,p3,p4,za,zb) result(scpmpmC12x34m0_res)
      use double_precision
      use sprod_dp
      include 'Inc/scpmpmC12x34m0_inc.f'
      end function scpmpmC12x34m0

      function scpmpmC12x34m0_qp(p1,p2,p3,p4,za,zb) result(scpmpmC12x34m0_res)
      use quad_precision
      use sprod_qp
      include 'Inc/scpmpmC12x34m0_inc.f'
      end function scpmpmC12x34m0_qp

      end module scpmpmC12x34m0_generic

