!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module scpmpmC12x34m0_unsym_generic
      implicit none
      public scpmpmC12x34m0_unsym,scpmpmC12x34m0_unsym_qp

      interface scpmpmC12x34m0_unsym
      module procedure scpmpmC12x34m0_unsym,scpmpmC12x34m0_unsym_qp
      end interface

      contains

      function scpmpmC12x34m0_unsym(p1,p2,p3,p4,za,zb) result(scpmpmC12x34m0_unsym_res)
      use double_precision
      use sprod_dp
      include 'Inc/scpmpmC12x34m0_unsym_inc.f'
      end function scpmpmC12x34m0_unsym

      function scpmpmC12x34m0_unsym_qp(p1,p2,p3,p4,za,zb) result(scpmpmC12x34m0_unsym_res)
      use quad_precision
      use sprod_qp
      include 'Inc/scpmpmC12x34m0_unsym_inc.f'
      end function scpmpmC12x34m0_unsym_qp

      end module scpmpmC12x34m0_unsym_generic

