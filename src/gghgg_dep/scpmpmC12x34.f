!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module scpmpmC12x34_generic
      implicit none
      public scpmpmC12x34,scpmpmC12x34_qp

      interface scpmpmC12x34
      module procedure scpmpmC12x34,scpmpmC12x34_qp
      end interface

      contains

      function scpmpmC12x34(p1,p2,p3,p4,mtsq,za,zb) result(scpmpmC12x34_res)
      use double_precision
      use sprod_dp
      include 'Inc/scpmpmC12x34_inc.f'
      end function scpmpmC12x34

      function scpmpmC12x34_qp(p1,p2,p3,p4,mtsq,za,zb) result(scpmpmC12x34_res)
      use quad_precision
      use sprod_qp
      include 'Inc/scpmpmC12x34_inc.f'
      end function scpmpmC12x34_qp

      end module scpmpmC12x34_generic

