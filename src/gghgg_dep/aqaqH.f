!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqaqH_generic
      implicit none
      public aqaqH,aqaqH_qp

      interface aqaqH
      module procedure aqaqH,aqaqH_qp
      end interface

      contains

      subroutine aqaqH(p1,p2,p3,p4,mtsq,za,zb,amp,idamp)
      use double_precision
      use sprod_dp
      include 'Inc/aqaqH_inc.f'
      end subroutine aqaqH

      subroutine aqaqH_qp(p1,p2,p3,p4,mtsq,za,zb,amp,idamp)
      use quad_precision
      use sprod_qp
      include 'Inc/aqaqH_inc.f'
      end subroutine aqaqH_qp

      end module aqaqH_generic

