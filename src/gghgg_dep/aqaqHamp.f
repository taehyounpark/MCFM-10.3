!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module aqaqHamp_generic
      implicit none
      public aqaqHamp,aqaqHamp_qp

      interface aqaqHamp
      module procedure aqaqHamp,aqaqHamp_qp
      end interface

      contains

      subroutine aqaqHamp(p1,p2,p3,p4,mtsq,za,zb,amp)
      use double_precision
      use sprod_dp
      include 'Inc/aqaqHamp_inc.f'
      end subroutine aqaqHamp

      subroutine aqaqHamp_qp(p1,p2,p3,p4,mtsq,za,zb,amp)
      use quad_precision
      use sprod_qp
      include 'Inc/aqaqHamp_inc.f'
      end subroutine aqaqHamp_qp

      end module aqaqHamp_generic

