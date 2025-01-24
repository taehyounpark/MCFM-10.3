!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module fillformfactor_generic
      implicit none
      public fillformfactor,fillformfactor_qp

      interface fillformfactor
      module procedure fillformfactor,fillformfactor_qp
      end interface

      contains

      subroutine fillformfactor(p1,p2,p3,p4,mtsq,FL,FT)
      use double_precision
      use sprod_dp
      include 'Inc/fillformfactor_inc.f'
      end subroutine fillformfactor

      subroutine fillformfactor_qp(p1,p2,p3,p4,mtsq,FL,FT)
      use quad_precision
      use sprod_qp
      include 'Inc/fillformfactor_qp_inc.f'
      end subroutine fillformfactor_qp

      end module fillformfactor_generic

