!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module haqgg_pmmm_generic
      implicit none
      public haqgg_pmmm,haqgg_pmmm_qp

      interface haqgg_pmmm
      module procedure haqgg_pmmm,haqgg_pmmm_qp
      end interface

      contains

      subroutine haqgg_pmmm(p1,p2,p3,p4,mtsq,za,zb,amp,Dint,Cint,Bint)
      use double_precision
      use sprod_dp
      include 'Inc/haqgg_pmmm_inc.f'
      end subroutine haqgg_pmmm

      subroutine haqgg_pmmm_qp(p1,p2,p3,p4,mtsq,za,zb,amp,Dint,Cint,Bint)
      use quad_precision
      use sprod_qp
      include 'Inc/haqgg_pmmm_inc.f'
      end subroutine haqgg_pmmm_qp

      end module haqgg_pmmm_generic

