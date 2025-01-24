!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module haqgg_pmpm_generic
      implicit none
      public haqgg_pmpm,haqgg_pmpm_qp

      interface haqgg_pmpm
      module procedure haqgg_pmpm,haqgg_pmpm_qp
      end interface

      contains

      subroutine haqgg_pmpm(p1,p2,p3,p4,mtsq,za,zb,amp,Dint,Cint,Bint)
      use double_precision
      use sprod_dp
      include 'Inc/haqgg_pmpm_inc.f'
      end subroutine haqgg_pmpm

      subroutine haqgg_pmpm_qp(p1,p2,p3,p4,mtsq,za,zb,amp,Dint,Cint,Bint)
      use quad_precision
      use sprod_qp
      include 'Inc/haqgg_pmpm_inc.f'
      end subroutine haqgg_pmpm_qp

      end module haqgg_pmpm_generic

