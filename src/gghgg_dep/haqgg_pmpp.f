!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module haqgg_pmpp_generic
      implicit none
      public haqgg_pmpp,haqgg_pmpp_qp

      interface haqgg_pmpp
      module procedure haqgg_pmpp,haqgg_pmpp_qp
      end interface

      contains

      subroutine haqgg_pmpp(p1,p2,p3,p4,mtsq,za,zb,amp,Dint,Cint,Bint)
      use double_precision
      use sprod_dp
      include 'Inc/haqgg_pmpp_inc.f'
      end subroutine haqgg_pmpp

      subroutine haqgg_pmpp_qp(p1,p2,p3,p4,mtsq,za,zb,amp,Dint,Cint,Bint)
      use quad_precision
      use sprod_qp
      include 'Inc/haqgg_pmpp_inc.f'
      end subroutine haqgg_pmpp_qp

      end module haqgg_pmpp_generic

