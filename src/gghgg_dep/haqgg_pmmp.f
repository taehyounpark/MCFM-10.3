!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module haqgg_pmmp_generic
      implicit none
      public haqgg_pmmp,haqgg_pmmp_qp

      interface haqgg_pmmp
      module procedure haqgg_pmmp,haqgg_pmmp_qp
      end interface

      contains

      subroutine haqgg_pmmp(p1,p2,p3,p4,mtsq,za,zb,amp,Dint,Cint,Bint)
      use double_precision
      use sprod_dp
      include 'Inc/haqgg_pmmp_inc.f'
      end subroutine haqgg_pmmp

      subroutine haqgg_pmmp_qp(p1,p2,p3,p4,mtsq,za,zb,amp,Dint,Cint,Bint)
      use quad_precision
      use sprod_qp
      include 'Inc/haqgg_pmmp_inc.f'
      end subroutine haqgg_pmmp_qp

      end module haqgg_pmmp_generic

