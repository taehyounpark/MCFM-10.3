!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module hgggg_pmpm_generic
      implicit none
      public hgggg_pmpm,hgggg_pmpm_qp

      interface hgggg_pmpm
      module procedure hgggg_pmpm,hgggg_pmpm_qp
      end interface

      contains

      subroutine hgggg_pmpm(p1,p2,p3,p4,mtsq,za,zb,Dcoeff,Ccoeff,Bcoeff,Rat,Cred,I5to4)
      use double_precision
      use sprod_dp
      include 'Inc/hgggg_pmpm_inc.f'
      end subroutine hgggg_pmpm

      subroutine hgggg_pmpm_qp(p1,p2,p3,p4,mtsq,za,zb,Dcoeff,Ccoeff,Bcoeff,Rat,Cred,I5to4)
      use quad_precision
      use sprod_qp
      include 'Inc/hgggg_pmpm_inc.f'
      end subroutine hgggg_pmpm_qp

      end module hgggg_pmpm_generic

