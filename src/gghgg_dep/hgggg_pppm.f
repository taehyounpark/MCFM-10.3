!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module hgggg_pppm_generic
      implicit none
      public hgggg_pppm,hgggg_pppm_qp

      interface hgggg_pppm
      module procedure hgggg_pppm,hgggg_pppm_qp
      end interface

      contains

      subroutine hgggg_pppm(p1,p2,p3,p4,mtsq,za,zb,Dcoeff,Ccoeff,Bcoeff,Rat,Cred,I5to4)
      use double_precision
      use sprod_dp
      include 'Inc/hgggg_pppm_inc.f'
      end subroutine hgggg_pppm

      subroutine hgggg_pppm_qp(p1,p2,p3,p4,mtsq,za,zb,Dcoeff,Ccoeff,Bcoeff,Rat,Cred,I5to4)
      use quad_precision
      use sprod_qp
      include 'Inc/hgggg_pppm_inc.f'
      end subroutine hgggg_pppm_qp

      end module hgggg_pppm_generic

