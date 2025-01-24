!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module hgggg_ppmm_generic
      implicit none
      public hgggg_ppmm,hgggg_ppmm_qp

      interface hgggg_ppmm
      module procedure hgggg_ppmm,hgggg_ppmm_qp
      end interface

      contains

      subroutine hgggg_ppmm(p1,p2,p3,p4,mtsq,za,zb,Dcoeff,Ccoeff,Bcoeff,Rat,Cred,I5to4)
      use double_precision
      use sprod_dp
      include 'Inc/hgggg_ppmm_inc.f'
      end subroutine hgggg_ppmm

      subroutine hgggg_ppmm_qp(p1,p2,p3,p4,mtsq,za,zb,Dcoeff,Ccoeff,Bcoeff,Rat,Cred,I5to4)
      use quad_precision
      use sprod_qp
      include 'Inc/hgggg_ppmm_inc.f'
      end subroutine hgggg_ppmm_qp

      end module hgggg_ppmm_generic

