!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module hgggg_pppp_generic
      implicit none
      public hgggg_pppp,hgggg_pppp_qp

      interface hgggg_pppp
      module procedure hgggg_pppp,hgggg_pppp_qp
      end interface

      contains

      subroutine hgggg_pppp(p1,p2,p3,p4,mtsq,za,zb,Dcoeff,Ccoeff,Bcoeff,Rat,Cred,I5to4)
      use double_precision
      use sprod_dp
      include 'Inc/hgggg_pppp_inc.f'
      end subroutine hgggg_pppp

      subroutine hgggg_pppp_qp(p1,p2,p3,p4,mtsq,za,zb,Dcoeff,Ccoeff,Bcoeff,Rat,Cred,I5to4)
      use quad_precision
      use sprod_qp
      include 'Inc/hgggg_pppp_inc.f'
      end subroutine hgggg_pppp_qp

      end module hgggg_pppp_generic

