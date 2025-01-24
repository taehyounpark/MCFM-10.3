!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module setreal_mcfm_generic
      implicit none
      public setreal_mcfm,setreal_mcfm_qp

      interface setreal_mcfm
      module procedure setreal_mcfm,setreal_mcfm_qp
      end interface

      contains

      subroutine setreal_mcfm(p,rflav,amp2)
      use double_precision
      use sprod_dp
      include 'Inc/setreal_mcfm_inc.f'
      end subroutine setreal_mcfm

      subroutine setreal_mcfm_qp(p,rflav,amp2)
      use quad_precision
      use sprod_qp
      include 'Inc/setreal_mcfm_inc.f'
      end subroutine setreal_mcfm_qp

      end module setreal_mcfm_generic

