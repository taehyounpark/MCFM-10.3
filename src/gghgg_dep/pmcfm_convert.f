!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pmcfm_convert_generic
      implicit none
      public pmcfm_convert,pmcfm_convert_qp

      interface pmcfm_convert
      module procedure pmcfm_convert,pmcfm_convert_qp
      end interface

      contains

      subroutine pmcfm_convert(p,pmcfm)
      use double_precision
      use sprod_dp
      include 'Inc/pmcfm_convert_inc.f'
      end subroutine pmcfm_convert

      subroutine pmcfm_convert_qp(p,pmcfm)
      use quad_precision
      use sprod_qp
      include 'Inc/pmcfm_convert_inc.f'
      end subroutine pmcfm_convert_qp

      end module pmcfm_convert_generic

