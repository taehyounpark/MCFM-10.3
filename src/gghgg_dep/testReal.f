!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module testReal_generic
      implicit none
      public testReal,testReal_qp

      interface testReal
      module procedure testReal,testReal_qp
      end interface

      contains

      subroutine testReal(p)
      use double_precision
      use sprod_dp
      include 'Inc/testReal_inc.f'
      end subroutine testReal

      subroutine testReal_qp(p)
      use quad_precision
      use sprod_qp
      include 'Inc/testReal_inc.f'
      end subroutine testReal_qp

      end module testReal_generic

