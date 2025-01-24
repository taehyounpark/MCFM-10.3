!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module double_precision
      implicit none
      integer, parameter :: dp=selected_real_kind(precision(1.0d0))
      include 'Inc/consts.f'
      end module double_precision

      module quad_precision
c      use iso_fortran_env
      implicit none
      integer, parameter :: dp=selected_real_kind(2*precision(1.0d0))
      include 'Inc/consts.f'
      end module quad_precision
