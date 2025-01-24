!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      module bbfrac_m
        use types
        implicit none

        private
        include 'maxd.f'

        ! 40 = maxd
        real(dp), save, public :: bbfrac(0:maxd)
!$omp threadprivate(bbfrac)

      end module
