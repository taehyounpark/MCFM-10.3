!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      ! initialize additional calculated parameters, maybe eft-modified
      module eftcouple
          use types
          implicit none

          public

          include 'ewcouple.f'
          include 'qcdcouple.f'
          include 'masses.f'

          real(dp), public, save :: ecossin
          real(dp), public, save :: eftgw, gb

          public :: eftcouple_init

          contains

          subroutine eftcouple_init()
              implicit none

              gb = sqrt(esq)/sqrt(1-xw)
              ecossin = sqrt(gw**2 + gb**2)

          end subroutine

      end module
