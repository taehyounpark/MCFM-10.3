!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      module scaleset_m
          use types
          implicit none

          private
          public :: set_scales
          logical, save, public :: use_resummation_recoil = .false.
          real(dp), save, public :: resummation_recoil = 0._dp

          contains

          subroutine set_scales(renscale_in,facscale_in)
              use constants, only: pi
              implicit none

              include 'couple.f'
              include 'qcdcouple.f'
              include 'nlooprun.f'
              include 'scale.f'
              include 'facscale.f'

              real(dp), intent(in) :: renscale_in, facscale_in

              real(dp) :: alphas

              scale = renscale_in
              facscale = facscale_in

              as = alphas(scale,amz,nlooprun)

              ason2pi = as/2/pi
              ason4pi = as/4/pi
              gsq = 4*pi*as
              musq = scale**2
          end subroutine
      end module
