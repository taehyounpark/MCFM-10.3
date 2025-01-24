!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module haqgg_mass_generic
      implicit none
      public haqgg_mass,haqgg_mass_qp

      interface haqgg_mass
      module procedure haqgg_mass,haqgg_mass_qp
      end interface

      contains

      subroutine haqgg_mass(i1,i2,i3,i4,za,zb,msq)
      use double_precision
      use sprod_dp
      include 'Inc/haqgg_mass_inc.f'
      end subroutine haqgg_mass

      subroutine haqgg_mass_qp(i1,i2,i3,i4,za,zb,msq)
      use quad_precision
      use sprod_qp
      include 'Inc/haqgg_mass_inc.f'
      end subroutine haqgg_mass_qp

      end module haqgg_mass_generic

