!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module hgggg_mass_generic
      implicit none
      public hgggg_mass,hgggg_mass_qp

      interface hgggg_mass
      module procedure hgggg_mass,hgggg_mass_qp
      end interface

      contains

      subroutine hgggg_mass(za,zb,msq)
      use double_precision
      use sprod_dp
      include 'Inc/hgggg_mass_inc.f'
      end subroutine hgggg_mass

      subroutine hgggg_mass_qp(za,zb,msq)
      use quad_precision
      use sprod_qp
      include 'Inc/hgggg_mass_inc.f'
      end subroutine hgggg_mass_qp

      end module hgggg_mass_generic

