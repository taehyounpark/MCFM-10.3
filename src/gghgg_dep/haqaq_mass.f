!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module haqaq_mass_generic
      implicit none
      public haqaq_mass,haqaq_mass_qp

      interface haqaq_mass
      module procedure haqaq_mass,haqaq_mass_qp
      end interface

      contains

      subroutine haqaq_mass(i1,i2,i3,i4,za,zb,msq,idmsq)
      use double_precision
      use sprod_dp
      include 'Inc/haqaq_mass_inc.f'
      end subroutine haqaq_mass

      subroutine haqaq_mass_qp(i1,i2,i3,i4,za,zb,msq,idmsq)
      use quad_precision
      use sprod_qp
      include 'Inc/haqaq_mass_inc.f'
      end subroutine haqaq_mass_qp

      end module haqaq_mass_generic

