!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module haqaq_mass_tb_generic
      implicit none
      public haqaq_mass_tb,haqaq_mass_tb_qp

      interface haqaq_mass_tb
      module procedure haqaq_mass_tb,haqaq_mass_tb_qp
      end interface

      contains

      subroutine haqaq_mass_tb(i1,i2,i3,i4,za,zb,msq,idmsq)
      use double_precision
      use sprod_dp
      include 'Inc/haqaq_mass_tb_inc.f'
      end subroutine haqaq_mass_tb

      subroutine haqaq_mass_tb_qp(i1,i2,i3,i4,za,zb,msq,idmsq)
      use quad_precision
      use sprod_qp
      include 'Inc/haqaq_mass_tb_inc.f'
      end subroutine haqaq_mass_tb_qp

      end module haqaq_mass_tb_generic

