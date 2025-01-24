!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module loopI1_generic
      implicit none
      public loopI1,loopI1c

      interface loopI1
      module procedure loopI1,loopI1c
      end interface

      contains

      function loopI1(m1,mu2,ep) result(loopI1_res)
      use avh_olo
      use mod_qcdloop_c
      implicit none
      include 'types.f'
      include 'scalarselect.f'
      real(dp):: m1
      include 'loopI1_inc.f'
      end function loopI1
      
      function loopI1c(m1,mu2,ep) result(loopI1_res)
      use avh_olo
      use mod_qcdloop_c
      implicit none
      include 'types.f'
      include 'scalarselect.f'
      complex(dp)::m1

      include 'loopI1c_inc.f'
      
      return
      end function loopI1c

      end module loopI1_generic

