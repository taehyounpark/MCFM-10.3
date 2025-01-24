!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module loopI2_generic
      implicit none
      public loopI2,loopI2c

      interface loopI2
      module procedure loopI2,loopI2c
      end interface

      contains

      function loopI2(p1,m1,m2,mu2,ep) result(loopI2_res)
      use avh_olo
      use mod_qcdloop_c
      implicit none
      include 'types.f'
      include 'scalarselect.f'
      real(dp):: p1,m1,m2
      include 'loopI2_inc.f'
      end function loopI2
      
      function loopI2c(p1,m1,m2,mu2,ep) result(loopI2_res)
      use avh_olo
      use mod_qcdloop_c
      implicit none
      include 'types.f'
      include 'scalarselect.f'
      real(dp)::p1
      complex(dp)::m1,m2

      include 'loopI2c_inc.f'
      
      return
      end function loopI2c

      end module loopI2_generic

