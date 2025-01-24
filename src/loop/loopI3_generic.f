!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module loopI3_generic
      implicit none
      public loopI3,loopI3c

      interface loopI3
      module procedure loopI3,loopI3c
      end interface

      contains

      function loopI3(p1,p2,p3,m1,m2,m3,mu2,ep) result(loopI3_res)
      use avh_olo
      use mod_qcdloop_c
      implicit none
      include 'types.f'
      include 'scalarselect.f'
      real(dp):: p1,p2,p3,m1,m2,m3
      include 'loopI3_inc.f'
      end function loopI3
      
      function loopI3c(p1,p2,p3,m1,m2,m3,mu2,ep) result(loopI3_res)
      use avh_olo
      use mod_qcdloop_c
      implicit none
      include 'types.f'
      include 'scalarselect.f'
      real(dp)::p1,p2,p3
      complex(dp)::m1,m2,m3
      include 'loopI3c_inc.f'
      return
      end function loopI3c

      end module loopI3_generic

