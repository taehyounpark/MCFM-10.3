!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module loopI4_generic
      implicit none
      public loopI4,loopI4c

      interface loopI4
      module procedure loopI4,loopI4c
      end interface

      contains

      function loopI4(p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,mu2,ep) 
     & result(loopI4_res)
      use avh_olo
      use mod_qcdloop_c
      implicit none
      include 'types.f'
      include 'scalarselect.f'
      real(dp):: p1,p2,p3,p4,s12,s23,m1,m2,m3,m4
      include 'loopI4_inc.f'
      end function loopI4
      
      function loopI4c(p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,mu2,ep) 
     & result(loopI4_res)
      use avh_olo
      use mod_qcdloop_c
      implicit none
      include 'types.f'
      include 'scalarselect.f'
      real(dp)::p1,p2,p3,p4,s12,s23
      complex(dp)::m1,m2,m3,m4
      include 'loopI4c_inc.f'
      return
      end function loopI4c

      end module loopI4_generic

