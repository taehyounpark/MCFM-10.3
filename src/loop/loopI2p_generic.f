!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module loopI2p_generic
      implicit none
      public loopI2p,loopI2pc

      interface loopI2p
      module procedure loopI2p,loopI2pc
      end interface

      contains

      function loopI2p(p1,m1,m2,mu2,ep) result(loopI2p_res)
      use avh_olo
      use mod_qcdloop_c
      implicit none
      include 'types.f'
      include 'scalarselect.f'
      real(dp):: p1,m1,m2
      include 'loopI2p_inc.f'
      end function loopI2p
      
      function loopI2pc(p1,m1,m2,mu2,ep) result(loopI2p_res)
      use avh_olo
      use mod_qcdloop_c
      implicit none
      include 'types.f'
      include 'scalarselect.f'
      real(dp)::p1
      complex(dp)::m1,m2

      include 'loopI2pc_inc.f'
      
      return
      end function loopI2pc

      end module loopI2p_generic

