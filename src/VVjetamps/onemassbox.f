!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function onemassbox(p1sq,p2sq,p3sq,p4sq,
     & s12,s23,m1sq,m2sq,m3sq,m4sq)
      implicit none
      include'types.f'
      real(dp):: p1sq,p2sq,p3sq,p4sq,s12,s23,m1sq,m2sq,m3sq,m4sq
      complex(dp):: onemassbox,Lsm1
      onemassbox=2d0*lsm1(-s12,-p4sq,-s23,-p4sq)/(s23*s12)
      return
      end
