!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function A1Hggggmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A1Hggggmpmp
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      complex(dp):: A1phiggggmpmp
c--- 0704.3914v3 Eqs. (2.4) and (2.6)
c--- Note: c.c. is equivalent to interchanging za and zb
      A1Hggggmpmp=A1phiggggmpmp(j1,j2,j3,j4,za,zb)
     &           +A1phiggggmpmp(j2,j3,j4,j1,zb,za)
      return
      end
