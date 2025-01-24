!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function A0Hggggpmmm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0Hggggpmmm

      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      complex(dp):: A0phiggggpmmm
      A0Hggggpmmm=A0phiggggpmmm(j1,j2,j3,j4,za,zb)
c     &           +A0phiggggmppp(j1,j2,j3,j4,zb,za) ! This term is zero
      return
      end
