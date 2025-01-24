!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function A0Hggggmmpp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0Hggggmmpp

      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      complex(dp):: A0phiggggmmpp
      A0Hggggmmpp=A0phiggggmmpp(j1,j2,j3,j4,za,zb)
     &           +A0phiggggmmpp(j3,j4,j1,j2,zb,za)
      return
      end
