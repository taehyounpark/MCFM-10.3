!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function A0phiggggmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0phiggggmpmp

c----Expresssion of Eq. (2.14) of hep-ph/0804.4149v3
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      A0phiggggmpmp=za(j1,j3)**4
     &             /(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j1))
      return
      end

