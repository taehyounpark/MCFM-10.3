!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function A0Hggggpppp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0Hggggpppp
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      complex(dp):: A0phiggggpppp,A0phiggggmmmm
      A0Hggggpppp=A0phiggggpppp(j1,j2,j3,j4,za,zb)
     &           +A0phiggggmmmm(j1,j2,j3,j4,zb,za)
      return
      end
