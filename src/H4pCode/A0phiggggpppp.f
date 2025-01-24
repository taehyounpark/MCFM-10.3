!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function A0phiggggpppp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0phiggggpppp

c----Expresssion of Eq. (4) of hep-ph/0607139v2
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      A0phiggggpppp=czip
      return
      end

