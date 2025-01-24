!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function Afphiqarbmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: Afphiqarbmpmp
      include 'mxpart.f'
      include 'zprods_decl.f'
c---arXIv:09060008v1 Eq.(4.15)
      integer:: j1,j2,j3,j4
      complex(dp):: Afphiqarbmppm
      Afphiqarbmpmp=-Afphiqarbmppm(j1,j2,j4,j3,za,zb)
      return
      end

