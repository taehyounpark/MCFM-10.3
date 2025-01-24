!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function t2(j1,j2,j3,j4)
      implicit none
      include 'types.f'
      complex(dp):: t2
      include 'mxpart.f'
      include 'zprods_com.f'
      integer:: j1,j2,j3,j4
      t2=zb(j1,j2)*za(j4,j2)+zb(j1,j3)*za(j4,j3)
      return
      end

