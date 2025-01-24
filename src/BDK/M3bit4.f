!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function M3bit4(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq.7.5
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: M3bit4,zab2,zab
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c     End statement functions

      M3bit4=-0.5_dp*za(j2,j3)*zb(j2,j4)*zab2(j5,j3,j4,j2)
     & /(zb(j2,j3)*za(j5,j6)*zab2(j1,j3,j4,j2))
     & *(zab(j5,j2,j4))

      return
      end



