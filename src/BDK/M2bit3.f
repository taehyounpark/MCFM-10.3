!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function M2bit3(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq.7.3
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: M2bit3,zab2,zab
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
      t(j1,j2,j3)=s(j1,j2)+s(j2,j3)+s(j3,j1)
c     End statement functions

      M2bit3=
     & -zab(j5,j4,j1)*zab2(j5,j1,j2,j3)*t(j1,j2,j3)
     & /(zb(j2,j3)*za(j5,j6)*zab2(j4,j1,j2,j3)**2)
      return
      end
