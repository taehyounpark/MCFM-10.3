!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function M3bit1(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq.7.5
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: M3bit1,zab2
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      t(j1,j2,j3)=s(j1,j2)+s(j2,j3)+s(j3,j1)
c     End statement functions


      M3bit1=+0.5_dp*za(j1,j5)**2*zb(j1,j4)*zb(j2,j4)*t(j2,j3,j4)
     & /(zb(j2,j3)*zb(j3,j4)*za(j5,j6)*zab2(j1,j3,j4,j2))
     & *(-2*zb(j2,j4)/zab2(j1,j3,j4,j2))

      return
      end



