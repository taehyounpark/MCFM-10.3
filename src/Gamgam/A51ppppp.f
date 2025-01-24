!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function A51ppppp(j1,j2,j3,j4,j5,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      complex(dp):: A51ppppp
      integer:: j1,j2,j3,j4,j5
      complex(dp):: asym,sym
      asym=+zb(j1,j2)*za(j2,j3)*zb(j3,j4)*za(j4,j1)
     &     -za(j1,j2)*zb(j2,j3)*za(j3,j4)*zb(j4,j1)
      sym=cplx1(s(j1,j2)*s(j2,j3)+s(j2,j3)*s(j3,j4)
     &          +s(j3,j4)*s(j4,j5)+s(j4,j5)*s(j5,j1)+s(j5,j1)*s(j1,j2))
c----Eq.(4) of hep-ph/9302280v1 of BDK multiplied by 16*pi^2*(i)
c----in order to match the normalization of A51mmppp.f and A51mpmpp.f
c---- (note additional minus sign from the relation A[1/2]=-A[1])
c--- It therefore gives (16*pi^2)*(-i)*A^{[1/2]}_{5;1}
      A51ppppp=-(sym+asym)
     & /(6._dp*za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j5)*za(j5,j1))
      return
      end
