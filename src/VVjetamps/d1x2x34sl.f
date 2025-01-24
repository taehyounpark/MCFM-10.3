!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function xd1x2x34sl(p1,p2,p3,p4,p5,p6,p7,za,zb)
c---   u(1) ubar(2) nu(3) e+(4) e-(5) nubar(6) g(7)
c      d34x567sl=14, formerly box9sl
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp)::xd1x2x34sl
      include 'zaa22.f'

c--- Box nine
c      coeff(4,d34x567sl)
      xd1x2x34sl=(s(p1,p2)*zb(p5,p6)*
     & (za(p3,p4)*za(p1,p5)**2*zb(p2,p4)**2*zaa22(p2,p3,p4,p5,p6,p7)
     & *zab2(p7,p5,p6,p1)*zab2(p2,p3,p4,p1)**2
     & +za(p1,p7)*zab2(p7,p3,p4,p2)*zaa22(p2,p3,p4,p6,p7,p5)**2
     & *zab2(p3,p2,p4,p1)**2*zb(p3,p4)))
     & /(2._dp*za(p1,p7)*zab2(p7,p3,p4,p2)
     & *zab2(p7,p5,p6,p1)*zab2(p2,p3,p4,p1)**2
     & *zaa22(p2,p3,p4,p5,p6,p7)*s(p3,p4)*s(p5,p6))
      return
      end
