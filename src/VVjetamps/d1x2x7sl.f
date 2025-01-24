!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function xd1x2x7sl(p1,p2,p3,p4,p5,p6,p7,za,zb)
c---   u(1) ubar(2) nu(3) e+(4) e-(5) nubar(6) g(7)
c     d3456slb=8, formerly box7sl
      implicit none
      include'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp):: xd1x2x7sl
      include 'zaa22.f'

c--- Box seven
c      coeff(4,d3456slb)
      xd1x2x7sl=zab2(p1,p2,p7,p4)*za(p1,p5)*zb(p2,p7)*s(p1,p2)
     &  *(za(p1,p7)*zab2(p3,p5,p1,p6)-za(p1,p3)*za(p2,p7)*zb(p2,p6))
     & /(2._dp*s(p3,p4)*s(p5,p6)*za(p1,p7)*zaa22(p1,p5,p6,p3,p4,p7))
      return
      end
