!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function xd2x1x7sl(p1,p2,p3,p4,p5,p6,p7,za,zb)
c---   u(1) ubar(2) nu(3) e+(4) e-(5) nubar(6) g(7)
c     d3456sla=7, formerly box8sl
      implicit none
      include'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp):: xd2x1x7sl
      include 'zaa22.f'

c--- Box eight
c      coeff(4,d3456sla)
      xd2x1x7sl=
     & -za(p1,p2)**2*s(p1,p2)*zb(p1,p7)*za(p5,p7)*zab2(p7,p2,p1,p4)
     &  *(za(p3,p7)*za(p2,p1)*zb(p1,p6)+za(p2,p7)*zab2(p3,p5,p7,p6))
     & /(2d0*s(p3,p4)*s(p5,p6)*za(p2,p7)**3*zaa22(p2,p4,p3,p6,p5,p7))
      return
      end

