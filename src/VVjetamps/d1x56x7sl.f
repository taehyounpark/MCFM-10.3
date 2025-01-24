!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function xd1x56x7sl(p1,p2,p3,p4,p5,p6,p7,za,zb)
c---   u(1) ubar(2) nu(3) e+(4) e-(5) nubar(6) g(7)
c     d56x234sla=10,formerly box1sl
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp)::xd1x56x7sl
      real(dp):: s234
      include 'zab2.f'

      s234=s(p2,p3)+s(p2,p4)+s(p3,p4)

c--- Box one
c      coeff(4,d56x234sla)
      xd1x56x7sl=0.5_dp*(za(p1,p5)*zb(p2,p4))**2
     & *zab2(p1,p5,p6,p7)*zab2(p7,p5,p6,p1)
     & /(zb(p3,p4)*za(p5,p6)*za(p1,p7)*zab2(p7,p3,p4,p2)*s234)
      return
      end
