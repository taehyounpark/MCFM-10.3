!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function xd1x7x2(p1,p2,p3,p4,p5,p6,p7,za,zb)
c---   u(1) ubar(2) nu(3) e+(4) e-(5) nubar(6) g(7)
c      d127=5    formerly box5
      implicit none
      include'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp)::xd1x7x2
      include 'zaa22.f'

c--- in my notation, corresponds to (s3456 = s127)
c      coeff(4,d127)
      xd1x7x2=0.5_dp*za(p1,p5)*zb(p2,p7)*zb(p1,p7)*zab2(p1,p2,p7,p4)
     & *(za(p2,p3)*za(p1,p5)*zb(p5,p6)+za(p1,p3)*zab2(p2,p4,p3,p6))
     & /(zaa22(p1,p6,p5,p4,p3,p2)*s(p3,p4)*s(p5,p6))
      return
      end
