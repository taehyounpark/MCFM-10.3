!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function xd1x27x34sl(p1,p2,p3,p4,p5,p6,p7,za,zb)
c---   u(1) ubar(2) nu(3) e+(4) e-(5) nubar(6) g(7)
c     d34x56x27sl=6, formerly box3sl

      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp)::xd1x27x34sl
      real(dp):: s127,s27,s56,S156
      include 'zaa22.f'

      s127=s(p1,p2)+s(p1,p7)+s(p2,p7)
      s156=s(p1,p5)+s(p1,p6)+s(p5,p6)
      s27=s(p2,p7)
      s56=s(p5,p6)
c--- Box three

c      coeff(4,d34x56x27sl)
      xd1x27x34sl=-0.5_dp*(za(p1,p5)*zab2(p1,p2,p7,p4))**2
     & *(s127*s156-s56*s27)
     & /(zb(p3,p4)*za(p5,p6)*za(p2,p7)
     & *zaa22(p1,p5,p6,p3,p4,p7)*zaa22(p1,p5,p6,p3,p4,p1))

      return
      end
