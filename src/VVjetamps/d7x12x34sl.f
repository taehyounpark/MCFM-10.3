!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function xd7x12x34sl(p1,p2,p3,p4,p5,p6,p7,za,zb)
c---   u(1) ubar(2) nu(3) e+(4) e-(5) nubar(6) g(7)
c     This is known as d34x12x56slb=16 formerly box12sl
      implicit none
      include'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp):: xd7x12x34sl
      real(dp)::s127,s56,s12,s567
      include 'zaa22.f'

      s127=s(p1,p2)+s(p1,p7)+s(p2,p7)
      s567=s(p5,p6)+s(p5,p7)+s(p6,p7)
      s56=s(p5,p6)
      s12=s(p1,p2)


c      coeff(4,d34x12x56slb)
      xd7x12x34sl=zab2(p7,p1,p2,p4)**2
     & *zaa22(p1,p3,p4,p5,p6,p7)**2*za(p7,p5)**2
     & *(s127*s567-s56*s12)
     & /(2d0*za(p1,p7)*zb(p3,p4)*za(p5,p6)
     & *zaa22(p7,p5,p6,p3,p4,p2)*zaa22(p7,p5,p6,p3,p4,p7)**3)

      return
      end
