!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function xd1x7x56sl(p1,p2,p3,p4,p5,p6,p7,za,zb)
c---   u(1) ubar(2) nu(3) e+(4) e-(5) nubar(6) g(7)
c     d56x234slb=11,formerly box2sl
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp)::xd1x7x56sl
      real(dp):: s234,s567
      include 'zaa22.f'

      s234=s(p2,p3)+s(p2,p4)+s(p3,p4)
      s567=s(p5,p6)+s(p5,p7)+s(p6,p7)

c      coeff(4,d56x234slb)
      xd1x7x56sl=s567**3*s(p1,p7)*za(p5,p7)**2
     & *zab2(p3,p2,p4,p1)**2
     & /(2._dp*s234*zaa22(p2,p3,p4,p5,p6,p7)*zab2(p7,p5,p6,p1)**3
     & *za(p3,p4)*za(p5,p6))
      return
      end
