!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function xd2x7x34sl(p1,p2,p3,p4,p5,p6,p7,za,zb)
c---   u(1) ubar(2) nu(3) e+(4) e-(5) nubar(6) g(7)
c     d34x156slb=13, i.e. formerly box5sl
      implicit none
      include'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp)::xd2x7x34sl
      real(dp):: s156,s347
      include 'zaa22.f'

      s347=s(p3,p4)+s(p3,p7)+s(p4,p7)
      s156=s(p1,p5)+s(p1,p6)+s(p5,p6)

c--- Box five
      xd2x7x34sl=s347*s(p2,p7)*za(p1,p5)**2*zab2(p3,p4,p7,p2)**2
     & /(2._dp*s156*zab2(p7,p3,p4,p2)*za(p3,p4)*za(p5,p6)
     & *zaa22(p7,p3,p4,p5,p6,p1))
      return
      end
