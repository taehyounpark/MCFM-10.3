!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function xd2x1x56sl(p1,p2,p3,p4,p5,p6,p7,za,zb)
c---   u(1) ubar(2) nu(3) e+(4) e-(5) nubar(6) g(7)
c     This is d56x347sl=17 formerly box11sl
      implicit none
      include'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp):: xd2x1x56sl
      include 'zaa22.f'

c--- Box eleven
c      coeff(4,d56x347sl)
      xd2x1x56sl=s(p1,p2)*zb(p3,p4)/(2._dp*s(p3,p4)*s(p5,p6))
     & *((-za(p1,p5)**2*zab2(p3,p4,p7,p2)**2*zb(p5,p6))
     & /(zab2(p7,p3,p4,p2)*zaa22(p1,p5,p6,p3,p4,p7))

     & -zab2(p3,p5,p6,p1)**2*zab2(p2,p1,p5,p6)**2*za(p5,p6)
     & /(zab2(p7,p5,p6,p1)*za(p2,p7)*zab2(p2,p5,p6,p1)**2))

      return
      end
