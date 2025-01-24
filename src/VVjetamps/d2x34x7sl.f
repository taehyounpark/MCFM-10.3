!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function xd2x34x7sl(p1,p2,p3,p4,p5,p6,p7,za,zb)
c---   u(1) ubar(2) nu(3) e+(4) e-(5) nubar(6) g(7)
c     d34x156sla=12, formerly box4sl
      implicit none
      include'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp)::xd2x34x7sl
      real(dp):: s156
      include 'zab2.f'

      s156=s(p1,p5)+s(p1,p6)+s(p5,p6)
c      coeff(4,d34x156sla)
      xd2x34x7sl=-0.5_dp*zab2(p2,p3,p4,p7)
     & *(za(p3,p7)*zab2(p7,p3,p4,p2)*zab2(p2,p1,p5,p6))**2
     & /(s156*za(p2,p7)**3*zab2(p7,p3,p4,p2)*zab2(p7,p5,p6,p1)
     & *za(p3,p4)*zb(p5,p6))
      return
      end
