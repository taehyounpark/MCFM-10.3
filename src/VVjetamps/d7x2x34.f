!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function xd7x2x34(p1,p2,p3,p4,p5,p6,p7,za,zb)
c---   u(1) ubar(2) nu(3) e+(4) e-(5) nubar(6) g(7)
c     d34x156=3  formerly box3
      implicit none
      include'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp)::xd7x2x34
      real(dp):: s156
      include 'zaa22.f'

      s156=s(p1,p5)+s(p1,p6)+s(p5,p6)

c--- in my notation, corresponds to (s34,s156)
c      coeff(4,d34x156)=(
      xd7x2x34=0.5_dp*za(p1,p5)**2*s(p2,p7)
     & *(zb(p2,p4)**2/(zb(p3,p4)*za(p5,p6)*za(p1,p7)*zab2(p7,p4,p3,p2))
     & -zab2(p3,p2,p4,p7)**2
     & /(za(p3,p4)*za(p5,p6)*zaa22(p1,p6,p5,p4,p3,p2)*s156))
      return
      end
