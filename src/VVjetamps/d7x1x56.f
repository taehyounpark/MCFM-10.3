!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function xd7x1x56(p1,p2,p3,p4,p5,p6,p7,za,zb)
c---   u(1) ubar(2) nu(3) e+(4) e-(5) nubar(6) g(7)
c     d234x56=1     formerly box1
      implicit none
      include'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp):: xd7x1x56
      real(dp):: s234
      include 'zaa22.f'

      s234=s(p2,p3)+s(p2,p4)+s(p3,p4)

c--- in my notation, corresponds to (d7x1x56)
c      coeff(4,d234x56)=
      xd7x1x56=0.5_dp*s(p1,p7)
     & *(za(p1,p5)**2*zab2(p3,p2,p4,p7)**2
     & /(za(p3,p4)*za(p5,p6)*zaa22(p2,p4,p3,p6,p5,p1)*s234)
     & -(zab2(p3,p6,p5,p1)*zab2(p7,p5,p1,p6))**2
     & /(za(p3,p4)*zb(p5,p6)*za(p2,p7)*zab2(p7,p6,p5,p1)**3))
      return
      end
