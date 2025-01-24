!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function xd7x12x56sl(p1,p2,p3,p4,p5,p6,p7,za,zb)
c---   u(1) ubar(2) nu(3) e+(4) e-(5) nubar(6) g(7)
c     d12x56x34sl=18, formerly box13sl
      implicit none
      include'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp):: xd7x12x56sl
      real(dp):: s347,s127,s12,s34
      include 'zaa22.f'

      s127=s(p1,p2)+s(p1,p7)+s(p2,p7)
      s347=s(p3,p4)+s(p3,p7)+s(p4,p7)
      s34=s(p3,p4)
      s12=s(p1,p2)

c--- Box thirteen
c      coeff(4,d12x56x34sl)=

      xd7x12x56sl=za(p1,p7)**2*(s347*s127-s12*s34)
     & *zab2(p7,p3,p4,p6)**2*zaa22(p7,p1,p2,p5,p6,p3)**2
     & /(2._dp*za(p3,p4)*zb(p5,p6)*za(p2,p7)
     & *zaa22(p7,p3,p4,p1,p2,p7)**3*zaa22(p7,p3,p4,p5,p6,p1))

      return
      end
