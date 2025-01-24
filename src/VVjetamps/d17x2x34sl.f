!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function xd17x2x34sl(p1,p2,p3,p4,p5,p6,p7,za,zb)
c---   u(1) ubar(2) nu(3) e+(4) e-(5) nubar(6) g(7)
c      d56x34x17sl=9, formerly box6sl
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp)::xd17x2x34sl
      real(dp):: s234,s127
      include 'zaa22.f'

      s234=s(p2,p3)+s(p2,p4)+s(p3,p4)
      s127=s(p1,p2)+s(p1,p7)+s(p2,p7)

c      coeff(4,d56x34x17sl)=
      xd17x2x34sl=(s127*s234-s(p3,p4)*s(p1,p7))*za(p1,p2)**2
     & *zab2(p2,p3,p4,p6)**2*zaa22(p2,p1,p7,p5,p6,p3)**2
     & /(2._dp*za(p3,p4)*zb(p5,p6)*za(p1,p7)
     & *zaa22(p2,p3,p4,p5,p6,p2)**3*zaa22(p2,p3,p4,p5,p6,p7))
      return
      end

