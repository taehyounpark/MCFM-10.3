!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function xd2x17x56(p1,p2,p3,p4,p5,p6,p7,za,zb)
c---   u(1) ubar(2) nu(3) e+(4) e-(5) nubar(6) g(7)
c     d34x56x27=4   formerly box4
      implicit none
      include'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp):: xd2x17x56
      real(dp):: s234,s127,s34,s17
      include 'zaa22.f'

      s127=s(p1,p2)+s(p1,p7)+s(p2,p7)
      s234=s(p2,p3)+s(p2,p4)+s(p3,p4)
      s34=s(p3,p4)
      s17=s(p1,p7)
c      coeff(4,d34x56x27)
      xd2x17x56=0.5_dp*za(p1,p2)**3
     & *(zab2(p2,p3,p4,p6)*zaa22(p2,p1,p7,p5,p6,p3))**2
     & *(s127*s234-s17*s34)/(za(p4,p3)*zb(p5,p6)*za(p1,p7)*za(p2,p7)
     & *zaa22(p2,p3,p4,p5,p6,p1)*zaa22(p2,p3,p4,p5,p6,p2)**3)

      return
      end
