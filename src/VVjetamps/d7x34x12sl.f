!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function xd7x34x12sl(p1,p2,p3,p4,p5,p6,p7,za,zb)
c---   u(1) ubar(2) nu(3) e+(4) e-(5) nubar(6) g(7)
c     d34x12x56sla=15, formerly box10sl
      implicit none
      include'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp):: xd7x34x12sl
      real(dp):: s567,s56,s34,s347
      include 'zaa22.f'

      s347=s(p3,p4)+s(p3,p7)+s(p4,p7)
      s567=s(p5,p6)+s(p5,p7)+s(p6,p7)
      s56=s(p5,p6)
      s34=s(p3,p4)

c--- Box ten
c      coeff(4,d34x12x56sla)
      xd7x34x12sl=za(p3,p7)**2*zab2(p7,p5,p6,p2)**2
     & *zaa22(p5,p1,p2,p3,p4,p7)**2*(s347*s567-s34*s56)
     &  /2._dp/za(p3,p4)/za(p5,p6)
     &  /zab2(p7,p5,p6,p1)/zab2(p7,p3,p4,p2)/zaa22(p7,p3,p4,p5,p6,p7)**3

      return
      end
