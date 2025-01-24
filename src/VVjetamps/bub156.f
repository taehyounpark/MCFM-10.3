!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function bub156(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer p1,p2,p3,p4,p5,p6,p7
      real(dp):: s156,s7x156,s15p16
      complex(dp):: res,bub156,zab7165,zab7561,zab3561
c--- statement functions
      include 'zab2.f'
c--- end statement functions
      s15p16=s(p5,p1)+s(p6,p1)
      s156=s(p5,p6)+s(p5,p1)+s(p6,p1)
      s7x156=s(p5,p7)+s(p6,p7)+s(p1,p7)
      zab7165=zab2(p7,p1,p6,p5)
      zab7561=zab2(p7,p5,p6,p1)
      zab3561=zab2(p3,p5,p6,p1)
      res=
     & za(p6,p1)*zab3561**2*zb(p6,p5)
     & *(zab7165/zab7561-za(p6,p1)*zb(p6,p5)/(2._dp*s15p16))
     & /(s156*s15p16*zab7561)
     & +(za(p3,p7)**2*zab7165*zb(p7,p5)*
     & (zab3561/(za(p3,p7)*zab7561)
     & +(s156*zb(p7,p5)/(2._dp*zab7165)-zab2(p3,p2,p4,p7)/za(p3,p7))
     & /s7x156))/(zab7561*s7x156)

      bub156=res/(za(p2,p7)*za(p4,p3)*zb(p6,p5))

      return
      end

