!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function bub156sl(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer p2,p4,p3,p5,p6,p1,p7
      complex(dp):: res,bub156sl,
     & zab2561,zab2165,zab7561,zab7165,zab3247
      real(dp)::s15p16,s456,s2x156,s7x156

c--- statement functions
      include 'zab.f'
c--- end statement functions
      zab2561=zab(p2,p5,p1)+zab(p2,p6,p1)
      zab2165=(zab(p2,p6,p5)+zab(p2,p1,p5))
      zab3247=(zab(p3,p2,p7)+zab(p3,p4,p7))
      zab7561=(zab(p7,p5,p1) + zab(p7,p6,p1))
      zab7165=(zab(p7,p6,p5) + zab(p7,p1,p5))
      s15p16=s(p5,p1)+s(p6,p1)
      s456=s(p5,p6)+s(p5,p1)+s(p6,p1)
      s2x156=s(p2,p5)+s(p2,p6)+s(p2,p1)
      s7x156=s(p5,p7)+s(p6,p7)+s(p1,p7)
      res=-zab2165*zab(p1,p6,p5)*za(p2,p3)**2
     & /(s456*zab2561*za(p2,p1)*za(p2,p7)**2)

     & +zab(p1,p6,p5)**2*za(p2,p3)**2
     & /(s456*s15p16*za(p2,p1)*za(p2,p7)**2)

     & +zab(p1,p6,p5)**2*za(p3,p1)*za(p3,p7)
     & /(s456*s15p16*za(p2,p7)*za(p1,p7)**2)

     & -zab(p1,p6,p5)*zab7165*za(p3,p1)*za(p3,p7)
     & /(s456*zab7561*za(p2,p7)*za(p1,p7)**2)

     & -zab(p1,p6,p5)**2*za(p2,p3)*za(p3,p7)
     & /(s456*s15p16*za(p2,p7)**2*za(p1,p7))

     & +zab(p1,p6,p5)*zab7165*za(p2,p3)*za(p3,p7)
     & /(s456*zab7561*za(p2,p7)**2*za(p1,p7))

     & +zab2165*za(p2,p3)**2*zb(p5,p2)
     & /(zab2561*s2x156*za(p2,p7)**2)

     & -zab(p1,p6,p5)*za(p3,p7)**2*zb(p1,p5)
     & /(zab7561**2*za(p2,p7)*za(p1,p7))

     & -zab7165*za(p2,p3)*za(p3,p7)*zb(p7,p5)
     & /(zab7561*s7x156*za(p2,p7)**2)

     & +zab3247*zab7165*za(p3,p7)*zb(p7,p5)
     & /(zab7561*s7x156**2*za(p2,p7))

     & +s456*za(p3,p7)**2*zb(p1,p5)*zb(p7,p5)
     & /(zab7561**2*s7x156*za(p2,p7))

      bub156sl=-res/(za(p4,p3)*zb(p5,p6))


      return
      end

