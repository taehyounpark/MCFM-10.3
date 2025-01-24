!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine rat_t8(j1,j2,j3,j4,j5,j6,j7,za,zb,res)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4,j5,j6,j7
      complex(dp)::res,zab2,zaba,zbab,zab3246,zab3247,zaba73,
     & zab7561,zab7562,zab7563,zab7564,zab1567,zab3567,zab4567
      real(dp):: s57p67,s234

c--- statement function
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c--- end statement function
      s57p67=s(j5,j7)+s(j6,j7)
      s234=s(j2,j3)+s(j2,j4)+s(j3,j4)
      zab7561=zab2(j7,j5,j6,j1)
      zab7562=zab2(j7,j5,j6,j2)
      zab7563=zab2(j7,j5,j6,j3)
      zab7564=zab2(j7,j5,j6,j4)
      zab1567=zab2(j1,j5,j6,j7)
      zab3567=zab2(j3,j5,j6,j7)
      zab4567=zab2(j4,j5,j6,j7)
      zab3247=zab2(j3,j2,j4,j7)
      zab3246=zab2(j3,j2,j4,j6)
      zaba=za(j4,j7)*zab7564+za(j3,j7)*zab7563
      zbab=zb(j7,j4)*zab4567+zb(j7,j3)*zab3567
      zaba73=zab7562*za(j2,j3)+zab7564*za(j4,j3)

      res=(za(j3,j7)*zb(j5,j6)*za(j5,j7)**2*zab7562*zab7564*zbab)
     & /(s57p67*zab7561*zaba**2)

     & +(zb(j2,j7)*zb(j4,j7)*za(j5,j6)*zb(j6,j7)**2*zab1567*zab3567)
     & /(zb(j1,j7)*s57p67*zab1567*zbab)

     & -((zab7564*zb(j6,j7)-zb(j4,j7)*zb(j5,j6)*za(j5,j7))
     &  *za(j5,j6)*za(j3,j7)*zb(j6,j7)*zab7562
     &  -zb(j2,j7)*zb(j5,j6)*za(j5,j7)**2*zab3567*zab7564)
     & /(s57p67*zab7561*zaba)

     & -za(j1,j5)*za(j3,j7)*za(j5,j7)*zb(j5,j6)*zab7564*zb(j2,j7)
     & /(za(j1,j7)*zab7561*zaba)

     & -zb(j1,j7)*za(j3,j7)*zb(j5,j6)*za(j5,j7)**2
     & *zab1567*zab7562*zab7564/(za(j1,j7)*s57p67*zab7561**2*zaba)

     & +zb(j2,j4)*za(j3,j5)*zb(j5,j6)*za(j5,j7)*s57p67/(zab7561*zaba)

     & +za(j1,j5)*zb(j2,j4)*(zaba73*zb(j6,j7)-zab3246*s57p67)
     & /(za(j1,j7)*s234*zab7561)

     & -zb(j2,j4)*za(j5,j6)*zb(j6,j7)
     & *(zaba73*zb(j6,j7)+zab3247*zb(j5,j6)*za(j5,j7))
     & /(s234*s57p67*zab7561)

     & -zb(j1,j7)*zb(j2,j4)*zb(j5,j6)*za(j5,j7)**2*zab1567*zaba73
     & /(za(j1,j7)*s234*s57p67*zab7561**2)

     & +zb(j2,j4)*za(j5,j6)*zb(j6,j7)**2*zab1567*zab3247
     & /(zb(j1,j7)*s234*s57p67*zab1567)

      res=res/(s(j3,j4)*s(j5,j6))

      return
      end

