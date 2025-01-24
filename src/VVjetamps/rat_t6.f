!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine rat_t6(j1,j2,j3,j4,j5,j6,j7,za,zb,res)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4,j5,j6,j7
      complex(dp):: res,zab2,zab3561,zab4561,zab7561,
     & zab1564,zab1563,zab1562,zab1567,zaba,zbab,zab3241,
     & zab1274,zab3247,zab3246
      real(dp):: s15p16,s234,s156

c--- statement function
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c--- end statement function
      s15p16=s(j1,j5)+s(j1,j6)
      s234=s(j2,j3)+s(j2,j4)+s(j3,j4)
      s156=s(j1,j5)+s(j1,j6)+s(j5,j6)
      zab3561=zab2(j3,j5,j6,j1)
      zab4561=zab2(j4,j5,j6,j1)
      zab7561=zab2(j7,j5,j6,j1)
      zab1567=zab2(j1,j5,j6,j7)
      zab1562=zab2(j1,j5,j6,j2)
      zab1563=zab2(j1,j5,j6,j3)
      zab1564=zab2(j1,j5,j6,j4)
      zab3241=zab2(j3,j2,j4,j1)
      zab3247=zab2(j3,j2,j4,j7)
      zab3246=zab2(j3,j2,j4,j6)
      zab1274=zab2(j1,j2,j7,j4)

      zbab=+(zb(j1,j4)*zab4561+zb(j1,j3)*zab3561)
      zaba=-(za(j1,j4)*zab1564+za(j1,j3)*zab1563)
      res=
     & (zb(j1,j2)*za(j1,j5)*zb(j1,j6)*zab1564*zab3561
     & -za(j1,j5)*zb(j1,j4)*zb(j2,j6)*s15p16*zab3561
     & -za(j1,j3)*zb(j1,j4)*zb(j1,j6)**2*za(j5,j6)*zab1562
     & +zb(j1,j6)*zb(j2,j4)*za(j3,j5)*s15p16**2)
     & /(za(j1,j7)*zab7561*zbab)
      res=res
     & +(zb(j1,j2)*zb(j1,j4)*zb(j1,j6)**2*za(j5,j6)*zab1567*zab3561)
     & /(zb(j1,j7)*zab7561**2*zbab)
      res=res
     & +(zb(j1,j6)*(za(j1,j5)*zb(j1,j7)*zab1564*zab3561
     & -za(j1,j3)*zb(j1,j4)*zb(j1,j6)*za(j5,j6)*zab1567))
     & /(za(j2,j7)*s15p16*zbab)
      res=res
     & +(za(j1,j2)*zb(j1,j6)
     & *(zb(j1,j2)*za(j1,j5)*zab1564*zab3561
     & -za(j1,j3)*zb(j1,j4)*zb(j1,j6)*za(j5,j6)*zab1562))
     & /(za(j1,j7)*za(j2,j7)*s15p16*zbab)
      res=res
     & +(za(j1,j5)*zb(j1,j4)*zb(j6,j7)*zab3561
     &  -zb(j1,j6)*za(j3,j5)*zb(j4,j7)*s15p16)
     & /(za(j2,j7)*zbab)
      res=res
     & -(za(j1,j2)*(za(j1,j5)*zb(j1,j4)*zb(j2,j6)*zab3561
     & -zb(j1,j6)*zb(j2,j4)*za(j3,j5)*s15p16))
     & /(za(j1,j7)*za(j2,j7)*zbab)
      res=res
     & -(zb(j1,j2)*zb(j1,j4)*zb(j1,j6)**2*za(j5,j6)*zaba*zab3561)
     & /(za(j1,j7)*zab7561*zbab**2)
      res=res
     & -(zb(j1,j4)*zb(j1,j6)**2*zb(j1,j7)*za(j5,j6)*zaba*zab3561)
     & /(za(j2,j7)*s15p16*zbab**2)
      res=res
     & -(za(j1,j2)*zb(j1,j2)*zb(j1,j4)*zb(j1,j6)**2
     & *za(j5,j6)*zaba*zab3561)/(za(j1,j7)*za(j2,j7)*s15p16*zbab**2)
      res=res
     & -(zb(j2,j4)*(za(j1,j5)*zb(j1,j6)*zb(j1,j7)*za(j3,j4)
     & *zab1564*zab7561
     & -za(j1,j5)*zb(j1,j7)*zab7561
     & *(zb(j1,j6)*za(j2,j3)*zab1562+zab3246*s15p16)
     & -zb(j1,j6)**2*za(j1,j7)*zab3241*za(j5,j6)*zab1567))
     & /(za(j1,j7)*zb(j1,j7)*s234*zab7561**2)
      res=res
     & -(za(j1,j3)*za(j1,j5)*zb(j6,j7)*zab1564)
     & /(za(j2,j7)*zaba)
      res=res
     & -(za(j1,j3)*za(j1,j5)*zb(j1,j6)*zab1564*zab1567)
     & /(za(j2,j7)*s15p16*zaba)
      res=res
     & -(za(j1,j2)*za(j1,j3)*za(j1,j5)*zb(j1,j6)*zab1562*zab1564)
     & /(za(j1,j7)*za(j2,j7)*s15p16*zaba)
      res=res
     & +(za(j1,j2)*za(j1,j3)*za(j1,j5)
     & *(zb(j4,j6)*zab1562+za(j1,j5)*zb(j2,j4)*zb(j5,j6)))
     & /(za(j1,j7)*za(j2,j7)*zaba)
      res=res
     & +(2*za(j1,j5)**2*zb(j2,j4)*zab3247*zb(j5,j6))
     & /(za(j1,j7)*s156*s234)
      res=res
     & -(2*za(j1,j5)*zab1274*zb(j5,j6)
     & *(za(j3,j5)*s15p16+za(j1,j3)*za(j6,j5)*zb(j1,j6)))
     & /(za(j1,j7)*za(j2,j7)*s156*s15p16)

      res=res/(2d0*s(j3,j4)*s(j5,j6))
      return
      end

