!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine rat_t11(j1,j2,j3,j4,j5,j6,j7,za,zb,res)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4,j5,j6,j7
      complex(dp)::zab2,res,zab7156,zab7342,zab7561,zab3247,
     & zab2347,zab7234,zab5167,zab1567,zab3246,zab3156,zbab,zaba,
     & zababis
      real(dp):: s234,s156,s156x7

c--- statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c--- end statement functions

      zab7561=zab2(j7,j5,j6,j1)
      zab7156=zab2(j7,j1,j5,j6)
      zab3156=zab2(j3,j1,j5,j6)
      zab5167=zab2(j5,j1,j6,j7)
      zab1567=zab2(j1,j5,j6,j7)
      zab7342=zab2(j7,j3,j4,j2)
      zab3246=zab2(j3,j2,j4,j6)
      zab3247=zab2(j3,j2,j4,j7)
      zab2347=zab2(j2,j3,j4,j7)
      zab7234=zab2(j7,j2,j3,j4)
      zbab=zb(j6,j5)*zab5167+zb(j6,j1)*zab1567
      zaba=zab7342*za(j2,j3)+zab7234*za(j4,j3)
      zababis=zab3247*za(j7,j2)-zab2347*za(j7,j3)
      s234=s(j2,j3)+s(j2,j4)+s(j3,j4)
      s156=s(j1,j5)+s(j1,j6)+s(j5,j6)
      s156x7=s(j1,j7)+s(j5,j7)+s(j6,j7)

      res=zb(j2,j4)*za(j5,j7)*zb(j6,j7)*zab1567*zaba
     & /(za(j1,j7)*s156x7*s234*zab7561)
      res=res+za(j5,j7)*zb(j2,j7)*zab7156*zab7234*zababis
     & /(za(j2,j7)**2*s156x7*zab7342*zab7561)
      res=res
     & +za(j1,j5)*zb(j2,j4)*za(j3,j7)*s156x7*zab7156
     & /(za(j1,j7)*za(j2,j7)*zab7342*zab7561)
      res=res
     & -zb(j2,j4)*za(j3,j5)*s156x7*zab7156/(za(j2,j7)*zab7342*zab7561)
      res=res
     & -za(j3,j7)*za(j5,j7)*zb(j6,j7)*zab1567*zab7234
     & /(za(j1,j7)*za(j2,j7)*s156x7*zab7561)
      res=res
     & +zb(j2,j4)*zab3247*zab5167*zab7156
     & /(s156x7*s234*zab7561)
      res=res
     & +za(j3,j7)*zb(j4,j7)*zab5167*zab7156/(za(j2,j7)*s156x7*zab7561)
      res=res
     & -za(j1,j5)*zb(j2,j4)*zab3246*s156x7/(za(j1,j7)*s234*zab7561)
      res=res
     & -za(j1,j5)*za(j3,j7)*zb(j4,j6)*s156x7
     & /(za(j1,j7)*za(j2,j7)*zab7561)
      res=res
     & -zb(j1,j7)*zb(j2,j4)*za(j5,j7)*zab1567*zab7156*zaba
     & /(za(j1,j7)*s156x7*s234*zab7561**2)
      res=res
     & +zb(j1,j7)*za(j3,j7)*za(j5,j7)*zab1567*zab7156*zab7234
     & /(za(j1,j7)*za(j2,j7)*s156x7*zab7561**2)
      res=res
     & +za(j1,j5)*zb(j2,j7)*zab7234*zab7156*zababis
     & /(za(j2,j7)**2*s156*s156x7*zab7342)
      res=res
     & +za(j1,j5)*zb(j2,j4)*zab3156*s156x7/(za(j2,j7)*s156*zab7342)
      res=res
     & +za(j1,j5)*zb(j4,j7)*zab3247*zbab/(s156*s156x7*zab2347)
      res=res
     & +2*za(j1,j5)*zb(j2,j4)*zab3247*zbab/(s156*s156x7*s234)
      res=res
     & +za(j1,j5)*za(j3,j7)*zb(j4,j7)*zbab/(za(j2,j7)*s156*s156x7)
      res=res
     & +zb(j4,j7)*zb(j6,j7)*zab3247*zab5167/(zb(j1,j7)*s156x7*zab2347)
      res=res
     & +zb(j2,j4)*zab3247*zb(j6,j7)*zab5167/(zb(j1,j7)*s156x7*s234)

      res=-res/(2d0*s(j4,j3)*s(j6,j5))
      return
      end

