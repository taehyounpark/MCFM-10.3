!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine rat_t9(j1,j2,j3,j4,j5,j6,j7,za,zb,res)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4,j5,j6,j7
      complex(dp):: zab2,res,zab1347,zab2347,zab5347,zab6347,
     & zab7342,zab7345,zab7346,zaba,zbab,zab7156,zab3156,zab4156
      real(dp):: s37p47,s156

c--- statement function
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c--- end statement function
      s37p47=s(j3,j7)+s(j4,j7)
      s156=s(j1,j5)+s(j1,j6)+s(j5,j6)
      zab7342=zab2(j7,j3,j4,j2)
      zab7345=zab2(j7,j3,j4,j5)
      zab7346=zab2(j7,j3,j4,j6)
      zab1347=zab2(j1,j3,j4,j7)
      zab2347=zab2(j2,j3,j4,j7)
      zab5347=zab2(j5,j3,j4,j7)
      zab6347=zab2(j6,j3,j4,j7)
      zab7156=zab2(j7,j1,j5,j6)
      zab3156=zab2(j3,j1,j5,j6)
      zab4156=zab2(j4,j1,j5,j6)

      zaba=zab7345*za(j5,j7)+zab7346*za(j6,j7)
      zbab=zb(j7,j5)*zab5347+zb(j7,j6)*zab6347

      res=za(j1,j7)*zb(j3,j4)*za(j3,j7)**2*za(j5,j7)*zab7346*zbab
     & /(za(j2,j7)*s37p47*zaba**2)
      res=res
     & +za(j3,j4)*zb(j4,j7)**2*zb(j6,j7)*zab1347*zab5347
     & /(s37p47*zab2347*zbab)
      res=res
     & -za(j5,j7)*(za(j1,j7)*za(j3,j4)*zb(j4,j7)**2*zab7346
     & -zb(j3,j4)*za(j3,j7)**2*zb(j6,j7)*zab1347)
     & /(za(j2,j7)*s37p47*zaba)
      res=res
     & +(za(j1,j3)*za(j5,j7)*zb(j4,j7)*zab7346
     & +za(j1,j5)*zb(j3,j4)*za(j3,j7)**2*zb(j6,j7))
     & /(za(j2,j7)*zaba)
      res=res
     & -(za(j1,j7)*zb(j2,j7)*zb(j3,j4)*za(j3,j7)**2*za(j5,j7)
     & *zab2347*zab7346)/(za(j2,j7)**2*s37p47*zab7342*zaba)
      res=res
     & +(zb(j2,j4)*za(j1,j3)*za(j5,j7)*s37p47*zab7346)
     & /(za(j2,j7)*zab7342*zaba)
      res=res
     & -(za(j1,j7)*zb(j2,j7)*s(j3,j4)*za(j3,j7)*zb(j4,j7)*za(j5,j7)
     & *zab7346)/(za(j2,j7)*s37p47*zab7342*zaba)
      res=res
     & +(za(j1,j5)*zb(j2,j7)*zb(j3,j4)*za(j3,j7)**2*zab2347*zab7156)
     & /(za(j2,j7)**2*s156*s37p47*zab7342)
      res=res
     & +za(j1,j5)*zb(j2,j7)*s(j3,j4)*za(j3,j7)*zb(j4,j7)*zab7156
     & /(za(j2,j7)*s156*s37p47*zab7342)
      res=res
     & -za(j1,j5)*zb(j2,j4)*s37p47*zab3156/(za(j2,j7)*s156*zab7342)
      res=res
     & +(za(j1,j5)*za(j3,j4)*zb(j4,j7)**2*zab7156)
     & /(za(j2,j7)*s156*s37p47)
      res=res
     & +za(j1,j5)*za(j3,j4)*zb(j4,j7)**2
     & *(zab4156*zb(j4,j7)+zb(j3,j7)*zab3156)
     & /(s156*s37p47*zab2347)
      res=res
     & -(za(j1,j5)*zb(j4,j7)*zab3156)/(za(j2,j7)*s156)

      res=res/(s(j3,j4)*s(j5,j6))
      return
      end

