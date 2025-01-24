!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine rat_t5(j1,j2,j3,j4,j5,j6,j7,za,zb,res)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4,j5,j6,j7
      complex(dp):: res,zab2,zab3172,zab4172,zab5172,
     & zab2176,zab2173,zab2174,zaba,zbab
      real(dp):: s12p27

c--- statement function
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c--- end statement function
      s12p27=s(j1,j2)+s(j2,j7)
      zab5172=zab2(j5,j1,j7,j2)
      zab4172=zab2(j4,j1,j7,j2)
      zab3172=zab2(j3,j1,j7,j2)
      zab2173=zab2(j2,j1,j7,j3)
      zab2174=zab2(j2,j1,j7,j4)
      zab2176=zab2(j2,j1,j7,j6)
      zaba=zab2174*za(j4,j2)+zab2173*za(j3,j2)
      zbab=zb(j2,j4)*zab4172+zb(j2,j3)*zab3172


      res=zb(j2,j4)*zb(j2,j7)*zb(j2,j6)
     & *(za(j1,j7)*zb(j2,j7)*(za(j2,j3)*zab5172+za(j3,j5)*s12p27)
     & +za(j1,j5)*s12p27*zab3172)/(zb(j1,j2)*zbab)
      res=res
     & +za(j1,j2)*za(j2,j3)
     & *(za(j1,j5)*s12p27*(zb(j2,j4)*zab2176-zb(j6,j4)*s12p27)
     & +za(j1,j7)*za(j2,j5)*zb(j2,j7)*zb(j2,j6)*zab2174)
     & /(za(j2,j7)*zaba)
      res=res
     & -(za(j1,j7)*za(j2,j3)*za(j2,j5)*zb(j2,j7)**2*zab2174*zab2176)
     & /(zb(j1,j7)*za(j2,j7)*zaba)
      res=res
     & +s(j1,j2)*za(j1,j7)*za(j2,j3)*za(j2,j5)*zb(j2,j7)
     & *zab2174*zab2176/(zb(j1,j7)*za(j2,j7)**2*zaba)
      res=res
     & +za(j1,j2)**2*zab2176
     & *(za(j2,j3)*(zb(j2,j4)*zaba*zab5172+za(j2,j5)*zab2174*zbab)
     & +zb(j2,j4)*za(j3,j5)*s12p27*zaba)
     & /(za(j2,j7)*zaba**2)

      res=-res/(2d0*s12p27*s(j3,j4)*s(j5,j6)*za(j1,j7))

      return
      end

