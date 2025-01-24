!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine rat_t1(j1,j2,j3,j4,j5,j6,j7,za,zb,res)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4,j5,j6,j7
      complex(dp):: res,zab2
      real(dp):: s234

c--- statement function
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      s234=s(j2,j3)+s(j2,j4)+s(j3,j4)

      res=-zb(j1,j6)*zb(j2,j4)*za(j5,j7)*zab2(j1,j5,j6,j7)
     & *zab2(j3,j2,j4,j1)/(zab2(j7,j5,j6,j1)**2*s234)
     & +za(j1,j5)*zb(j2,j4)*zb(j6,j7)*zab2(j3,j2,j4,j1)
     & /(zab2(j7,j5,j6,j1)*s234)
     & +zb(j1,j4)*zb(j1,j6)*za(j3,j7)*za(j5,j7)*zab2(j1,j5,j6,j7)
     & /(za(j2,j7)*zab2(j7,j5,j6,j1)**2)
     & +za(j1,j2)*zb(j1,j4)*zb(j1,j6)*zb(j2,j7)*za(j3,j7)*za(j5,j7)
     & /(zb(j1,j2)*za(j2,j7)**2*zab2(j7,j5,j6,j1))

     & +(
     &   +za(j1,j3)*za(j3,j5)*zb(j1,j6)*zb(j3,j4)
     &   -za(j1,j2)*za(j3,j5)*zb(j1,j4)*zb(j2,j6)
     &   -za(j1,j3)*zb(j4,j6)*zab2(j5,j2,j6,j1)
     &   -za(j1,j5)*zb(j1,j4)*zab2(j3,j4,j5,j6)
     & -2*za(j1,j3)*zb(j1,j4)*zab2(j5,j1,j2,j6))
     & /(za(j2,j7)*zab2(j7,j5,j6,j1))

      res=res/(2._dp*s(j3,j4)*s(j5,j6))

      return
      end

