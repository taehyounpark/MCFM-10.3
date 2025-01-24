!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine rat_t2(j1,j2,j3,j4,j5,j6,j7,za,zb,res)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4,j5,j6,j7
      complex(dp):: res,zab2,zab2156
      real(dp)::s156
c--- statement function
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      s156=s(j1,j5)+s(j1,j6)+s(j5,j6)
      zab2156=zab2(j2,j1,j5,j6)
      res=
     & -2*za(j1,j5)*zb(j4,j7)*zab2(j3,j1,j5,j6)/(za(j2,j7)*s156)
     & -za(j1,j5)*za(j2,j3)*zb(j2,j7)*zb(j4,j7)*zab2156
     & /(za(j2,j7)*zab2(j2,j3,j4,j7)*s156)

     & +za(j2,j3)*za(j2,j5)*zb(j2,j7)*zb(j4,j7)*zb(j6,j7)
     & /(zab2(j2,j3,j4,j7)*za(j2,j7)*zb(j1,j7))

     & -za(j1,j5)*za(j3,j7)*zb(j2,j4)*zb(j2,j7)*zab2156
     & /(za(j2,j7)*zab2(j7,j3,j4,j2)*s156)

     & +za(j3,j7)*za(j5,j7)*za(j1,j2)*zb(j2,j4)*zb(j2,j6)*zb(j2,j7)
     & /(za(j1,j7)*za(j2,j7)*zb(j1,j2)*zab2(j7,j3,j4,j2))
      res=-res/(2d0*s(j3,j4)*s(j5,j6))
      return
      end

