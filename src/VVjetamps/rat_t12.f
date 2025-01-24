!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine rat_t12(j1,j2,j3,j4,j5,j6,j7,za,zb,res)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4,j5,j6,j7
      complex(dp):: res,zab2,zaba,zbab,
     & zab7124,zab7126,zab3127,zab5127
      real(dp):: s17p27

c--- statement function
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c--- end statement function
      s17p27=s(j1,j7)+s(j2,j7)
      zaba=zab2(j7,j1,j2,j4)*za(j4,j7)+zab2(j7,j1,j2,j3)*za(j3,j7)
      zbab=zb(j7,j4)*zab2(j4,j1,j2,j7)+zb(j7,j3)*zab2(j3,j1,j2,j7)

      zab7126=zab2(j7,j1,j2,j6)
      zab7124=zab2(j7,j1,j2,j4)
      zab3127=zab2(j3,j1,j2,j7)
      zab5127=zab2(j5,j1,j2,j7)

       res= +za(j1,j3)*zb(j2,j4)*za(j5,j7)*zab7126*s17p27
     & /(zaba*zb(j1,j2)*za(j2,j7)*za(j7,j1))

     & +za(j3,j7)*za(j5,j1)*zb(j6,j4)*s17p27/(zaba*za(j2,j7))

     & +za(j1,j2)*zb(j2,j7)*zab3127*zb(j7,j4)*zb(j7,j6)*zab5127
     & /(zbab*zb(j1,j7)*za(j2,j1)*s17p27)

     & +za(j3,j7)*(za(j1,j2)*zb(j2,j7)*za(j5,j7)*zab7124*zb(j7,j6)
     & +za(j1,j7)*zab5127*zab7126*zb(j7,j4))
     & /(zaba*za(j2,j7)*s17p27)

     & +s(j1,j7)*za(j1,j2)*za(j3,j7)*za(j5,j7)*zab7124*zab7126*zb(j2,j7)
     & /(zaba*zb(j1,j2)*za(j2,j7)**2*za(j7,j1)*s17p27)

     & -za(j1,j7)*zab3127*za(j5,j7)*zab7124*zab7126*zb(j2,j7)
     &  /(zaba*zb(j1,j2)*za(j2,j7)*za(j7,j1)*s17p27)

     & -(zbab*za(j1,j7)*za(j3,j7)*za(j5,j7)*zab7124*zab7126)
     &  /(zaba**2*za(j2,j7)*s17p27)

      res=res/(s(j3,j4)*s(j5,j6))
      return
      end

