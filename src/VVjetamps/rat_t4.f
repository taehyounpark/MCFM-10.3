!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine rat_t4(j1,j2,j3,j4,j5,j6,j7,za,zb,res)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4,j5,j6,j7
      complex(dp):: res,zab2,den1,den2,
     & zab3271,zab4271,zab5271,zab1273,zab1274,zab1276
      real(dp):: s12p17

c--- statement function
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c--- end statement function

      s12p17=s(j1,j2)+s(j1,j7)
      den1=zab2(j1,j2,j7,j4)*za(j4,j1)+zab2(j1,j2,j7,j3)*za(j3,j1)
      den2=zb(j1,j4)*zab2(j4,j2,j7,j1)+zb(j1,j3)*zab2(j3,j2,j7,j1)

      zab3271=zab2(j3,j2,j7,j1)
      zab4271=zab2(j4,j2,j7,j1)
      zab5271=zab2(j5,j2,j7,j1)
      zab1273=zab2(j1,j2,j7,j3)
      zab1274=zab2(j1,j2,j7,j4)
      zab1276=zab2(j1,j2,j7,j6)

      res=
     & (-((zb(j1,j4)*zab4271+zb(j1,j3)*zab3271)
     & *zb(j7,j2)*zb(j7,j4)
     & *(zab1274*zb(j1,j6)*za(j3,j5)*za(j4,j1)*s12p17
     &  +zab1273*zb(j1,j6)*za(j3,j1)*za(j3,j5)*s12p17
     &  +zab1276*zb(j1,j4)*za(j3,j1)*zab4271*za(j5,j1)
     &  +zab1274*zb(j1,j6)*zab3271*za(j4,j1)*za(j5,j1)
     &  +zab1273*zb(j1,j6)*za(j3,j1)*zab3271*za(j5,j1)
     &  +zab1276*zb(j1,j3)*za(j3,j1)*zab3271*za(j5,j1)))
     & /(den1*den2**2*s(j2,j7)))
      res=res
     & +(zb(j1,j6)*zb(j4,j2)*zb(j7,j1)
     & *(za(j3,j5)*s12p17+zab3271*za(j5,j1)))
     & /(den2*zb(j2,j1)*za(j7,j2))

      res=res
     & +(zb(j1,j6)*zb(j7,j1)
     & *(zb(j1,j4)*za(j1,j7)*za(j3,j1)*zab5271*zb(j7,j2)
     & +zb(j1,j2)*zab1274*zab3271*za(j5,j1)))
     & /(den2*zb(j2,j1)*za(j7,j2)*s12p17)

      res=res
     & +(zb(j1,j6)*zab3271
     & *(zab1274*za(j4,j1)+zab1273*za(j3,j1))
     & *(zab1274*zb(j1,j4)*za(j4,j1)*zab5271
     & +zab1273*zb(j1,j4)*za(j3,j1)*zab5271
     & -zab1274*zb(j1,j4)*zab4271*za(j5,j1)
     & -zab1274*zb(j1,j3)*zab3271*za(j5,j1))
     & *zb(j7,j1)*zb(j7,j2))
     & /(den1*den2**2*s(j2,j7)*s12p17)

      res=res
     & +(za(j1,j2)*zb(j1,j6)*za(j3,j1)
     & *(zb(j1,j4)*zab4271+zb(j1,j3)*zab3271)
     & *(zab1274*zb(j1,j4)*za(j4,j1)*zab5271
     & +zab1273*zb(j1,j4)*za(j3,j1)*zab5271
     & +zab1274*zb(j1,j4)*zab4271*za(j5,j1)
     & +zab1274*zb(j1,j3)*zab3271*za(j5,j1))
     & *zb(j7,j2))
     & /(den1*den2**2*za(j7,j2)*s12p17)

      res=res
     & -(zab1274*zb(j1,j6)*za(j2,j1)
     & *za(j3,j1)*za(j5,j1)*zb(j7,j2))
     & /(den1*za(j2,j7)*s12p17)

      res=res
     & -(den1*zb(j1,j2)*zb(j1,j4)*zb(j1,j6)*zab3271*zab5271
     & *zb(j7,j1))
     & /(den2**2*zb(j2,j1)*za(j7,j2)*s12p17)

      res=res
     & +(s(j1,j2)*zb(j1,j4)*zb(j1,j6)*za(j2,j1)*zb(j2,j7)
     & *zab3271*zab5271)
     & /(den2*zb(j2,j1)*za(j2,j7)*zb(j7,j1)*za(j7,j2)*s12p17)

      res=res
     & +s(j1,j2)*zb(j1,j4)*zb(j1,j6)*za(j1,j7)*zb(j2,j7)*zab3271*zab5271
     & /(den2*zb(j2,j1)**2*za(j7,j2)**2*s12p17)
     & -(zab1276*za(j2,j1)*za(j3,j1)*zb(j4,j2)*za(j5,j1))
     & /(den1*za(j1,j7)*za(j2,j7))
      res=res/(2d0*s(j4,j3)*s(j6,j5))
      return
      end

