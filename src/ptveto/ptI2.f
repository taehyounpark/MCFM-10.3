!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine pTI2(z,R,Lperp,LQ,I2)
c     Elaborated using the form results of 2207.07037v1
      implicit none
      include 'types.f'
      include 'distributions.f'
c     dmin=-1,dmax=2,delt=-1,plus=0,lpls=1,rglr=2
      real(dp):: z,R,Lperp,LQ,I2(0:6,dmin:dmax),
     & DeltaI2(0:6,dmin:dmax),I2P(0:6,dmin:dmax)
c      I2, index one transition flavor
c      I2, index2,distribution type
c      I2, index3,power of Lb

      call I2perp(z,Lperp,LQ,I2P)
      call DeltaI2R(z,LQ,R,DeltaI2)

      I2(:,:)=I2P(:,:)+DeltaI2(:,:)

      return
      end
