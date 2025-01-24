!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function wpwmamp(j1,j2,j3,j4,j5,j6,j7,j8)
      implicit none
c----Author R.K.Ellis June 2006
c----Amplitude for
c     u(-p1)+s(-p2)-->W^+ + W^- d(p7)+c(p8)
c                     |     |
c                     |     --> e^-(p5)+nu~(p6)
c                     |
c                     |---> num(p3)+mu^+(p4)
c    with a factor of 8*gsq*gw**4/4 and color factor V/4 or -V/4/xn removed

      include 'types.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'

      real(dp):: s34,s56,t134,t347,t256,t568,t1347
      complex(dp):: diag1a,diag2a,diag3a,diag4a
      complex(dp):: propwp1,propwp2,zba2,wpwmamp
      integer j1,j2,j3,j4,j5,j6,j7,j8
      zba2(j1,j2,j3,j4)=zb(j1,j2)*za(j2,j4)+zb(j1,j3)*za(j3,j4)

      s34=s(j3,j4)
      s56=s(j5,j6)
      t134=s(j1,j3)+s(j3,j4)+s(j4,j1)
      t347=s(j3,j4)+s(j4,j7)+s(j7,j3)
      t256=s(j2,j5)+s(j5,j6)+s(j6,j2)
      t568=s(j5,j6)+s(j6,j8)+s(j8,j5)

      t1347=s(j1,j3)+s(j1,j4)+s(j1,j7)
     &     +s(j3,j4)+s(j3,j7)+s(j4,j7)

c      write(6,*) 's34',s34
c      write(6,*) 's56',s56

      diag1a=-za(j7,j3)*zb(j6,j2)*zba2(j4,j3,j7,j8)
     & *(za(j2,j5)*zb(j2,j1)+za(j6,j5)*zb(j6,j1))
     & /(s34*s56*t347*t256*t1347)

      diag2a=+za(j7,j3)*za(j8,j5)*zb(j2,j1)
     & *(zb(j6,j8)*zba2(j4,j3,j7,j8)+zb(j6,j5)*zba2(j4,j3,j7,j5))
     & /(s34*s56*t347*t568*t1347)

      diag3a=+za(j7,j8)*zb(j4,j1)*zb(j6,j2)
     & *(za(j2,j5)*zba2(j2,j1,j4,j3)+za(j6,j5)*zba2(j6,j1,j4,j3))
     & /(s34*s56*t134*t256*t1347)

      diag4a=-za(j8,j5)*zb(j4,j1)*zba2(j2,j1,j4,j3)
     & *(za(j7,j5)*zb(j6,j5)+za(j7,j8)*zb(j6,j8))
     & /(s34*s56*t134*t568*t1347)

      propwp1=cmplx(s34,kind=dp)/cmplx(s34-wmass**2,wmass*wwidth,dp)
      propwp2=cmplx(s56,kind=dp)/cmplx(s56-wmass**2,wmass*wwidth,dp)
      wpwmamp=(diag1a+diag2a+diag3a+diag4a)*propwp1*propwp2

c      ampsq(1)=cdabs(diag1a*propwp1*propwp2)**2*0.019644281d0
c      ampsq(2)=cdabs(diag2a*propwp1*propwp2)**2*0.019644281d0
c      ampsq(3)=cdabs(diag3a*propwp1*propwp2)**2*0.019644281d0
c      ampsq(4)=cdabs(diag4a*propwp1*propwp2)**2*0.019644281d0
      return
      end


