!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine diag1pm2(k1,k2,k3,k4,k5,k6,k7,za,zb,gl2,gr2,glr)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: k1,k2,k3,k4,k5,k6,k7
      complex(dp):: gl2,gr2,glr,zab2,zaba22,zbab22,
     & ang7x56x12x7,sqr7x12x56x7
      real(dp):: sd12x56,s17p27,s127,s567,s12,s34,s56,s3
c--- statement functions
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
      zaba22(k1,k2,k3,k4,k5,k6)=
     & zab2(k1,k2,k3,k4)*za(k4,k6)+zab2(k1,k2,k3,k5)*za(k5,k6)
      zbab22(k1,k2,k3,k4,k5,k6)=
     & zb(k1,k2)*zab2(k2,k4,k5,k6)+zb(k1,k3)*zab2(k3,k4,k5,k6)
      s3(k1,k2,k3)=s(k1,k2)+s(k2,k3)+s(k3,k1)
c--- end statement functions
      s17p27=s(k1,k7)+s(k2,k7)
      s127=s3(k1,k2,k7)
      s567=s3(k5,k6,k7)
      s12=s(k1,k2)
      s34=s(k3,k4)
      s56=s(k5,k6)
      ang7x56x12x7=zaba22(k7,k5,k6,k1,k2,k7)
      sqr7x12x56x7=zbab22(k7,k1,k2,k5,k6,k7)
      sd12x56=s127*s567-s12*s56
      gl2=1.D0/2.D0*za(k1,k2)*za(k1,k3)*za(k5,k6)*zb(k1,k2)*zb(k2,k7)*
     & zb(k4,k6)*zb(k6,k7)*s17p27**(-1) + 1.D0/2.D0*za(k1,k2)*za(k1,k3)
     & *za(k5,k7)*zb(k1,k2)*zb(k2,k7)*zb(k4,k7)*zb(k6,k7)*s17p27**(-2)*
     & s12 + 1.D0/2.D0*za(k1,k2)*za(k1,k3)*za(k5,k7)*zb(k1,k2)*zb(k2,k7
     & )*zb(k4,k7)*zb(k6,k7)*s17p27**(-1) - 1.D0/4.D0*za(k1,k2)*za(k1,
     & k3)*za(k5,k7)*zb(k1,k2)*zb(k2,k7)*zb(k6,k7)*zab2(k7,k1,k2,k4)*
     & sd12x56*s17p27**(-2)*ang7x56x12x7**(-1) + 1.D0/4.D0*za(k1,k2)*
     & za(k1,k3)*zb(k1,k2)*zb(k2,k7)*zb(k4,k7)*zb(k6,k7)*zab2(k5,k1,k2,
     & k7)*sd12x56*s17p27**(-2)*sqr7x12x56x7**(-1) - 1.D0/2.D0*za(k1,k2
     & )*za(k1,k7)*za(k3,k5)*za(k5,k7)*zb(k1,k2)*zb(k2,k4)*zb(k5,k6)*
     & s17p27**(-1)*ang7x56x12x7**(-1)*sqr7x12x56x7 + 1.D0/4.D0*za(k1,
     & k2)*za(k1,k7)*za(k3,k5)*za(k5,k7)*zb(k1,k2)*zb(k2,k7)*zb(k4,k7)*
     & zb(k5,k6)*sd12x56*s17p27**(-2)*ang7x56x12x7**(-1) - 1.D0/4.D0*
     & za(k1,k2)*za(k1,k7)*za(k3,k7)*za(k5,k6)*zb(k1,k2)*zb(k2,k7)*zb(
     & k4,k6)*zb(k6,k7)*sd12x56*s17p27**(-2)*ang7x56x12x7**(-1) + 1.D0/
     & 2.D0*za(k1,k2)*za(k1,k7)*za(k3,k7)*za(k5,k7)*zb(k1,k2)*zb(k2,k4)
     & *zb(k6,k7)*s17p27**(-2)*ang7x56x12x7**(-1)*sqr7x12x56x7*s12
      gl2 = gl2 + 1.D0/2.D0*za(k1,k2)*za(k1,k7)*za(k3,k7)*za(k5,k7)*zb(
     & k1,k2)*zb(k2,k4)*zb(k6,k7)*s17p27**(-1)*ang7x56x12x7**(-1)*
     & sqr7x12x56x7 - 3.D0/4.D0*za(k1,k2)*za(k1,k7)*za(k3,k7)*za(k5,k7)
     & *zb(k1,k2)*zb(k2,k4)*zab2(k7,k1,k2,k6)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-2)*sqr7x12x56x7 - 1.D0/2.D0*za(k1,k2)*za(k1,k7)*
     & za(k3,k7)*za(k5,k7)*zb(k1,k2)*zb(k2,k7)*zb(k4,k7)*zb(k6,k7)*
     & sd12x56*s17p27**(-3)*ang7x56x12x7**(-1)*s12 - 1.D0/2.D0*za(k1,k2
     & )*za(k1,k7)*za(k3,k7)*za(k5,k7)*zb(k1,k2)*zb(k2,k7)*zb(k4,k7)*
     & zb(k6,k7)*sd12x56*s17p27**(-2)*ang7x56x12x7**(-1) + 1.D0/4.D0*
     & za(k1,k2)*za(k1,k7)*za(k3,k7)*za(k5,k7)*zb(k1,k2)*zb(k2,k7)*zb(
     & k4,k7)*zab2(k7,k1,k2,k6)*sd12x56**2*s17p27**(-3)*
     & ang7x56x12x7**(-2) + 1.D0/4.D0*za(k1,k2)*za(k1,k7)*za(k3,k7)*za(
     & k5,k7)*zb(k1,k2)*zb(k2,k7)*zb(k6,k7)*zab2(k7,k1,k2,k4)*
     & sd12x56**2*s17p27**(-3)*ang7x56x12x7**(-2) - 1.D0/4.D0*za(k1,k2)
     & *za(k1,k7)*za(k5,k7)*zb(k1,k2)*zb(k2,k4)*zb(k6,k7)*zab2(k3,k1,k2
     & ,k7)*sd12x56*s17p27**(-2)*ang7x56x12x7**(-1)
      gl2 = gl2 - 1.D0/4.D0*za(k1,k2)*za(k3,k5)*za(k5,k7)*zb(k2,k4)*zb(
     & k2,k7)*zb(k5,k6)*zab2(k7,k1,k2,k7)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-1) + 1.D0/2.D0*za(k1,k2)*za(k3,k5)*za(k5,k7)*zb(
     & k2,k7)**2*zb(k4,k7)*zb(k5,k6)*zab2(k7,k1,k2,k7)*s17p27**(-3)*s12
     &  + 1.D0/2.D0*za(k1,k2)*za(k3,k5)*zb(k2,k4)*zb(k2,k7)*zb(k5,k6)*
     & zab2(k5,k1,k2,k7)*s17p27**(-1) + 1.D0/4.D0*za(k1,k2)*za(k3,k5)*
     & zb(k2,k7)**2*zb(k4,k7)*zb(k5,k6)*zab2(k5,k1,k2,k7)*sd12x56*
     & s17p27**(-2)*sqr7x12x56x7**(-1) - 1.D0/2.D0*za(k1,k2)*za(k3,k7)*
     & za(k5,k6)*zb(k2,k7)**2*zb(k4,k6)*zb(k6,k7)*s17p27**(-2)*s12 + 1.D
     & 0/4.D0*za(k1,k2)*za(k3,k7)*za(k5,k7)*zb(k2,k4)*zb(k2,k7)*zb(k6,
     & k7)*zab2(k7,k1,k2,k7)*sd12x56*s17p27**(-3)*ang7x56x12x7**(-1)*
     & s12 + 1.D0/4.D0*za(k1,k2)*za(k3,k7)*za(k5,k7)*zb(k2,k4)*zb(k2,k7
     & )*zb(k6,k7)*zab2(k7,k1,k2,k7)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-1) - 1.D0/2.D0*za(k1,k2)*za(k3,k7)*za(k5,k7)*zb(
     & k2,k4)*zb(k2,k7)*zab2(k7,k1,k2,k6)*s17p27**(-2)*
     & ang7x56x12x7**(-1)*sqr7x12x56x7*s12
      gl2 = gl2 - 1.D0/4.D0*za(k1,k2)*za(k3,k7)*za(k5,k7)*zb(k2,k4)*zb(
     & k2,k7)*zab2(k7,k1,k2,k6)*sd12x56**2*s17p27**(-2)*
     & ang7x56x12x7**(-2) - 1.D0/2.D0*za(k1,k2)*za(k3,k7)*za(k5,k7)*zb(
     & k2,k7)**2*zb(k4,k7)*zb(k6,k7)*s17p27**(-3)*s12**2 - 1.D0/2.D0*
     & za(k1,k2)*za(k3,k7)*za(k5,k7)*zb(k2,k7)**2*zb(k4,k7)*zb(k6,k7)*
     & s17p27**(-2)*s12 - 1.D0/2.D0*za(k1,k2)*za(k3,k7)*za(k5,k7)*zb(k2
     & ,k7)**2*zb(k4,k7)*zb(k6,k7)*zab2(k7,k1,k2,k7)*s17p27**(-4)*
     & s12**2 - 1.D0/2.D0*za(k1,k2)*za(k3,k7)*za(k5,k7)*zb(k2,k7)**2*
     & zb(k4,k7)*zb(k6,k7)*zab2(k7,k1,k2,k7)*s17p27**(-3)*s12 + 1.D0/4.D
     & 0*za(k1,k2)*za(k3,k7)*za(k5,k7)*zb(k2,k7)**2*zb(k4,k7)*zab2(k7,
     & k1,k2,k6)*sd12x56*s17p27**(-3)*ang7x56x12x7**(-1)*s12 + 1.D0/4.D0
     & *za(k1,k2)*za(k3,k7)*za(k5,k7)*zb(k2,k7)**2*zb(k6,k7)*zab2(k7,k1
     & ,k2,k4)*sd12x56*s17p27**(-3)*ang7x56x12x7**(-1)*s12 - 1.D0/2.D0*
     & za(k1,k2)*za(k3,k7)*zb(k2,k4)*zb(k2,k7)*zb(k6,k7)*zab2(k5,k1,k2,
     & k7)*s17p27**(-2)*s12 - 1.D0/2.D0*za(k1,k2)*za(k3,k7)*zb(k2,k4)*
     & zb(k2,k7)*zb(k6,k7)*zab2(k5,k1,k2,k7)*s17p27**(-1)
      gl2 = gl2 + 1.D0/4.D0*za(k1,k2)*za(k3,k7)*zb(k2,k4)*zb(k2,k7)*
     & zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k6)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-1) - 1.D0/2.D0*za(k1,k2)*za(k3,k7)*zb(k2,k7)**2*
     & zb(k4,k7)*zb(k6,k7)*zab2(k5,k1,k2,k7)*sd12x56*s17p27**(-3)*
     & sqr7x12x56x7**(-1)*s12 - 1.D0/4.D0*za(k1,k2)*za(k3,k7)*zb(k2,k7)
     & **2*zb(k4,k7)*zb(k6,k7)*zab2(k5,k1,k2,k7)*sd12x56*s17p27**(-2)*
     & sqr7x12x56x7**(-1) - 1.D0/4.D0*za(k1,k2)*za(k5,k6)*zb(k2,k7)**2*
     & zb(k4,k6)*zb(k6,k7)*zab2(k3,k1,k2,k7)*sd12x56*s17p27**(-2)*
     & sqr7x12x56x7**(-1) - 1.D0/2.D0*za(k1,k2)*za(k5,k7)*zb(k2,k4)*zb(
     & k2,k7)*zb(k6,k7)*zab2(k3,k1,k2,k7)*zab2(k7,k1,k2,k7)*
     & s17p27**(-3)*s12 - 1.D0/4.D0*za(k1,k2)*za(k5,k7)*zb(k2,k7)**2*
     & zb(k4,k7)*zb(k6,k7)*zab2(k3,k1,k2,k7)*sd12x56*s17p27**(-3)*
     & sqr7x12x56x7**(-1)*s12 - 1.D0/4.D0*za(k1,k2)*za(k5,k7)*zb(k2,k7)
     & **2*zb(k4,k7)*zb(k6,k7)*zab2(k3,k1,k2,k7)*sd12x56*s17p27**(-2)*
     & sqr7x12x56x7**(-1) - 1.D0/4.D0*za(k1,k2)*za(k5,k7)*zb(k2,k7)**2*
     & zb(k4,k7)*zb(k6,k7)*zab2(k3,k1,k2,k7)*zab2(k7,k1,k2,k7)*sd12x56*
     & s17p27**(-4)*sqr7x12x56x7**(-1)*s12
      gl2 = gl2 - 1.D0/4.D0*za(k1,k2)*zb(k2,k4)*zb(k2,k7)*zb(k6,k7)*
     & zab2(k3,k1,k2,k7)*zab2(k5,k1,k2,k7)*sd12x56*s17p27**(-2)*
     & sqr7x12x56x7**(-1) - 1.D0/2.D0*za(k1,k2)*zb(k2,k7)**2*zb(k4,k7)*
     & zb(k6,k7)*zab2(k3,k1,k2,k7)*zab2(k5,k1,k2,k7)*sd12x56**2*
     & s17p27**(-3)*sqr7x12x56x7**(-2) - 1.D0/4.D0*za(k1,k3)*za(k1,k7)*
     & za(k5,k6)*zb(k1,k2)*zb(k4,k6)*zb(k6,k7)*sd12x56*s17p27**(-1)*
     & ang7x56x12x7**(-1) + 1.D0/2.D0*za(k1,k3)*za(k1,k7)*za(k5,k6)*zb(
     & k1,k2)*zb(k4,k6)*zab2(k7,k1,k2,k6)*s17p27**(-1)*
     & ang7x56x12x7**(-1)*sqr7x12x56x7 - 1.D0/4.D0*za(k1,k3)*za(k1,k7)*
     & za(k5,k7)*zb(k1,k2)*zb(k4,k7)*zb(k6,k7)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-1)*s12 - 1.D0/4.D0*za(k1,k3)*za(k1,k7)*za(k5,k7)
     & *zb(k1,k2)*zb(k4,k7)*zb(k6,k7)*sd12x56*s17p27**(-1)*
     & ang7x56x12x7**(-1) + 1.D0/2.D0*za(k1,k3)*za(k1,k7)*za(k5,k7)*zb(
     & k1,k2)*zb(k4,k7)*zab2(k7,k1,k2,k6)*s17p27**(-2)*
     & ang7x56x12x7**(-1)*sqr7x12x56x7*s12 + 1.D0/2.D0*za(k1,k3)*za(k1,
     & k7)*za(k5,k7)*zb(k1,k2)*zb(k4,k7)*zab2(k7,k1,k2,k6)*s17p27**(-1)
     & *ang7x56x12x7**(-1)*sqr7x12x56x7
      gl2 = gl2 + 1.D0/2.D0*za(k1,k3)*za(k1,k7)*za(k5,k7)*zb(k1,k2)*zb(
     & k6,k7)*zab2(k7,k1,k2,k4)*s17p27**(-2)*ang7x56x12x7**(-1)*
     & sqr7x12x56x7*s12 + 1.D0/2.D0*za(k1,k3)*za(k1,k7)*za(k5,k7)*zb(k1
     & ,k2)*zb(k6,k7)*zab2(k7,k1,k2,k4)*s17p27**(-1)*ang7x56x12x7**(-1)
     & *sqr7x12x56x7 + 1.D0/4.D0*za(k1,k3)*za(k1,k7)*za(k5,k7)*zb(k1,k2
     & )*zb(k6,k7)*zab2(k7,k1,k2,k4)*sd12x56**2*s17p27**(-2)*
     & ang7x56x12x7**(-2) - 3.D0/4.D0*za(k1,k3)*za(k1,k7)*za(k5,k7)*zb(
     & k1,k2)*zab2(k7,k1,k2,k4)*zab2(k7,k1,k2,k6)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-2)*sqr7x12x56x7 + 1.D0/2.D0*za(k1,k3)*za(k1,k7)*
     & zb(k1,k2)*zb(k4,k7)*zb(k6,k7)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k7
     & )*s17p27**(-3)*s12 + 1.D0/2.D0*za(k1,k3)*za(k1,k7)*zb(k1,k2)*zb(
     & k4,k7)*zb(k6,k7)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k7)*
     & s17p27**(-2) - 1.D0/4.D0*za(k1,k3)*za(k1,k7)*zb(k1,k2)*zb(k4,k7)
     & *zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k6)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-1) + 1.D0/2.D0*za(k1,k3)*za(k5,k6)*zb(k2,k7)*zb(
     & k4,k6)*zb(k6,k7)
      gl2 = gl2 - 1.D0/4.D0*za(k1,k3)*za(k5,k6)*zb(k2,k7)*zb(k4,k6)*
     & zab2(k7,k1,k2,k6)*sd12x56*s17p27**(-1)*ang7x56x12x7**(-1) + 1.D0/
     & 2.D0*za(k1,k3)*za(k5,k7)*zb(k2,k7)*zb(k4,k7)*zb(k6,k7) + 1.D0/2.D
     & 0*za(k1,k3)*za(k5,k7)*zb(k2,k7)*zb(k4,k7)*zb(k6,k7)*s17p27**(-1)
     & *s12 - 1.D0/4.D0*za(k1,k3)*za(k5,k7)*zb(k2,k7)*zb(k4,k7)*zab2(k7
     & ,k1,k2,k6)*sd12x56*s17p27**(-2)*ang7x56x12x7**(-1)*s12 - 1.D0/4.D
     & 0*za(k1,k3)*za(k5,k7)*zb(k2,k7)*zb(k4,k7)*zab2(k7,k1,k2,k6)*
     & sd12x56*s17p27**(-1)*ang7x56x12x7**(-1) - 1.D0/4.D0*za(k1,k3)*
     & za(k5,k7)*zb(k2,k7)*zb(k6,k7)*zab2(k7,k1,k2,k4)*sd12x56*
     & s17p27**(-1)*ang7x56x12x7**(-1) + 1.D0/4.D0*za(k1,k3)*za(k5,k7)*
     & zb(k2,k7)*zab2(k7,k1,k2,k4)*zab2(k7,k1,k2,k6)*sd12x56**2*
     & s17p27**(-2)*ang7x56x12x7**(-2) + 1.D0/4.D0*za(k1,k3)*zb(k2,k7)*
     & zb(k4,k7)*zb(k6,k7)*zab2(k5,k1,k2,k7)*sd12x56*s17p27**(-1)*
     & sqr7x12x56x7**(-1) - 1.D0/2.D0*za(k1,k7)**2*za(k3,k5)*za(k5,k7)*
     & zb(k1,k2)*zb(k4,k7)*zb(k5,k6)*s17p27**(-2)*ang7x56x12x7**(-1)*
     & sqr7x12x56x7*s12
      gl2 = gl2 + 3.D0/4.D0*za(k1,k7)**2*za(k3,k5)*za(k5,k7)*zb(k1,k2)*
     & zb(k5,k6)*zab2(k7,k1,k2,k4)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-2)*sqr7x12x56x7 + 1.D0/2.D0*za(k1,k7)**2*za(k3,
     & k7)*za(k5,k6)*zb(k1,k2)*zb(k4,k6)*zb(k6,k7)*s17p27**(-2)*
     & ang7x56x12x7**(-1)*sqr7x12x56x7*s12 + 1.D0/2.D0*za(k1,k7)**2*za(
     & k3,k7)*za(k5,k6)*zb(k1,k2)*zb(k4,k6)*zb(k6,k7)*s17p27**(-1)*
     & ang7x56x12x7**(-1)*sqr7x12x56x7 - 3.D0/4.D0*za(k1,k7)**2*za(k3,
     & k7)*za(k5,k6)*zb(k1,k2)*zb(k4,k6)*zab2(k7,k1,k2,k6)*sd12x56*
     & s17p27**(-2)*ang7x56x12x7**(-2)*sqr7x12x56x7 + za(k1,k7)**2*za(
     & k3,k7)*za(k5,k7)*zb(k1,k2)*zb(k4,k7)*zb(k6,k7)*s17p27**(-3)*
     & ang7x56x12x7**(-1)*sqr7x12x56x7*s12**2 + 3.D0/2.D0*za(k1,k7)**2*
     & za(k3,k7)*za(k5,k7)*zb(k1,k2)*zb(k4,k7)*zb(k6,k7)*s17p27**(-2)*
     & ang7x56x12x7**(-1)*sqr7x12x56x7*s12 + 1.D0/2.D0*za(k1,k7)**2*za(
     & k3,k7)*za(k5,k7)*zb(k1,k2)*zb(k4,k7)*zb(k6,k7)*s17p27**(-1)*
     & ang7x56x12x7**(-1)*sqr7x12x56x7 - 3.D0/2.D0*za(k1,k7)**2*za(k3,
     & k7)*za(k5,k7)*zb(k1,k2)*zb(k4,k7)*zab2(k7,k1,k2,k6)*sd12x56*
     & s17p27**(-3)*ang7x56x12x7**(-2)*sqr7x12x56x7*s12
      gl2 = gl2 - 3.D0/4.D0*za(k1,k7)**2*za(k3,k7)*za(k5,k7)*zb(k1,k2)*
     & zb(k4,k7)*zab2(k7,k1,k2,k6)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-2)*sqr7x12x56x7 - 3.D0/2.D0*za(k1,k7)**2*za(k3,
     & k7)*za(k5,k7)*zb(k1,k2)*zb(k6,k7)*zab2(k7,k1,k2,k4)*sd12x56*
     & s17p27**(-3)*ang7x56x12x7**(-2)*sqr7x12x56x7*s12 - 3.D0/2.D0*za(
     & k1,k7)**2*za(k3,k7)*za(k5,k7)*zb(k1,k2)*zb(k6,k7)*zab2(k7,k1,k2,
     & k4)*sd12x56*s17p27**(-2)*ang7x56x12x7**(-2)*sqr7x12x56x7 + 2*za(
     & k1,k7)**2*za(k3,k7)*za(k5,k7)*zb(k1,k2)*zab2(k7,k1,k2,k4)*zab2(
     & k7,k1,k2,k6)*sd12x56**2*s17p27**(-3)*ang7x56x12x7**(-3)*
     & sqr7x12x56x7 - 1.D0/4.D0*za(k1,k7)**2*za(k3,k7)*zb(k1,k2)*zb(k4,
     & k7)*zb(k6,k7)*zab2(k5,k1,k2,k7)*sd12x56*s17p27**(-3)*
     & ang7x56x12x7**(-1)*s12 - 1.D0/4.D0*za(k1,k7)**2*za(k3,k7)*zb(k1,
     & k2)*zb(k4,k7)*zb(k6,k7)*zab2(k5,k1,k2,k7)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-1) + 1.D0/4.D0*za(k1,k7)**2*za(k3,k7)*zb(k1,k2)*
     & zb(k4,k7)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k6)*sd12x56**2*
     & s17p27**(-3)*ang7x56x12x7**(-2)
      gl2 = gl2 - 1.D0/4.D0*za(k1,k7)**2*za(k5,k7)*zb(k1,k2)*zb(k4,k7)*
     & zb(k6,k7)*zab2(k3,k1,k2,k7)*sd12x56*s17p27**(-3)*
     & ang7x56x12x7**(-1)*s12 + 1.D0/4.D0*za(k1,k7)**2*za(k5,k7)*zb(k1,
     & k2)*zb(k6,k7)*zab2(k3,k1,k2,k7)*zab2(k7,k1,k2,k4)*sd12x56**2*
     & s17p27**(-3)*ang7x56x12x7**(-2) + 1.D0/2.D0*za(k1,k7)*za(k3,k5)*
     & za(k5,k7)*zb(k2,k7)*zb(k5,k6)*zab2(k7,k1,k2,k4)*s17p27**(-2)*
     & ang7x56x12x7**(-1)*sqr7x12x56x7*s12 + 1.D0/4.D0*za(k1,k7)*za(k3,
     & k5)*zb(k2,k4)*zb(k5,k6)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k7)*
     & sd12x56*s17p27**(-2)*ang7x56x12x7**(-1) + 1.D0/2.D0*za(k1,k7)*
     & za(k3,k5)*zb(k2,k7)*zb(k4,k7)*zb(k5,k6)*zab2(k5,k1,k2,k7)*
     & s17p27**(-2)*s12 - 1.D0/4.D0*za(k1,k7)*za(k3,k5)*zb(k2,k7)*zb(k5
     & ,k6)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k4)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-1) - 1.D0/2.D0*za(k1,k7)*za(k3,k7)*za(k5,k6)*zb(
     & k2,k7)*zb(k4,k6)*zab2(k7,k1,k2,k6)*s17p27**(-2)*
     & ang7x56x12x7**(-1)*sqr7x12x56x7*s12 - za(k1,k7)*za(k3,k7)*za(k5,
     & k7)*zb(k2,k7)*zb(k4,k7)*zab2(k7,k1,k2,k6)*s17p27**(-3)*
     & ang7x56x12x7**(-1)*sqr7x12x56x7*s12**2
      gl2 = gl2 - 1.D0/2.D0*za(k1,k7)*za(k3,k7)*za(k5,k7)*zb(k2,k7)*zb(
     & k4,k7)*zab2(k7,k1,k2,k6)*s17p27**(-2)*ang7x56x12x7**(-1)*
     & sqr7x12x56x7*s12 - za(k1,k7)*za(k3,k7)*za(k5,k7)*zb(k2,k7)*zb(k6
     & ,k7)*zab2(k7,k1,k2,k4)*s17p27**(-3)*ang7x56x12x7**(-1)*
     & sqr7x12x56x7*s12**2 - za(k1,k7)*za(k3,k7)*za(k5,k7)*zb(k2,k7)*
     & zb(k6,k7)*zab2(k7,k1,k2,k4)*s17p27**(-2)*ang7x56x12x7**(-1)*
     & sqr7x12x56x7*s12 + 3.D0/2.D0*za(k1,k7)*za(k3,k7)*za(k5,k7)*zb(k2
     & ,k7)*zab2(k7,k1,k2,k4)*zab2(k7,k1,k2,k6)*sd12x56*s17p27**(-3)*
     & ang7x56x12x7**(-2)*sqr7x12x56x7*s12 - 1.D0/4.D0*za(k1,k7)*za(k3,
     & k7)*zb(k2,k4)*zb(k6,k7)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k7)*
     & sd12x56*s17p27**(-3)*ang7x56x12x7**(-1)*s12 - 1.D0/4.D0*za(k1,k7
     & )*za(k3,k7)*zb(k2,k4)*zb(k6,k7)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,
     & k7)*sd12x56*s17p27**(-2)*ang7x56x12x7**(-1) + 1.D0/4.D0*za(k1,k7
     & )*za(k3,k7)*zb(k2,k4)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k6)*
     & sd12x56**2*s17p27**(-2)*ang7x56x12x7**(-2) - 1.D0/2.D0*za(k1,k7)
     & *za(k3,k7)*zb(k2,k7)*zb(k4,k7)*zb(k6,k7)*zab2(k5,k1,k2,k7)*
     & s17p27**(-3)*s12**2
      gl2 = gl2 - 1.D0/2.D0*za(k1,k7)*za(k3,k7)*zb(k2,k7)*zb(k4,k7)*zb(
     & k6,k7)*zab2(k5,k1,k2,k7)*s17p27**(-2)*s12 - 1.D0/2.D0*za(k1,k7)*
     & za(k3,k7)*zb(k2,k7)*zb(k4,k7)*zb(k6,k7)*zab2(k5,k1,k2,k7)*zab2(
     & k7,k1,k2,k7)*s17p27**(-4)*s12**2 - 1.D0/2.D0*za(k1,k7)*za(k3,k7)
     & *zb(k2,k7)*zb(k4,k7)*zb(k6,k7)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,
     & k7)*s17p27**(-3)*s12 + 1.D0/2.D0*za(k1,k7)*za(k3,k7)*zb(k2,k7)*
     & zb(k4,k7)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k6)*sd12x56*
     & s17p27**(-3)*ang7x56x12x7**(-1)*s12 + 1.D0/4.D0*za(k1,k7)*za(k3,
     & k7)*zb(k2,k7)*zb(k6,k7)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k4)*
     & sd12x56*s17p27**(-3)*ang7x56x12x7**(-1)*s12 + 1.D0/4.D0*za(k1,k7
     & )*za(k3,k7)*zb(k2,k7)*zb(k6,k7)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,
     & k4)*sd12x56*s17p27**(-2)*ang7x56x12x7**(-1) - 1.D0/4.D0*za(k1,k7
     & )*za(k3,k7)*zb(k2,k7)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k4)*zab2(
     & k7,k1,k2,k6)*sd12x56**2*s17p27**(-3)*ang7x56x12x7**(-2) - 1.D0/2.
     & D0*za(k1,k7)*za(k5,k6)*zb(k2,k7)*zb(k4,k6)*zb(k6,k7)*zab2(k3,k1,
     & k2,k7)*zab2(k7,k1,k2,k7)*s17p27**(-3)*s12
      gl2 = gl2 - 1.D0/2.D0*za(k1,k7)*za(k5,k6)*zb(k2,k7)*zb(k4,k6)*zb(
     & k6,k7)*zab2(k3,k1,k2,k7)*zab2(k7,k1,k2,k7)*s17p27**(-2) + 1.D0/4.
     & D0*za(k1,k7)*za(k5,k6)*zb(k2,k7)*zb(k4,k6)*zab2(k3,k1,k2,k7)*
     & zab2(k7,k1,k2,k6)*sd12x56*s17p27**(-2)*ang7x56x12x7**(-1) - za(
     & k1,k7)*za(k5,k7)*zb(k2,k7)*zb(k4,k7)*zb(k6,k7)*zab2(k3,k1,k2,k7)
     & *zab2(k7,k1,k2,k7)*s17p27**(-4)*s12**2 - za(k1,k7)*za(k5,k7)*zb(
     & k2,k7)*zb(k4,k7)*zb(k6,k7)*zab2(k3,k1,k2,k7)*zab2(k7,k1,k2,k7)*
     & s17p27**(-3)*s12 - 1.D0/2.D0*za(k1,k7)*za(k5,k7)*zb(k2,k7)*zb(k4
     & ,k7)*zb(k6,k7)*zab2(k3,k1,k2,k7)*zab2(k7,k1,k2,k7)*s17p27**(-2)
     &  + 1.D0/4.D0*za(k1,k7)*za(k5,k7)*zb(k2,k7)*zb(k4,k7)*zab2(k3,k1,
     & k2,k7)*zab2(k7,k1,k2,k6)*sd12x56*s17p27**(-3)*ang7x56x12x7**(-1)
     & *s12 + 1.D0/4.D0*za(k1,k7)*za(k5,k7)*zb(k2,k7)*zb(k4,k7)*zab2(k3
     & ,k1,k2,k7)*zab2(k7,k1,k2,k6)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-1) + 1.D0/2.D0*za(k1,k7)*za(k5,k7)*zb(k2,k7)*zb(
     & k6,k7)*zab2(k3,k1,k2,k7)*zab2(k7,k1,k2,k4)*sd12x56*s17p27**(-3)*
     & ang7x56x12x7**(-1)*s12
      gl2 = gl2 + 1.D0/4.D0*za(k1,k7)*za(k5,k7)*zb(k2,k7)*zb(k6,k7)*
     & zab2(k3,k1,k2,k7)*zab2(k7,k1,k2,k4)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-1) - 1.D0/4.D0*za(k1,k7)*za(k5,k7)*zb(k2,k7)*
     & zab2(k3,k1,k2,k7)*zab2(k7,k1,k2,k4)*zab2(k7,k1,k2,k6)*sd12x56**2
     & *s17p27**(-3)*ang7x56x12x7**(-2) - 1.D0/4.D0*za(k1,k7)*zb(k2,k7)
     & *zb(k4,k7)*zb(k6,k7)*zab2(k3,k1,k2,k7)*zab2(k5,k1,k2,k7)*sd12x56
     & *s17p27**(-3)*sqr7x12x56x7**(-1)*s12 - 1.D0/4.D0*za(k1,k7)*zb(k2
     & ,k7)*zb(k4,k7)*zb(k6,k7)*zab2(k3,k1,k2,k7)*zab2(k5,k1,k2,k7)*
     & zab2(k7,k1,k2,k7)*sd12x56*s17p27**(-4)*sqr7x12x56x7**(-1)*s12 -
     & 1.D0/4.D0*za(k1,k7)*zb(k2,k7)*zb(k4,k7)*zb(k6,k7)*zab2(k3,k1,k2,
     & k7)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k7)*sd12x56*s17p27**(-3)*
     & sqr7x12x56x7**(-1)
      gr2=1.D0/2.D0*za(k1,k2)*za(k1,k3)*za(k5,k6)*zb(k1,k2)*zb(k2,k7)*
     & zb(k4,k6)*zb(k6,k7)*s17p27**(-1) + 1.D0/2.D0*za(k1,k2)*za(k1,k3)
     & *za(k5,k7)*zb(k1,k2)*zb(k2,k7)*zb(k4,k7)*zb(k6,k7)*s17p27**(-2)*
     & s12 + 1.D0/2.D0*za(k1,k2)*za(k1,k3)*za(k5,k7)*zb(k1,k2)*zb(k2,k7
     & )*zb(k4,k7)*zb(k6,k7)*s17p27**(-1) - 1.D0/4.D0*za(k1,k2)*za(k1,
     & k3)*za(k5,k7)*zb(k1,k2)*zb(k2,k7)*zb(k6,k7)*zab2(k7,k1,k2,k4)*
     & sd12x56*s17p27**(-2)*ang7x56x12x7**(-1) + 1.D0/4.D0*za(k1,k2)*
     & za(k1,k3)*zb(k1,k2)*zb(k2,k7)*zb(k4,k7)*zb(k6,k7)*zab2(k5,k1,k2,
     & k7)*sd12x56*s17p27**(-2)*sqr7x12x56x7**(-1) - 1.D0/2.D0*za(k1,k2
     & )*za(k1,k7)*za(k3,k5)*za(k5,k7)*zb(k1,k2)*zb(k2,k4)*zb(k5,k6)*
     & s17p27**(-1)*ang7x56x12x7**(-1)*sqr7x12x56x7 + 1.D0/4.D0*za(k1,
     & k2)*za(k1,k7)*za(k3,k5)*za(k5,k7)*zb(k1,k2)*zb(k2,k7)*zb(k4,k7)*
     & zb(k5,k6)*sd12x56*s17p27**(-2)*ang7x56x12x7**(-1) - 1.D0/4.D0*
     & za(k1,k2)*za(k1,k7)*za(k3,k7)*za(k5,k6)*zb(k1,k2)*zb(k2,k7)*zb(
     & k4,k6)*zb(k6,k7)*sd12x56*s17p27**(-2)*ang7x56x12x7**(-1) + 1.D0/
     & 2.D0*za(k1,k2)*za(k1,k7)*za(k3,k7)*za(k5,k7)*zb(k1,k2)*zb(k2,k4)
     & *zb(k6,k7)*s17p27**(-2)*ang7x56x12x7**(-1)*sqr7x12x56x7*s12
      gr2 = gr2 + 1.D0/2.D0*za(k1,k2)*za(k1,k7)*za(k3,k7)*za(k5,k7)*zb(
     & k1,k2)*zb(k2,k4)*zb(k6,k7)*s17p27**(-1)*ang7x56x12x7**(-1)*
     & sqr7x12x56x7 - 3.D0/4.D0*za(k1,k2)*za(k1,k7)*za(k3,k7)*za(k5,k7)
     & *zb(k1,k2)*zb(k2,k4)*zab2(k7,k1,k2,k6)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-2)*sqr7x12x56x7 - 1.D0/2.D0*za(k1,k2)*za(k1,k7)*
     & za(k3,k7)*za(k5,k7)*zb(k1,k2)*zb(k2,k7)*zb(k4,k7)*zb(k6,k7)*
     & sd12x56*s17p27**(-3)*ang7x56x12x7**(-1)*s12 - 1.D0/2.D0*za(k1,k2
     & )*za(k1,k7)*za(k3,k7)*za(k5,k7)*zb(k1,k2)*zb(k2,k7)*zb(k4,k7)*
     & zb(k6,k7)*sd12x56*s17p27**(-2)*ang7x56x12x7**(-1) + 1.D0/4.D0*
     & za(k1,k2)*za(k1,k7)*za(k3,k7)*za(k5,k7)*zb(k1,k2)*zb(k2,k7)*zb(
     & k4,k7)*zab2(k7,k1,k2,k6)*sd12x56**2*s17p27**(-3)*
     & ang7x56x12x7**(-2) + 1.D0/4.D0*za(k1,k2)*za(k1,k7)*za(k3,k7)*za(
     & k5,k7)*zb(k1,k2)*zb(k2,k7)*zb(k6,k7)*zab2(k7,k1,k2,k4)*
     & sd12x56**2*s17p27**(-3)*ang7x56x12x7**(-2) - 1.D0/4.D0*za(k1,k2)
     & *za(k1,k7)*za(k5,k7)*zb(k1,k2)*zb(k2,k4)*zb(k6,k7)*zab2(k3,k1,k2
     & ,k7)*sd12x56*s17p27**(-2)*ang7x56x12x7**(-1)
      gr2 = gr2 - 1.D0/4.D0*za(k1,k2)*za(k3,k5)*za(k5,k7)*zb(k2,k4)*zb(
     & k2,k7)*zb(k5,k6)*zab2(k7,k1,k2,k7)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-1) + 1.D0/2.D0*za(k1,k2)*za(k3,k5)*za(k5,k7)*zb(
     & k2,k7)**2*zb(k4,k7)*zb(k5,k6)*zab2(k7,k1,k2,k7)*s17p27**(-3)*s12
     &  + 1.D0/2.D0*za(k1,k2)*za(k3,k5)*zb(k2,k4)*zb(k2,k7)*zb(k5,k6)*
     & zab2(k5,k1,k2,k7)*s17p27**(-1) + 1.D0/4.D0*za(k1,k2)*za(k3,k5)*
     & zb(k2,k7)**2*zb(k4,k7)*zb(k5,k6)*zab2(k5,k1,k2,k7)*sd12x56*
     & s17p27**(-2)*sqr7x12x56x7**(-1) - 1.D0/2.D0*za(k1,k2)*za(k3,k7)*
     & za(k5,k6)*zb(k2,k7)**2*zb(k4,k6)*zb(k6,k7)*s17p27**(-2)*s12 + 1.D
     & 0/4.D0*za(k1,k2)*za(k3,k7)*za(k5,k7)*zb(k2,k4)*zb(k2,k7)*zb(k6,
     & k7)*zab2(k7,k1,k2,k7)*sd12x56*s17p27**(-3)*ang7x56x12x7**(-1)*
     & s12 + 1.D0/4.D0*za(k1,k2)*za(k3,k7)*za(k5,k7)*zb(k2,k4)*zb(k2,k7
     & )*zb(k6,k7)*zab2(k7,k1,k2,k7)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-1) - 1.D0/2.D0*za(k1,k2)*za(k3,k7)*za(k5,k7)*zb(
     & k2,k4)*zb(k2,k7)*zab2(k7,k1,k2,k6)*s17p27**(-2)*
     & ang7x56x12x7**(-1)*sqr7x12x56x7*s12
      gr2 = gr2 - 1.D0/4.D0*za(k1,k2)*za(k3,k7)*za(k5,k7)*zb(k2,k4)*zb(
     & k2,k7)*zab2(k7,k1,k2,k6)*sd12x56**2*s17p27**(-2)*
     & ang7x56x12x7**(-2) - 1.D0/2.D0*za(k1,k2)*za(k3,k7)*za(k5,k7)*zb(
     & k2,k7)**2*zb(k4,k7)*zb(k6,k7)*s17p27**(-3)*s12**2 - 1.D0/2.D0*
     & za(k1,k2)*za(k3,k7)*za(k5,k7)*zb(k2,k7)**2*zb(k4,k7)*zb(k6,k7)*
     & s17p27**(-2)*s12 - 1.D0/2.D0*za(k1,k2)*za(k3,k7)*za(k5,k7)*zb(k2
     & ,k7)**2*zb(k4,k7)*zb(k6,k7)*zab2(k7,k1,k2,k7)*s17p27**(-4)*
     & s12**2 - 1.D0/2.D0*za(k1,k2)*za(k3,k7)*za(k5,k7)*zb(k2,k7)**2*
     & zb(k4,k7)*zb(k6,k7)*zab2(k7,k1,k2,k7)*s17p27**(-3)*s12 + 1.D0/4.D
     & 0*za(k1,k2)*za(k3,k7)*za(k5,k7)*zb(k2,k7)**2*zb(k4,k7)*zab2(k7,
     & k1,k2,k6)*sd12x56*s17p27**(-3)*ang7x56x12x7**(-1)*s12 + 1.D0/4.D0
     & *za(k1,k2)*za(k3,k7)*za(k5,k7)*zb(k2,k7)**2*zb(k6,k7)*zab2(k7,k1
     & ,k2,k4)*sd12x56*s17p27**(-3)*ang7x56x12x7**(-1)*s12 - 1.D0/2.D0*
     & za(k1,k2)*za(k3,k7)*zb(k2,k4)*zb(k2,k7)*zb(k6,k7)*zab2(k5,k1,k2,
     & k7)*s17p27**(-2)*s12 - 1.D0/2.D0*za(k1,k2)*za(k3,k7)*zb(k2,k4)*
     & zb(k2,k7)*zb(k6,k7)*zab2(k5,k1,k2,k7)*s17p27**(-1)
      gr2 = gr2 + 1.D0/4.D0*za(k1,k2)*za(k3,k7)*zb(k2,k4)*zb(k2,k7)*
     & zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k6)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-1) - 1.D0/2.D0*za(k1,k2)*za(k3,k7)*zb(k2,k7)**2*
     & zb(k4,k7)*zb(k6,k7)*zab2(k5,k1,k2,k7)*sd12x56*s17p27**(-3)*
     & sqr7x12x56x7**(-1)*s12 - 1.D0/4.D0*za(k1,k2)*za(k3,k7)*zb(k2,k7)
     & **2*zb(k4,k7)*zb(k6,k7)*zab2(k5,k1,k2,k7)*sd12x56*s17p27**(-2)*
     & sqr7x12x56x7**(-1) - 1.D0/4.D0*za(k1,k2)*za(k5,k6)*zb(k2,k7)**2*
     & zb(k4,k6)*zb(k6,k7)*zab2(k3,k1,k2,k7)*sd12x56*s17p27**(-2)*
     & sqr7x12x56x7**(-1) - 1.D0/2.D0*za(k1,k2)*za(k5,k7)*zb(k2,k4)*zb(
     & k2,k7)*zb(k6,k7)*zab2(k3,k1,k2,k7)*zab2(k7,k1,k2,k7)*
     & s17p27**(-3)*s12 - 1.D0/4.D0*za(k1,k2)*za(k5,k7)*zb(k2,k7)**2*
     & zb(k4,k7)*zb(k6,k7)*zab2(k3,k1,k2,k7)*sd12x56*s17p27**(-3)*
     & sqr7x12x56x7**(-1)*s12 - 1.D0/4.D0*za(k1,k2)*za(k5,k7)*zb(k2,k7)
     & **2*zb(k4,k7)*zb(k6,k7)*zab2(k3,k1,k2,k7)*sd12x56*s17p27**(-2)*
     & sqr7x12x56x7**(-1) - 1.D0/4.D0*za(k1,k2)*za(k5,k7)*zb(k2,k7)**2*
     & zb(k4,k7)*zb(k6,k7)*zab2(k3,k1,k2,k7)*zab2(k7,k1,k2,k7)*sd12x56*
     & s17p27**(-4)*sqr7x12x56x7**(-1)*s12
      gr2 = gr2 - 1.D0/4.D0*za(k1,k2)*zb(k2,k4)*zb(k2,k7)*zb(k6,k7)*
     & zab2(k3,k1,k2,k7)*zab2(k5,k1,k2,k7)*sd12x56*s17p27**(-2)*
     & sqr7x12x56x7**(-1) - 1.D0/2.D0*za(k1,k2)*zb(k2,k7)**2*zb(k4,k7)*
     & zb(k6,k7)*zab2(k3,k1,k2,k7)*zab2(k5,k1,k2,k7)*sd12x56**2*
     & s17p27**(-3)*sqr7x12x56x7**(-2) - 1.D0/4.D0*za(k1,k3)*za(k1,k7)*
     & za(k5,k6)*zb(k1,k2)*zb(k4,k6)*zb(k6,k7)*sd12x56*s17p27**(-1)*
     & ang7x56x12x7**(-1) + 1.D0/2.D0*za(k1,k3)*za(k1,k7)*za(k5,k6)*zb(
     & k1,k2)*zb(k4,k6)*zab2(k7,k1,k2,k6)*s17p27**(-1)*
     & ang7x56x12x7**(-1)*sqr7x12x56x7 - 1.D0/4.D0*za(k1,k3)*za(k1,k7)*
     & za(k5,k7)*zb(k1,k2)*zb(k4,k7)*zb(k6,k7)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-1)*s12 - 1.D0/4.D0*za(k1,k3)*za(k1,k7)*za(k5,k7)
     & *zb(k1,k2)*zb(k4,k7)*zb(k6,k7)*sd12x56*s17p27**(-1)*
     & ang7x56x12x7**(-1) + 1.D0/2.D0*za(k1,k3)*za(k1,k7)*za(k5,k7)*zb(
     & k1,k2)*zb(k4,k7)*zab2(k7,k1,k2,k6)*s17p27**(-2)*
     & ang7x56x12x7**(-1)*sqr7x12x56x7*s12 + 1.D0/2.D0*za(k1,k3)*za(k1,
     & k7)*za(k5,k7)*zb(k1,k2)*zb(k4,k7)*zab2(k7,k1,k2,k6)*s17p27**(-1)
     & *ang7x56x12x7**(-1)*sqr7x12x56x7
      gr2 = gr2 + 1.D0/2.D0*za(k1,k3)*za(k1,k7)*za(k5,k7)*zb(k1,k2)*zb(
     & k6,k7)*zab2(k7,k1,k2,k4)*s17p27**(-2)*ang7x56x12x7**(-1)*
     & sqr7x12x56x7*s12 + 1.D0/2.D0*za(k1,k3)*za(k1,k7)*za(k5,k7)*zb(k1
     & ,k2)*zb(k6,k7)*zab2(k7,k1,k2,k4)*s17p27**(-1)*ang7x56x12x7**(-1)
     & *sqr7x12x56x7 + 1.D0/4.D0*za(k1,k3)*za(k1,k7)*za(k5,k7)*zb(k1,k2
     & )*zb(k6,k7)*zab2(k7,k1,k2,k4)*sd12x56**2*s17p27**(-2)*
     & ang7x56x12x7**(-2) - 3.D0/4.D0*za(k1,k3)*za(k1,k7)*za(k5,k7)*zb(
     & k1,k2)*zab2(k7,k1,k2,k4)*zab2(k7,k1,k2,k6)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-2)*sqr7x12x56x7 + 1.D0/2.D0*za(k1,k3)*za(k1,k7)*
     & zb(k1,k2)*zb(k4,k7)*zb(k6,k7)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k7
     & )*s17p27**(-3)*s12 + 1.D0/2.D0*za(k1,k3)*za(k1,k7)*zb(k1,k2)*zb(
     & k4,k7)*zb(k6,k7)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k7)*
     & s17p27**(-2) - 1.D0/4.D0*za(k1,k3)*za(k1,k7)*zb(k1,k2)*zb(k4,k7)
     & *zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k6)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-1) + 1.D0/2.D0*za(k1,k3)*za(k5,k6)*zb(k2,k7)*zb(
     & k4,k6)*zb(k6,k7)
      gr2 = gr2 - 1.D0/4.D0*za(k1,k3)*za(k5,k6)*zb(k2,k7)*zb(k4,k6)*
     & zab2(k7,k1,k2,k6)*sd12x56*s17p27**(-1)*ang7x56x12x7**(-1) + 1.D0/
     & 2.D0*za(k1,k3)*za(k5,k7)*zb(k2,k7)*zb(k4,k7)*zb(k6,k7) + 1.D0/2.D
     & 0*za(k1,k3)*za(k5,k7)*zb(k2,k7)*zb(k4,k7)*zb(k6,k7)*s17p27**(-1)
     & *s12 - 1.D0/4.D0*za(k1,k3)*za(k5,k7)*zb(k2,k7)*zb(k4,k7)*zab2(k7
     & ,k1,k2,k6)*sd12x56*s17p27**(-2)*ang7x56x12x7**(-1)*s12 - 1.D0/4.D
     & 0*za(k1,k3)*za(k5,k7)*zb(k2,k7)*zb(k4,k7)*zab2(k7,k1,k2,k6)*
     & sd12x56*s17p27**(-1)*ang7x56x12x7**(-1) - 1.D0/4.D0*za(k1,k3)*
     & za(k5,k7)*zb(k2,k7)*zb(k6,k7)*zab2(k7,k1,k2,k4)*sd12x56*
     & s17p27**(-1)*ang7x56x12x7**(-1) + 1.D0/4.D0*za(k1,k3)*za(k5,k7)*
     & zb(k2,k7)*zab2(k7,k1,k2,k4)*zab2(k7,k1,k2,k6)*sd12x56**2*
     & s17p27**(-2)*ang7x56x12x7**(-2) + 1.D0/4.D0*za(k1,k3)*zb(k2,k7)*
     & zb(k4,k7)*zb(k6,k7)*zab2(k5,k1,k2,k7)*sd12x56*s17p27**(-1)*
     & sqr7x12x56x7**(-1) - 1.D0/2.D0*za(k1,k7)**2*za(k3,k5)*za(k5,k7)*
     & zb(k1,k2)*zb(k4,k7)*zb(k5,k6)*s17p27**(-2)*ang7x56x12x7**(-1)*
     & sqr7x12x56x7*s12
      gr2 = gr2 + 3.D0/4.D0*za(k1,k7)**2*za(k3,k5)*za(k5,k7)*zb(k1,k2)*
     & zb(k5,k6)*zab2(k7,k1,k2,k4)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-2)*sqr7x12x56x7 + 1.D0/2.D0*za(k1,k7)**2*za(k3,
     & k7)*za(k5,k6)*zb(k1,k2)*zb(k4,k6)*zb(k6,k7)*s17p27**(-2)*
     & ang7x56x12x7**(-1)*sqr7x12x56x7*s12 + 1.D0/2.D0*za(k1,k7)**2*za(
     & k3,k7)*za(k5,k6)*zb(k1,k2)*zb(k4,k6)*zb(k6,k7)*s17p27**(-1)*
     & ang7x56x12x7**(-1)*sqr7x12x56x7 - 3.D0/4.D0*za(k1,k7)**2*za(k3,
     & k7)*za(k5,k6)*zb(k1,k2)*zb(k4,k6)*zab2(k7,k1,k2,k6)*sd12x56*
     & s17p27**(-2)*ang7x56x12x7**(-2)*sqr7x12x56x7 + za(k1,k7)**2*za(
     & k3,k7)*za(k5,k7)*zb(k1,k2)*zb(k4,k7)*zb(k6,k7)*s17p27**(-3)*
     & ang7x56x12x7**(-1)*sqr7x12x56x7*s12**2 + 3.D0/2.D0*za(k1,k7)**2*
     & za(k3,k7)*za(k5,k7)*zb(k1,k2)*zb(k4,k7)*zb(k6,k7)*s17p27**(-2)*
     & ang7x56x12x7**(-1)*sqr7x12x56x7*s12 + 1.D0/2.D0*za(k1,k7)**2*za(
     & k3,k7)*za(k5,k7)*zb(k1,k2)*zb(k4,k7)*zb(k6,k7)*s17p27**(-1)*
     & ang7x56x12x7**(-1)*sqr7x12x56x7 - 3.D0/2.D0*za(k1,k7)**2*za(k3,
     & k7)*za(k5,k7)*zb(k1,k2)*zb(k4,k7)*zab2(k7,k1,k2,k6)*sd12x56*
     & s17p27**(-3)*ang7x56x12x7**(-2)*sqr7x12x56x7*s12
      gr2 = gr2 - 3.D0/4.D0*za(k1,k7)**2*za(k3,k7)*za(k5,k7)*zb(k1,k2)*
     & zb(k4,k7)*zab2(k7,k1,k2,k6)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-2)*sqr7x12x56x7 - 3.D0/2.D0*za(k1,k7)**2*za(k3,
     & k7)*za(k5,k7)*zb(k1,k2)*zb(k6,k7)*zab2(k7,k1,k2,k4)*sd12x56*
     & s17p27**(-3)*ang7x56x12x7**(-2)*sqr7x12x56x7*s12 - 3.D0/2.D0*za(
     & k1,k7)**2*za(k3,k7)*za(k5,k7)*zb(k1,k2)*zb(k6,k7)*zab2(k7,k1,k2,
     & k4)*sd12x56*s17p27**(-2)*ang7x56x12x7**(-2)*sqr7x12x56x7 + 2*za(
     & k1,k7)**2*za(k3,k7)*za(k5,k7)*zb(k1,k2)*zab2(k7,k1,k2,k4)*zab2(
     & k7,k1,k2,k6)*sd12x56**2*s17p27**(-3)*ang7x56x12x7**(-3)*
     & sqr7x12x56x7 - 1.D0/4.D0*za(k1,k7)**2*za(k3,k7)*zb(k1,k2)*zb(k4,
     & k7)*zb(k6,k7)*zab2(k5,k1,k2,k7)*sd12x56*s17p27**(-3)*
     & ang7x56x12x7**(-1)*s12 - 1.D0/4.D0*za(k1,k7)**2*za(k3,k7)*zb(k1,
     & k2)*zb(k4,k7)*zb(k6,k7)*zab2(k5,k1,k2,k7)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-1) + 1.D0/4.D0*za(k1,k7)**2*za(k3,k7)*zb(k1,k2)*
     & zb(k4,k7)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k6)*sd12x56**2*
     & s17p27**(-3)*ang7x56x12x7**(-2)
      gr2 = gr2 - 1.D0/4.D0*za(k1,k7)**2*za(k5,k7)*zb(k1,k2)*zb(k4,k7)*
     & zb(k6,k7)*zab2(k3,k1,k2,k7)*sd12x56*s17p27**(-3)*
     & ang7x56x12x7**(-1)*s12 + 1.D0/4.D0*za(k1,k7)**2*za(k5,k7)*zb(k1,
     & k2)*zb(k6,k7)*zab2(k3,k1,k2,k7)*zab2(k7,k1,k2,k4)*sd12x56**2*
     & s17p27**(-3)*ang7x56x12x7**(-2) + 1.D0/2.D0*za(k1,k7)*za(k3,k5)*
     & za(k5,k7)*zb(k2,k7)*zb(k5,k6)*zab2(k7,k1,k2,k4)*s17p27**(-2)*
     & ang7x56x12x7**(-1)*sqr7x12x56x7*s12 + 1.D0/4.D0*za(k1,k7)*za(k3,
     & k5)*zb(k2,k4)*zb(k5,k6)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k7)*
     & sd12x56*s17p27**(-2)*ang7x56x12x7**(-1) + 1.D0/2.D0*za(k1,k7)*
     & za(k3,k5)*zb(k2,k7)*zb(k4,k7)*zb(k5,k6)*zab2(k5,k1,k2,k7)*
     & s17p27**(-2)*s12 - 1.D0/4.D0*za(k1,k7)*za(k3,k5)*zb(k2,k7)*zb(k5
     & ,k6)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k4)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-1) - 1.D0/2.D0*za(k1,k7)*za(k3,k7)*za(k5,k6)*zb(
     & k2,k7)*zb(k4,k6)*zab2(k7,k1,k2,k6)*s17p27**(-2)*
     & ang7x56x12x7**(-1)*sqr7x12x56x7*s12 - za(k1,k7)*za(k3,k7)*za(k5,
     & k7)*zb(k2,k7)*zb(k4,k7)*zab2(k7,k1,k2,k6)*s17p27**(-3)*
     & ang7x56x12x7**(-1)*sqr7x12x56x7*s12**2
      gr2 = gr2 - 1.D0/2.D0*za(k1,k7)*za(k3,k7)*za(k5,k7)*zb(k2,k7)*zb(
     & k4,k7)*zab2(k7,k1,k2,k6)*s17p27**(-2)*ang7x56x12x7**(-1)*
     & sqr7x12x56x7*s12 - za(k1,k7)*za(k3,k7)*za(k5,k7)*zb(k2,k7)*zb(k6
     & ,k7)*zab2(k7,k1,k2,k4)*s17p27**(-3)*ang7x56x12x7**(-1)*
     & sqr7x12x56x7*s12**2 - za(k1,k7)*za(k3,k7)*za(k5,k7)*zb(k2,k7)*
     & zb(k6,k7)*zab2(k7,k1,k2,k4)*s17p27**(-2)*ang7x56x12x7**(-1)*
     & sqr7x12x56x7*s12 + 3.D0/2.D0*za(k1,k7)*za(k3,k7)*za(k5,k7)*zb(k2
     & ,k7)*zab2(k7,k1,k2,k4)*zab2(k7,k1,k2,k6)*sd12x56*s17p27**(-3)*
     & ang7x56x12x7**(-2)*sqr7x12x56x7*s12 - 1.D0/4.D0*za(k1,k7)*za(k3,
     & k7)*zb(k2,k4)*zb(k6,k7)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k7)*
     & sd12x56*s17p27**(-3)*ang7x56x12x7**(-1)*s12 - 1.D0/4.D0*za(k1,k7
     & )*za(k3,k7)*zb(k2,k4)*zb(k6,k7)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,
     & k7)*sd12x56*s17p27**(-2)*ang7x56x12x7**(-1) + 1.D0/4.D0*za(k1,k7
     & )*za(k3,k7)*zb(k2,k4)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k6)*
     & sd12x56**2*s17p27**(-2)*ang7x56x12x7**(-2) - 1.D0/2.D0*za(k1,k7)
     & *za(k3,k7)*zb(k2,k7)*zb(k4,k7)*zb(k6,k7)*zab2(k5,k1,k2,k7)*
     & s17p27**(-3)*s12**2
      gr2 = gr2 - 1.D0/2.D0*za(k1,k7)*za(k3,k7)*zb(k2,k7)*zb(k4,k7)*zb(
     & k6,k7)*zab2(k5,k1,k2,k7)*s17p27**(-2)*s12 - 1.D0/2.D0*za(k1,k7)*
     & za(k3,k7)*zb(k2,k7)*zb(k4,k7)*zb(k6,k7)*zab2(k5,k1,k2,k7)*zab2(
     & k7,k1,k2,k7)*s17p27**(-4)*s12**2 - 1.D0/2.D0*za(k1,k7)*za(k3,k7)
     & *zb(k2,k7)*zb(k4,k7)*zb(k6,k7)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,
     & k7)*s17p27**(-3)*s12 + 1.D0/2.D0*za(k1,k7)*za(k3,k7)*zb(k2,k7)*
     & zb(k4,k7)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k6)*sd12x56*
     & s17p27**(-3)*ang7x56x12x7**(-1)*s12 + 1.D0/4.D0*za(k1,k7)*za(k3,
     & k7)*zb(k2,k7)*zb(k6,k7)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k4)*
     & sd12x56*s17p27**(-3)*ang7x56x12x7**(-1)*s12 + 1.D0/4.D0*za(k1,k7
     & )*za(k3,k7)*zb(k2,k7)*zb(k6,k7)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,
     & k4)*sd12x56*s17p27**(-2)*ang7x56x12x7**(-1) - 1.D0/4.D0*za(k1,k7
     & )*za(k3,k7)*zb(k2,k7)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k4)*zab2(
     & k7,k1,k2,k6)*sd12x56**2*s17p27**(-3)*ang7x56x12x7**(-2) - 1.D0/2.
     & D0*za(k1,k7)*za(k5,k6)*zb(k2,k7)*zb(k4,k6)*zb(k6,k7)*zab2(k3,k1,
     & k2,k7)*zab2(k7,k1,k2,k7)*s17p27**(-3)*s12
      gr2 = gr2 - 1.D0/2.D0*za(k1,k7)*za(k5,k6)*zb(k2,k7)*zb(k4,k6)*zb(
     & k6,k7)*zab2(k3,k1,k2,k7)*zab2(k7,k1,k2,k7)*s17p27**(-2) + 1.D0/4.
     & D0*za(k1,k7)*za(k5,k6)*zb(k2,k7)*zb(k4,k6)*zab2(k3,k1,k2,k7)*
     & zab2(k7,k1,k2,k6)*sd12x56*s17p27**(-2)*ang7x56x12x7**(-1) - za(
     & k1,k7)*za(k5,k7)*zb(k2,k7)*zb(k4,k7)*zb(k6,k7)*zab2(k3,k1,k2,k7)
     & *zab2(k7,k1,k2,k7)*s17p27**(-4)*s12**2 - za(k1,k7)*za(k5,k7)*zb(
     & k2,k7)*zb(k4,k7)*zb(k6,k7)*zab2(k3,k1,k2,k7)*zab2(k7,k1,k2,k7)*
     & s17p27**(-3)*s12 - 1.D0/2.D0*za(k1,k7)*za(k5,k7)*zb(k2,k7)*zb(k4
     & ,k7)*zb(k6,k7)*zab2(k3,k1,k2,k7)*zab2(k7,k1,k2,k7)*s17p27**(-2)
     &  + 1.D0/4.D0*za(k1,k7)*za(k5,k7)*zb(k2,k7)*zb(k4,k7)*zab2(k3,k1,
     & k2,k7)*zab2(k7,k1,k2,k6)*sd12x56*s17p27**(-3)*ang7x56x12x7**(-1)
     & *s12 + 1.D0/4.D0*za(k1,k7)*za(k5,k7)*zb(k2,k7)*zb(k4,k7)*zab2(k3
     & ,k1,k2,k7)*zab2(k7,k1,k2,k6)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-1) + 1.D0/2.D0*za(k1,k7)*za(k5,k7)*zb(k2,k7)*zb(
     & k6,k7)*zab2(k3,k1,k2,k7)*zab2(k7,k1,k2,k4)*sd12x56*s17p27**(-3)*
     & ang7x56x12x7**(-1)*s12
      gr2 = gr2 + 1.D0/4.D0*za(k1,k7)*za(k5,k7)*zb(k2,k7)*zb(k6,k7)*
     & zab2(k3,k1,k2,k7)*zab2(k7,k1,k2,k4)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-1) - 1.D0/4.D0*za(k1,k7)*za(k5,k7)*zb(k2,k7)*
     & zab2(k3,k1,k2,k7)*zab2(k7,k1,k2,k4)*zab2(k7,k1,k2,k6)*sd12x56**2
     & *s17p27**(-3)*ang7x56x12x7**(-2) - 1.D0/4.D0*za(k1,k7)*zb(k2,k7)
     & *zb(k4,k7)*zb(k6,k7)*zab2(k3,k1,k2,k7)*zab2(k5,k1,k2,k7)*sd12x56
     & *s17p27**(-3)*sqr7x12x56x7**(-1)*s12 - 1.D0/4.D0*za(k1,k7)*zb(k2
     & ,k7)*zb(k4,k7)*zb(k6,k7)*zab2(k3,k1,k2,k7)*zab2(k5,k1,k2,k7)*
     & zab2(k7,k1,k2,k7)*sd12x56*s17p27**(-4)*sqr7x12x56x7**(-1)*s12 -
     & 1.D0/4.D0*za(k1,k7)*zb(k2,k7)*zb(k4,k7)*zb(k6,k7)*zab2(k3,k1,k2,
     & k7)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k7)*sd12x56*s17p27**(-3)*
     & sqr7x12x56x7**(-1)
      glr= - 1.D0/2.D0*za(k1,k2)*za(k1,k3)*za(k5,k7)*zb(k1,k2)*zb(k2,k7
     & )*zb(k4,k6)*zab2(k7,k1,k2,k7)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-1) - 1.D0/2.D0*za(k1,k2)*za(k1,k7)*za(k3,k5)*zb(
     & k1,k2)*zb(k2,k4)*zb(k6,k7)*sd12x56*s17p27**(-1)*
     & ang7x56x12x7**(-1) + 1.D0/2.D0*za(k1,k2)*za(k1,k7)*za(k3,k7)*za(
     & k5,k7)*zb(k1,k2)*zb(k2,k7)*zb(k4,k6)*sd12x56**2*s17p27**(-2)*
     & ang7x56x12x7**(-2) + za(k1,k2)*za(k3,k5)*zb(k2,k4)*zb(k2,k7)*zb(
     & k6,k7) - 1.D0/2.D0*za(k1,k2)*za(k3,k5)*zb(k2,k4)*zb(k2,k7)*zab2(
     & k7,k1,k2,k6)*sd12x56*s17p27**(-1)*ang7x56x12x7**(-1) + 1.D0/2.D0
     & *za(k1,k2)*za(k3,k5)*zb(k2,k7)**2*zb(k4,k7)*zb(k6,k7)*sd12x56*
     & s17p27**(-1)*sqr7x12x56x7**(-1) + 1.D0/2.D0*za(k1,k2)*za(k3,k7)*
     & za(k5,k7)*zb(k2,k7)**2*zb(k4,k6)*zab2(k7,k1,k2,k7)*sd12x56*
     & s17p27**(-3)*ang7x56x12x7**(-1)*s12 + 1.D0/2.D0*za(k1,k3)*za(k1,
     & k7)*za(k5,k7)*zb(k1,k2)*zb(k4,k6)*sd12x56**2*s17p27**(-1)*
     & ang7x56x12x7**(-2) + 1.D0/2.D0*za(k1,k3)*za(k1,k7)*zb(k1,k2)*zb(
     & k4,k6)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k7)*sd12x56*s17p27**(-2)*
     & ang7x56x12x7**(-1)
      glr = glr + 1.D0/2.D0*za(k1,k3)*za(k5,k7)*zb(k2,k7)*zb(k4,k6)*
     & sd12x56*s17p27**(-1)*ang7x56x12x7**(-1)*s12 - 1.D0/2.D0*za(k1,k7
     & )**2*za(k3,k5)*zb(k1,k2)*zb(k4,k7)*zb(k6,k7)*sd12x56*
     & s17p27**(-2)*ang7x56x12x7**(-1)*s12 + 1.D0/2.D0*za(k1,k7)**2*za(
     & k3,k5)*zb(k1,k2)*zb(k6,k7)*zab2(k7,k1,k2,k4)*sd12x56**2*
     & s17p27**(-2)*ang7x56x12x7**(-2) - 1.D0/2.D0*za(k1,k7)**2*za(k3,
     & k7)*zb(k1,k2)*zb(k4,k6)*zab2(k5,k1,k2,k7)*sd12x56**2*
     & s17p27**(-2)*ang7x56x12x7**(-2) + 1.D0/2.D0*za(k1,k7)*za(k3,k5)*
     & zb(k2,k4)*zb(k6,k7)*sd12x56*s17p27**(-1)*ang7x56x12x7**(-1)*s12
     &  + 1.D0/2.D0*za(k1,k7)*za(k3,k5)*zb(k2,k4)*zb(k6,k7)*sd12x56*
     & ang7x56x12x7**(-1) - 1.D0/2.D0*za(k1,k7)*za(k3,k5)*zb(k2,k4)*
     & zab2(k7,k1,k2,k6)*sd12x56**2*s17p27**(-1)*ang7x56x12x7**(-2) +
     & za(k1,k7)*za(k3,k5)*zb(k2,k7)*zb(k4,k7)*zb(k6,k7)*s17p27**(-1)*
     & s12 - 1.D0/2.D0*za(k1,k7)*za(k3,k5)*zb(k2,k7)*zb(k4,k7)*zab2(k7,
     & k1,k2,k6)*sd12x56*s17p27**(-2)*ang7x56x12x7**(-1)*s12 - 1.D0/2.D0
     & *za(k1,k7)*za(k3,k5)*zb(k2,k7)*zb(k6,k7)*zab2(k7,k1,k2,k4)*
     & sd12x56*s17p27**(-1)*ang7x56x12x7**(-1)
      glr = glr + 1.D0/2.D0*za(k1,k7)*za(k3,k5)*zb(k2,k7)*zab2(k7,k1,k2
     & ,k4)*zab2(k7,k1,k2,k6)*sd12x56**2*s17p27**(-2)*
     & ang7x56x12x7**(-2) - 1.D0/2.D0*za(k1,k7)*za(k3,k7)*zb(k2,k7)*zb(
     & k4,k6)*zab2(k5,k1,k2,k7)*zab2(k7,k1,k2,k7)*sd12x56*s17p27**(-3)*
     & ang7x56x12x7**(-1)*s12
      gl2=gl2/(s12*s34*s56)
      gr2=gr2/(s12*s34*s56)
      glr=glr/(s12*s34*s56)
      return
      end
