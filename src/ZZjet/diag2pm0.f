!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine diag2pm0(k1,k2,k3,k4,k5,k6,k7,za,zb,gl2,gr2,glr)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: k1,k2,k3,k4,k5,k6,k7
      complex(dp):: gl2,gr2,glr,zab2,zaba22,ang7x56x34x7
      real(dp):: sd34x56,s37p47,s347,s567,s12,s34,s56,s3
c--- statement functions
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
      zaba22(k1,k2,k3,k4,k5,k6)=
     & zab2(k1,k2,k3,k4)*za(k4,k6)+zab2(k1,k2,k3,k5)*za(k5,k6)
      s3(k1,k2,k3)=s(k1,k2)+s(k2,k3)+s(k3,k1)
c--- end statement functions
      s37p47=s(k3,k7)+s(k4,k7)
      s347=s3(k3,k4,k7)
      s567=s3(k5,k6,k7)
      s12=s(k1,k2)
      s34=s(k3,k4)
      s56=s(k5,k6)
      ang7x56x34x7=zaba22(k7,k5,k6,k3,k4,k7)
      sd34x56=s347*s567-s34*s56
      gl2=1.D0/4.D0*za(k1,k3)*za(k3,k7)*za(k5,k6)*zb(k2,k6)*zb(k3,k4)*
     & zb(k6,k7)*sd34x56*s37p47**(-1)*ang7x56x34x7**(-1)*s34 + 1.D0/4.D0
     & *za(k1,k3)*za(k3,k7)*za(k5,k6)*zb(k2,k6)*zb(k3,k4)*zb(k6,k7)*
     & sd34x56*ang7x56x34x7**(-1) - 1.D0/4.D0*za(k1,k3)*za(k3,k7)*za(k5
     & ,k6)*zb(k2,k6)*zb(k3,k4)*zab2(k7,k3,k4,k6)*sd34x56**2*
     & s37p47**(-1)*ang7x56x34x7**(-2) + 1.D0/4.D0*za(k1,k3)*za(k3,k7)*
     & za(k5,k7)*zb(k2,k7)*zb(k3,k4)*zb(k6,k7)*sd34x56*s37p47**(-2)*
     & ang7x56x34x7**(-1)*s34**2 + 1.D0/2.D0*za(k1,k3)*za(k3,k7)*za(k5,
     & k7)*zb(k2,k7)*zb(k3,k4)*zb(k6,k7)*sd34x56*s37p47**(-1)*
     & ang7x56x34x7**(-1)*s34 + 1.D0/4.D0*za(k1,k3)*za(k3,k7)*za(k5,k7)
     & *zb(k2,k7)*zb(k3,k4)*zb(k6,k7)*sd34x56*ang7x56x34x7**(-1) - 1.D0/
     & 4.D0*za(k1,k3)*za(k3,k7)*za(k5,k7)*zb(k2,k7)*zb(k3,k4)*zab2(k7,
     & k3,k4,k6)*sd34x56**2*s37p47**(-2)*ang7x56x34x7**(-2)*s34 - 1.D0/
     & 4.D0*za(k1,k3)*za(k3,k7)*za(k5,k7)*zb(k2,k7)*zb(k3,k4)*zab2(k7,
     & k3,k4,k6)*sd34x56**2*s37p47**(-1)*ang7x56x34x7**(-2) - 1.D0/4.D0
     & *za(k1,k3)*za(k3,k7)*za(k5,k7)*zb(k3,k4)*zb(k6,k7)*zab2(k7,k3,k4
     & ,k2)*sd34x56**2*s37p47**(-2)*ang7x56x34x7**(-2)*s34
      gl2 = gl2 - 1.D0/4.D0*za(k1,k3)*za(k3,k7)*za(k5,k7)*zb(k3,k4)*zb(
     & k6,k7)*zab2(k7,k3,k4,k2)*sd34x56**2*s37p47**(-1)*
     & ang7x56x34x7**(-2) + 1.D0/4.D0*za(k1,k3)*za(k3,k7)*za(k5,k7)*zb(
     & k3,k4)*zab2(k7,k3,k4,k2)*zab2(k7,k3,k4,k6)*sd34x56**3*
     & s37p47**(-2)*ang7x56x34x7**(-3) + 1.D0/4.D0*za(k1,k5)*za(k3,k4)*
     & za(k3,k7)*za(k5,k7)*zb(k2,k4)*zb(k3,k4)*zb(k5,k6)*sd34x56**2*
     & s37p47**(-1)*ang7x56x34x7**(-2) + 1.D0/4.D0*za(k1,k5)*za(k3,k4)*
     & za(k5,k7)*zb(k2,k4)*zb(k4,k7)*zb(k5,k6)*sd34x56*s37p47**(-1)*
     & ang7x56x34x7**(-1)*s34 - 1.D0/4.D0*za(k1,k5)*za(k3,k7)**2*za(k5,
     & k7)*zb(k2,k7)*zb(k3,k4)*zb(k5,k6)*sd34x56**2*s37p47**(-2)*
     & ang7x56x34x7**(-2)*s34 + 1.D0/4.D0*za(k1,k5)*za(k3,k7)**2*za(k5,
     & k7)*zb(k3,k4)*zb(k5,k6)*zab2(k7,k3,k4,k2)*sd34x56**3*
     & s37p47**(-2)*ang7x56x34x7**(-3) - 1.D0/4.D0*za(k1,k5)*za(k3,k7)*
     & za(k5,k7)*zb(k2,k7)*zb(k4,k7)*zb(k5,k6)*sd34x56*s37p47**(-2)*
     & ang7x56x34x7**(-1)*s34**2 + 1.D0/4.D0*za(k1,k5)*za(k3,k7)*za(k5,
     & k7)*zb(k4,k7)*zb(k5,k6)*zab2(k7,k3,k4,k2)*sd34x56**2*
     & s37p47**(-2)*ang7x56x34x7**(-2)*s34
      gl2 = gl2 - 1.D0/4.D0*za(k1,k7)*za(k3,k4)*za(k3,k7)*za(k5,k7)*zb(
     & k2,k4)*zb(k3,k4)*zb(k6,k7)*sd34x56**2*s37p47**(-2)*
     & ang7x56x34x7**(-2)*s34 - 1.D0/4.D0*za(k1,k7)*za(k3,k4)*za(k3,k7)
     & *za(k5,k7)*zb(k2,k4)*zb(k3,k4)*zb(k6,k7)*sd34x56**2*s37p47**(-1)
     & *ang7x56x34x7**(-2) + 1.D0/4.D0*za(k1,k7)*za(k3,k4)*za(k3,k7)*
     & za(k5,k7)*zb(k2,k4)*zb(k3,k4)*zab2(k7,k3,k4,k6)*sd34x56**3*
     & s37p47**(-2)*ang7x56x34x7**(-3) - 1.D0/4.D0*za(k1,k7)*za(k3,k4)*
     & za(k5,k7)*zb(k2,k4)*zb(k4,k7)*zb(k6,k7)*sd34x56*s37p47**(-2)*
     & ang7x56x34x7**(-1)*s34**2 - 1.D0/4.D0*za(k1,k7)*za(k3,k4)*za(k5,
     & k7)*zb(k2,k4)*zb(k4,k7)*zb(k6,k7)*sd34x56*s37p47**(-1)*
     & ang7x56x34x7**(-1)*s34 + 1.D0/4.D0*za(k1,k7)*za(k3,k4)*za(k5,k7)
     & *zb(k2,k4)*zb(k4,k7)*zab2(k7,k3,k4,k6)*sd34x56**2*s37p47**(-2)*
     & ang7x56x34x7**(-2)*s34 + 1.D0/4.D0*za(k1,k7)*za(k3,k7)**2*za(k5,
     & k6)*zb(k2,k6)*zb(k3,k4)*zb(k6,k7)*sd34x56**2*s37p47**(-2)*
     & ang7x56x34x7**(-2)*s34 + 1.D0/4.D0*za(k1,k7)*za(k3,k7)**2*za(k5,
     & k6)*zb(k2,k6)*zb(k3,k4)*zb(k6,k7)*sd34x56**2*s37p47**(-1)*
     & ang7x56x34x7**(-2)
      gl2 = gl2 - 1.D0/4.D0*za(k1,k7)*za(k3,k7)**2*za(k5,k6)*zb(k2,k6)*
     & zb(k3,k4)*zab2(k7,k3,k4,k6)*sd34x56**3*s37p47**(-2)*
     & ang7x56x34x7**(-3) + 1.D0/2.D0*za(k1,k7)*za(k3,k7)**2*za(k5,k7)*
     & zb(k2,k7)*zb(k3,k4)*zb(k6,k7)*sd34x56**2*s37p47**(-3)*
     & ang7x56x34x7**(-2)*s34**2 + 3.D0/4.D0*za(k1,k7)*za(k3,k7)**2*za(
     & k5,k7)*zb(k2,k7)*zb(k3,k4)*zb(k6,k7)*sd34x56**2*s37p47**(-2)*
     & ang7x56x34x7**(-2)*s34 + 1.D0/4.D0*za(k1,k7)*za(k3,k7)**2*za(k5,
     & k7)*zb(k2,k7)*zb(k3,k4)*zb(k6,k7)*sd34x56**2*s37p47**(-1)*
     & ang7x56x34x7**(-2) - 1.D0/2.D0*za(k1,k7)*za(k3,k7)**2*za(k5,k7)*
     & zb(k2,k7)*zb(k3,k4)*zab2(k7,k3,k4,k6)*sd34x56**3*s37p47**(-3)*
     & ang7x56x34x7**(-3)*s34 - 1.D0/4.D0*za(k1,k7)*za(k3,k7)**2*za(k5,
     & k7)*zb(k2,k7)*zb(k3,k4)*zab2(k7,k3,k4,k6)*sd34x56**3*
     & s37p47**(-2)*ang7x56x34x7**(-3) - 1.D0/2.D0*za(k1,k7)*za(k3,k7)
     & **2*za(k5,k7)*zb(k3,k4)*zb(k6,k7)*zab2(k7,k3,k4,k2)*sd34x56**3*
     & s37p47**(-3)*ang7x56x34x7**(-3)*s34 - 1.D0/2.D0*za(k1,k7)*za(k3,
     & k7)**2*za(k5,k7)*zb(k3,k4)*zb(k6,k7)*zab2(k7,k3,k4,k2)*
     & sd34x56**3*s37p47**(-2)*ang7x56x34x7**(-3)
      gl2 = gl2 + 1.D0/2.D0*za(k1,k7)*za(k3,k7)**2*za(k5,k7)*zb(k3,k4)*
     & zab2(k7,k3,k4,k2)*zab2(k7,k3,k4,k6)*sd34x56**4*s37p47**(-3)*
     & ang7x56x34x7**(-4) + 1.D0/4.D0*za(k1,k7)*za(k3,k7)*za(k5,k6)*zb(
     & k2,k6)*zb(k4,k7)*zb(k6,k7)*sd34x56*s37p47**(-2)*
     & ang7x56x34x7**(-1)*s34**2 + 1.D0/4.D0*za(k1,k7)*za(k3,k7)*za(k5,
     & k6)*zb(k2,k6)*zb(k4,k7)*zb(k6,k7)*sd34x56*s37p47**(-1)*
     & ang7x56x34x7**(-1)*s34 - 1.D0/4.D0*za(k1,k7)*za(k3,k7)*za(k5,k6)
     & *zb(k2,k6)*zb(k4,k7)*zab2(k7,k3,k4,k6)*sd34x56**2*s37p47**(-2)*
     & ang7x56x34x7**(-2)*s34 + 1.D0/2.D0*za(k1,k7)*za(k3,k7)*za(k5,k7)
     & *zb(k2,k7)*zb(k4,k7)*zb(k6,k7)*sd34x56*s37p47**(-3)*
     & ang7x56x34x7**(-1)*s34**3 + 3.D0/4.D0*za(k1,k7)*za(k3,k7)*za(k5,
     & k7)*zb(k2,k7)*zb(k4,k7)*zb(k6,k7)*sd34x56*s37p47**(-2)*
     & ang7x56x34x7**(-1)*s34**2 + 1.D0/4.D0*za(k1,k7)*za(k3,k7)*za(k5,
     & k7)*zb(k2,k7)*zb(k4,k7)*zb(k6,k7)*sd34x56*s37p47**(-1)*
     & ang7x56x34x7**(-1)*s34 - 1.D0/2.D0*za(k1,k7)*za(k3,k7)*za(k5,k7)
     & *zb(k2,k7)*zb(k4,k7)*zab2(k7,k3,k4,k6)*sd34x56**2*s37p47**(-3)*
     & ang7x56x34x7**(-2)*s34**2
      gl2 = gl2 - 1.D0/4.D0*za(k1,k7)*za(k3,k7)*za(k5,k7)*zb(k2,k7)*zb(
     & k4,k7)*zab2(k7,k3,k4,k6)*sd34x56**2*s37p47**(-2)*
     & ang7x56x34x7**(-2)*s34 - 1.D0/2.D0*za(k1,k7)*za(k3,k7)*za(k5,k7)
     & *zb(k4,k7)*zb(k6,k7)*zab2(k7,k3,k4,k2)*sd34x56**2*s37p47**(-3)*
     & ang7x56x34x7**(-2)*s34**2 - 1.D0/2.D0*za(k1,k7)*za(k3,k7)*za(k5,
     & k7)*zb(k4,k7)*zb(k6,k7)*zab2(k7,k3,k4,k2)*sd34x56**2*
     & s37p47**(-2)*ang7x56x34x7**(-2)*s34 + 1.D0/2.D0*za(k1,k7)*za(k3,
     & k7)*za(k5,k7)*zb(k4,k7)*zab2(k7,k3,k4,k2)*zab2(k7,k3,k4,k6)*
     & sd34x56**3*s37p47**(-3)*ang7x56x34x7**(-3)*s34
      gr2=1.D0/4.D0*za(k1,k3)*za(k3,k7)*za(k5,k6)*zb(k2,k6)*zb(k3,k4)*
     & zb(k6,k7)*sd34x56*s37p47**(-1)*ang7x56x34x7**(-1)*s34 + 1.D0/4.D0
     & *za(k1,k3)*za(k3,k7)*za(k5,k6)*zb(k2,k6)*zb(k3,k4)*zb(k6,k7)*
     & sd34x56*ang7x56x34x7**(-1) - 1.D0/4.D0*za(k1,k3)*za(k3,k7)*za(k5
     & ,k6)*zb(k2,k6)*zb(k3,k4)*zab2(k7,k3,k4,k6)*sd34x56**2*
     & s37p47**(-1)*ang7x56x34x7**(-2) + 1.D0/4.D0*za(k1,k3)*za(k3,k7)*
     & za(k5,k7)*zb(k2,k7)*zb(k3,k4)*zb(k6,k7)*sd34x56*s37p47**(-2)*
     & ang7x56x34x7**(-1)*s34**2 + 1.D0/2.D0*za(k1,k3)*za(k3,k7)*za(k5,
     & k7)*zb(k2,k7)*zb(k3,k4)*zb(k6,k7)*sd34x56*s37p47**(-1)*
     & ang7x56x34x7**(-1)*s34 + 1.D0/4.D0*za(k1,k3)*za(k3,k7)*za(k5,k7)
     & *zb(k2,k7)*zb(k3,k4)*zb(k6,k7)*sd34x56*ang7x56x34x7**(-1) - 1.D0/
     & 4.D0*za(k1,k3)*za(k3,k7)*za(k5,k7)*zb(k2,k7)*zb(k3,k4)*zab2(k7,
     & k3,k4,k6)*sd34x56**2*s37p47**(-2)*ang7x56x34x7**(-2)*s34 - 1.D0/
     & 4.D0*za(k1,k3)*za(k3,k7)*za(k5,k7)*zb(k2,k7)*zb(k3,k4)*zab2(k7,
     & k3,k4,k6)*sd34x56**2*s37p47**(-1)*ang7x56x34x7**(-2) - 1.D0/4.D0
     & *za(k1,k3)*za(k3,k7)*za(k5,k7)*zb(k3,k4)*zb(k6,k7)*zab2(k7,k3,k4
     & ,k2)*sd34x56**2*s37p47**(-2)*ang7x56x34x7**(-2)*s34
      gr2 = gr2 - 1.D0/4.D0*za(k1,k3)*za(k3,k7)*za(k5,k7)*zb(k3,k4)*zb(
     & k6,k7)*zab2(k7,k3,k4,k2)*sd34x56**2*s37p47**(-1)*
     & ang7x56x34x7**(-2) + 1.D0/4.D0*za(k1,k3)*za(k3,k7)*za(k5,k7)*zb(
     & k3,k4)*zab2(k7,k3,k4,k2)*zab2(k7,k3,k4,k6)*sd34x56**3*
     & s37p47**(-2)*ang7x56x34x7**(-3) + 1.D0/4.D0*za(k1,k5)*za(k3,k4)*
     & za(k3,k7)*za(k5,k7)*zb(k2,k4)*zb(k3,k4)*zb(k5,k6)*sd34x56**2*
     & s37p47**(-1)*ang7x56x34x7**(-2) + 1.D0/4.D0*za(k1,k5)*za(k3,k4)*
     & za(k5,k7)*zb(k2,k4)*zb(k4,k7)*zb(k5,k6)*sd34x56*s37p47**(-1)*
     & ang7x56x34x7**(-1)*s34 - 1.D0/4.D0*za(k1,k5)*za(k3,k7)**2*za(k5,
     & k7)*zb(k2,k7)*zb(k3,k4)*zb(k5,k6)*sd34x56**2*s37p47**(-2)*
     & ang7x56x34x7**(-2)*s34 + 1.D0/4.D0*za(k1,k5)*za(k3,k7)**2*za(k5,
     & k7)*zb(k3,k4)*zb(k5,k6)*zab2(k7,k3,k4,k2)*sd34x56**3*
     & s37p47**(-2)*ang7x56x34x7**(-3) - 1.D0/4.D0*za(k1,k5)*za(k3,k7)*
     & za(k5,k7)*zb(k2,k7)*zb(k4,k7)*zb(k5,k6)*sd34x56*s37p47**(-2)*
     & ang7x56x34x7**(-1)*s34**2 + 1.D0/4.D0*za(k1,k5)*za(k3,k7)*za(k5,
     & k7)*zb(k4,k7)*zb(k5,k6)*zab2(k7,k3,k4,k2)*sd34x56**2*
     & s37p47**(-2)*ang7x56x34x7**(-2)*s34
      gr2 = gr2 - 1.D0/4.D0*za(k1,k7)*za(k3,k4)*za(k3,k7)*za(k5,k7)*zb(
     & k2,k4)*zb(k3,k4)*zb(k6,k7)*sd34x56**2*s37p47**(-2)*
     & ang7x56x34x7**(-2)*s34 - 1.D0/4.D0*za(k1,k7)*za(k3,k4)*za(k3,k7)
     & *za(k5,k7)*zb(k2,k4)*zb(k3,k4)*zb(k6,k7)*sd34x56**2*s37p47**(-1)
     & *ang7x56x34x7**(-2) + 1.D0/4.D0*za(k1,k7)*za(k3,k4)*za(k3,k7)*
     & za(k5,k7)*zb(k2,k4)*zb(k3,k4)*zab2(k7,k3,k4,k6)*sd34x56**3*
     & s37p47**(-2)*ang7x56x34x7**(-3) - 1.D0/4.D0*za(k1,k7)*za(k3,k4)*
     & za(k5,k7)*zb(k2,k4)*zb(k4,k7)*zb(k6,k7)*sd34x56*s37p47**(-2)*
     & ang7x56x34x7**(-1)*s34**2 - 1.D0/4.D0*za(k1,k7)*za(k3,k4)*za(k5,
     & k7)*zb(k2,k4)*zb(k4,k7)*zb(k6,k7)*sd34x56*s37p47**(-1)*
     & ang7x56x34x7**(-1)*s34 + 1.D0/4.D0*za(k1,k7)*za(k3,k4)*za(k5,k7)
     & *zb(k2,k4)*zb(k4,k7)*zab2(k7,k3,k4,k6)*sd34x56**2*s37p47**(-2)*
     & ang7x56x34x7**(-2)*s34 + 1.D0/4.D0*za(k1,k7)*za(k3,k7)**2*za(k5,
     & k6)*zb(k2,k6)*zb(k3,k4)*zb(k6,k7)*sd34x56**2*s37p47**(-2)*
     & ang7x56x34x7**(-2)*s34 + 1.D0/4.D0*za(k1,k7)*za(k3,k7)**2*za(k5,
     & k6)*zb(k2,k6)*zb(k3,k4)*zb(k6,k7)*sd34x56**2*s37p47**(-1)*
     & ang7x56x34x7**(-2)
      gr2 = gr2 - 1.D0/4.D0*za(k1,k7)*za(k3,k7)**2*za(k5,k6)*zb(k2,k6)*
     & zb(k3,k4)*zab2(k7,k3,k4,k6)*sd34x56**3*s37p47**(-2)*
     & ang7x56x34x7**(-3) + 1.D0/2.D0*za(k1,k7)*za(k3,k7)**2*za(k5,k7)*
     & zb(k2,k7)*zb(k3,k4)*zb(k6,k7)*sd34x56**2*s37p47**(-3)*
     & ang7x56x34x7**(-2)*s34**2 + 3.D0/4.D0*za(k1,k7)*za(k3,k7)**2*za(
     & k5,k7)*zb(k2,k7)*zb(k3,k4)*zb(k6,k7)*sd34x56**2*s37p47**(-2)*
     & ang7x56x34x7**(-2)*s34 + 1.D0/4.D0*za(k1,k7)*za(k3,k7)**2*za(k5,
     & k7)*zb(k2,k7)*zb(k3,k4)*zb(k6,k7)*sd34x56**2*s37p47**(-1)*
     & ang7x56x34x7**(-2) - 1.D0/2.D0*za(k1,k7)*za(k3,k7)**2*za(k5,k7)*
     & zb(k2,k7)*zb(k3,k4)*zab2(k7,k3,k4,k6)*sd34x56**3*s37p47**(-3)*
     & ang7x56x34x7**(-3)*s34 - 1.D0/4.D0*za(k1,k7)*za(k3,k7)**2*za(k5,
     & k7)*zb(k2,k7)*zb(k3,k4)*zab2(k7,k3,k4,k6)*sd34x56**3*
     & s37p47**(-2)*ang7x56x34x7**(-3) - 1.D0/2.D0*za(k1,k7)*za(k3,k7)
     & **2*za(k5,k7)*zb(k3,k4)*zb(k6,k7)*zab2(k7,k3,k4,k2)*sd34x56**3*
     & s37p47**(-3)*ang7x56x34x7**(-3)*s34 - 1.D0/2.D0*za(k1,k7)*za(k3,
     & k7)**2*za(k5,k7)*zb(k3,k4)*zb(k6,k7)*zab2(k7,k3,k4,k2)*
     & sd34x56**3*s37p47**(-2)*ang7x56x34x7**(-3)
      gr2 = gr2 + 1.D0/2.D0*za(k1,k7)*za(k3,k7)**2*za(k5,k7)*zb(k3,k4)*
     & zab2(k7,k3,k4,k2)*zab2(k7,k3,k4,k6)*sd34x56**4*s37p47**(-3)*
     & ang7x56x34x7**(-4) + 1.D0/4.D0*za(k1,k7)*za(k3,k7)*za(k5,k6)*zb(
     & k2,k6)*zb(k4,k7)*zb(k6,k7)*sd34x56*s37p47**(-2)*
     & ang7x56x34x7**(-1)*s34**2 + 1.D0/4.D0*za(k1,k7)*za(k3,k7)*za(k5,
     & k6)*zb(k2,k6)*zb(k4,k7)*zb(k6,k7)*sd34x56*s37p47**(-1)*
     & ang7x56x34x7**(-1)*s34 - 1.D0/4.D0*za(k1,k7)*za(k3,k7)*za(k5,k6)
     & *zb(k2,k6)*zb(k4,k7)*zab2(k7,k3,k4,k6)*sd34x56**2*s37p47**(-2)*
     & ang7x56x34x7**(-2)*s34 + 1.D0/2.D0*za(k1,k7)*za(k3,k7)*za(k5,k7)
     & *zb(k2,k7)*zb(k4,k7)*zb(k6,k7)*sd34x56*s37p47**(-3)*
     & ang7x56x34x7**(-1)*s34**3 + 3.D0/4.D0*za(k1,k7)*za(k3,k7)*za(k5,
     & k7)*zb(k2,k7)*zb(k4,k7)*zb(k6,k7)*sd34x56*s37p47**(-2)*
     & ang7x56x34x7**(-1)*s34**2 + 1.D0/4.D0*za(k1,k7)*za(k3,k7)*za(k5,
     & k7)*zb(k2,k7)*zb(k4,k7)*zb(k6,k7)*sd34x56*s37p47**(-1)*
     & ang7x56x34x7**(-1)*s34 - 1.D0/2.D0*za(k1,k7)*za(k3,k7)*za(k5,k7)
     & *zb(k2,k7)*zb(k4,k7)*zab2(k7,k3,k4,k6)*sd34x56**2*s37p47**(-3)*
     & ang7x56x34x7**(-2)*s34**2
      gr2 = gr2 - 1.D0/4.D0*za(k1,k7)*za(k3,k7)*za(k5,k7)*zb(k2,k7)*zb(
     & k4,k7)*zab2(k7,k3,k4,k6)*sd34x56**2*s37p47**(-2)*
     & ang7x56x34x7**(-2)*s34 - 1.D0/2.D0*za(k1,k7)*za(k3,k7)*za(k5,k7)
     & *zb(k4,k7)*zb(k6,k7)*zab2(k7,k3,k4,k2)*sd34x56**2*s37p47**(-3)*
     & ang7x56x34x7**(-2)*s34**2 - 1.D0/2.D0*za(k1,k7)*za(k3,k7)*za(k5,
     & k7)*zb(k4,k7)*zb(k6,k7)*zab2(k7,k3,k4,k2)*sd34x56**2*
     & s37p47**(-2)*ang7x56x34x7**(-2)*s34 + 1.D0/2.D0*za(k1,k7)*za(k3,
     & k7)*za(k5,k7)*zb(k4,k7)*zab2(k7,k3,k4,k2)*zab2(k7,k3,k4,k6)*
     & sd34x56**3*s37p47**(-3)*ang7x56x34x7**(-3)*s34
      glr= 0
      gl2=gl2/(s12*s34*s56)
      gr2=gr2/(s12*s34*s56)
      glr=glr/(s12*s34*s56)
      return
      end
