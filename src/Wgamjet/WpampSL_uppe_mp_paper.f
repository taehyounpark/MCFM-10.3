!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine WpampSL_uppe_mp_paper(p1,p2,p3,p4,p5,p6,za,zb,amp)
c     arXiv:2105.00954, Eq.(C.17)
      implicit none
      include 'types.f'
      include 'cplx.h'
      include 'constants.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
c      include 'verbose.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
c      include 'uselfunctions.f'
      integer p1,p2,p3,p4,p5,p6
      complex(dp):: amp,zab2,prp34,Alo,
     & L0,L1,Lsm1,Lsm1_2mht,Lsm1_2mh,I3m,Lnrat,Vpole,
     & Wpamp_Tsym,Wpamp_Tantisym,Wpamp_TSL,Wpamp_Lmp
      real(dp):: t,Delta3,DeltaTM1234,DeltaTM1526

c Statement functions
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      t(p1,p2,p3)=s(p1,p2)+s(p1,p3)+s(p2,p3)
      Delta3(p1,p2,p3,p4,p5,p6)=
     & (s(p1,p2)**2+s(p3,p4)**2+s(p5,p6)**2
     & -2._dp*(s(p1,p2)*s(p3,p4)+s(p1,p2)*s(p5,p6)+s(p3,p4)*s(p5,p6)))

      DeltaTM1234=Delta3(p1,p2,p3,p4,p5,p6)
      DeltaTM1526=Delta3(p1,p5,p2,p6,p3,p4)

      prp34=cone/cplx2(s(p3,p4)-wmass**2,wmass*wwidth)

      amp=
     & +(zb(p1,p2)**2*zab2(p3,p1,p2,p5)**2*zb(p3,p4))
     & /(zb(p1,p5)*zb(p2,p5)**2*zab2(p6,p1,p2,p5)*zab2(p6,p1,p5,p2))*Lsm1(-s(p1,p5),-t(p1,p2,p5),-s(p1,p2),-t(p1,p2,p5))
     & -(za(p1,p2)*za(p3,p6)*s(p1,p2)*zab2(p6,p1,p2,p4))
     & /(za(p1,p6)**3*zb(p1,p5)*zab2(p6,p1,p2,p5))*Lsm1(-s(p1,p2),-t(p1,p2,p6),-s(p2,p6),-t(p1,p2,p6))
     & +(za(p2,p3)*zb(p1,p6)*zab2(p2,p1,p6,p4))
     & /(za(p1,p6)*zb(p1,p5)*zab2(p2,p1,p6,p5))*Lsm1(-s(p1,p2),-t(p1,p2,p6),-s(p1,p6),-t(p1,p2,p6))
     & -(za(p5,p6)**2*zab2(p1,p5,p6,p4)**2*za(p3,p4))
     & /(za(p1,p6)**3*zab2(p6,p1,p5,p2)*t(p1,p5,p6))*Lsm1(-s(p1,p5),-t(p1,p5,p6),-s(p5,p6),-t(p1,p5,p6))
     & -(za(p2,p3)**2*zb(p1,p6)**2*zb(p3,p4))
     & /(zb(p1,p5)*zab2(p2,p1,p6,p5)*t(p1,p5,p6))*Lsm1(-s(p1,p6),-t(p1,p5,p6),-s(p5,p6),-t(p1,p5,p6))
     & +1._dp/(t(p1,p5,p6))*((za(p1,p5)*zab2(p1,p5,p6,p4)**2*zab2(p5,p1,p6,p2)**2*za(p3,p4))
     & /(za(p1,p6)*za(p5,p6)*zab2(p1,p5,p6,p2)**3)
     & -(za(p2,p3)**2*zb(p1,p6)**3*zb(p3,p4))/(zb(p1,p5)*zb(p5,p6)*zab2(p2,p5,p6,p1))
     & +(za(p2,p3)**2*zb(p1,p6)**2*zab2(p2,p1,p5,p6)*zb(p3,p4))/(zb(p5,p6)*zab2(p2,p1,p6,p5)*zab2(p2,p5,p6,p1))
     & -(zab2(p1,p5,p6,p4)**2*zab2(p5,p1,p6,p2)**3*za(p3,p4))/(za(p5,p6)*zab2(p1,p5,p6,p2)**3*zab2(p6,p1,p5,p2))
     & )
     &  *Lsm1_2mht(s(p1,p2),t(p1,p5,p6),s(p5,p6),s(p3,p4))
      amp=amp
     & +(zab2(p3,p1,p2,p5)**2*zab2(p6,p2,p5,p1)**2*zb(p3,p4))
     & /(zb(p1,p5)*zab2(p6,p1,p2,p5)**3*zab2(p6,p1,p5,p2))*Lsm1_2mht(s(p5,p6),t(p1,p2,p5),s(p1,p2),s(p3,p4))
     & +(za(p2,p6)**2*za(p3,p4)*zb(p4,p5)**2*t(p1,p2,p6)**2)
     & /(za(p1,p6)*zab2(p2,p1,p6,p5)*zab2(p6,p1,p2,p5)**3)*Lsm1_2mht(s(p5,p6),t(p1,p2,p6),s(p1,p2),s(p3,p4))
     & -(za(p1,p2)**2*za(p3,p4)*zb(p4,p5)**2*t(p1,p2,p6)**2)
     & /(za(p2,p6)*zab2(p1,p2,p6,p5)**3*zab2(p6,p1,p2,p5))*Lsm1_2mh(s(p1,p5),t(p1,p2,p6),s(p2,p6),s(p3,p4))
     & +(za(p3,p6)**2*zb(p1,p2)**2*t(p1,p2,p5)**2*zb(p3,p4))
     & /(zb(p1,p5)*zab2(p6,p1,p2,p5)*zab2(p6,p1,p5,p2)**3)*Lsm1_2mh(s(p2,p6),t(p1,p2,p5),s(p1,p5),s(p3,p4))

     & -(Wpamp_Tantisym(p1,p2,p3,p4,p5,p6,za,zb)-Wpamp_Tantisym(p2,p1,p4,p3,p6,p5,zb,za)
     &  +Wpamp_Tsym(p1,p2,p3,p4,p5,p6,za,zb)+Wpamp_Tsym(p2,p1,p4,p3,p6,p5,zb,za))*I3m(s(p1,p2),s(p3,p4),s(p5,p6))

     & -(Wpamp_TSL(p1,p2,p3,p4,p5,p6,za,zb)+Wpamp_TSL(p2,p1,p4,p3,p6,p5,zb,za))*I3m(s(p1,p5),s(p2,p6),s(p3,p4))

      amp=amp
     & +(za(p1,p3)*zb(p2,p6)*zb(p4,p3)*zab2(p3,p1,p4,p6))/(zb(p2,p5)*zab2(p1,p2,p6,p5)*zab2(p1,p3,p4,p2))
     &  *Lnrat(-s(p1,p2),-s(p1,p5))
     & +(za(p1,p3)*za(p2,p3)*zb(p4,p3)*t(p3,p4,p5))/(za(p1,p6)*zab2(p1,p2,p6,p5)*zab2(p6,p1,p2,p5))
     &  *Lnrat(-s(p1,p5),-s(p3,p4))
     & -((za(p1,p2)*za(p3,p6)*zb(p4,p6)*t(p3,p4,p5))/(za(p1,p6)*zab2(p1,p2,p6,p5)*zab2(p6,p1,p2,p5))
     & +(za(p1,p2)*za(p3,p5)*zb(p4,p6))/(za(p1,p6)*zab2(p1,p2,p6,p5))
     & -(za(p1,p3)*za(p2,p6)*zb(p1,p6)*zb(p4,p6))/(zb(p1,p5)*za(p1,p6)*zab2(p1,p2,p6,p5))
     & +(za(p3,p4)*zb(p2,p6)*zb(p3,p4)*zab2(p3,p1,p2,p4))/(zb(p2,p5)*zab2(p1,p3,p4,p2)*zab2(p6,p1,p2,p5))
     & +(za(p2,p6)*za(p3,p5)*zb(p4,p6))/(za(p1,p6)*zab2(p6,p1,p2,p5))
     & +(za(p3,p4)*zb(p4,p6)*zab2(p2,p1,p6,p5)*zab2(p6,p1,p2,p4))/(zab2(p1,p2,p6,p5)*zab2(p6,p1,p2,p5)**2)
     &  )*Lnrat(-s(p1,p2),-s(p3,p4))
      amp=amp
     & -(za(p1,p3)*za(p2,p6)*zb(p4,p6)*t(p3,p4,p5))/(za(p1,p6)*zab2(p1,p2,p6,p5)*zab2(p6,p1,p2,p5))
     &  *Lnrat(-s(p3,p4),-t(p1,p2,p6))
     & -(za(p1,p3)*zb(p1,p4)*t(p3,p4,p5)*zab2(p2,p1,p6,p5))/(zb(p1,p5)*za(p1,p6)*zab2(p1,p2,p6,p5)*zab2(p6,p1,p2,p5))
     &  *Lnrat(-s(p1,p2),-t(p1,p2,p6))
     & -2*Lnrat(-s(p3,p4),-s(p1,p5))*(za(p1,p2)*zb(p2,p6)*zab2(p3,p1,p5,p4)*zab2(p5,p2,p6,p1)*s(p3,p4))
     & /(DeltaTM1526*zab2(p1,p2,p6,p5)*zab2(p6,p1,p5,p2))
     & -Lnrat(-s(p1,p5),-s(p2,p6))*(zab2(p1,p2,p5,p6)*zab2(p3,p1,p5,p4)*zab2(p5,p2,p6,p1)*s(p3,p4))
     & /(DeltaTM1526*zab2(p1,p2,p6,p5)*zab2(p6,p1,p5,p2))
     & +Lnrat(-s(p1,p5),-s(p2,p6))*((s(p1,p5)-s(p2,p6))*zab2(p1,p2,p5,p6)*zab2(p3,p1,p5,p4)*zab2(p5,p2,p6,p1))
     & /(DeltaTM1526*zab2(p1,p2,p6,p5)*zab2(p6,p1,p5,p2))
     & -Lnrat(-s(p3,p4),-s(p2,p6))*(s(p3,p4)*(s(p3,p4)-s(p1,p5)-s(p2,p6))*zb(p5,p6)*zab2(p3,p1,p5,p4)*zab2(p5,p2,p6,p1))
     & /(zb(p1,p5)*DeltaTM1526*zab2(p1,p2,p6,p5)*zab2(p6,p1,p5,p2))
     & -Lnrat(-s(p5,p6),-s(p1,p2))*
     & (zab2(p3,p1,p2,p4)*zab2(p5,p1,p2,p6)*(2*s(p3,p4)*s(p5,p6)-(s(p3,p4)+s(p5,p6)-s(p1,p2))*t(p2,p3,p4)))
     & /(zab2(p1,p3,p4,p2)*zab2(p6,p1,p2,p5)*DeltaTM1234)
      amp=amp
     & +Lnrat(-s(p1,p2),-s(p3,p4))*(s(p3,p4)*(t(p2,p3,p4)-t(p1,p3,p4))*zab2(p3,p1,p2,p4)*zab2(p5,p1,p2,p6))
     & /(zab2(p1,p3,p4,p2)*zab2(p6,p1,p2,p5)*DeltaTM1234)
     & -((za(p1,p3)*za(p5,p6)*zb(p2,p4)*zab2(p2,p1,p5,p6))/(za(p1,p6)*zab2(p1,p2,p6,p5)*zab2(p6,p1,p5,p2))
     & -(za(p1,p3)*za(p2,p5)*zb(p4,p6))/(za(p1,p6)*zab2(p1,p2,p6,p5))
     & +(za(p1,p3)*za(p2,p3)*zb(p1,p2)*zb(p3,p4))/(zb(p2,p5)*za(p1,p6)*zab2(p6,p1,p2,p5))
     & +(za(p3,p5)*za(p5,p6)*zb(p4,p6))/(za(p1,p6)*zab2(p6,p1,p5,p2))
     &  )*Lnrat(-s(p1,p5),-s(p3,p4))
     & -((za(p1,p3)*za(p2,p6)*zb(p1,p2)*zb(p2,p4)*zab2(p2,p1,p5,p6))/(zb(p1,p5)*za(p1,p6)*zab2(p1,p2,p6,p5)*zab2(p6,p1,p5,p2))
     & +(za(p3,p6)*za(p5,p6)*zb(p1,p6)*zb(p4,p6))/(zb(p1,p5)*za(p1,p6)*zab2(p6,p1,p5,p2))
     & -(za(p3,p6)*za(p5,p6)*zb(p1,p2)*zb(p2,p4)*zab2(p2,p1,p5,p6))/(zb(p1,p5)*za(p1,p6)*zab2(p6,p1,p5,p2)**2)
     & )*Lnrat(-s(p2,p6),-s(p3,p4))
      amp=amp
     & +(za(p3,p6)*zab2(p5,p1,p2,p6)*(zab2(p1,p5,p6,p1)*zb(p4,p5)+zab2(p6,p1,p2,p4)*zb(p5,p6)))
     & /(zab2(p1,p3,p4,p2)*zab2(p6,p1,p2,p5)**2)
     &  *Lnrat(-s(p1,p2),-s(p5,p6))
     & -(za(p1,p2)*za(p2,p6)*za(p1,p3)*zb(p1,p6)*zb(p2,p4))/(zb(p1,p5)*za(p1,p6)**2*zab2(p1,p2,p6,p5))
     &  *Lnrat(-s(p1,p2),-s(p2,p6))
     & -(za(p1,p5)*za(p3,p4)*zab2(p1,p2,p3,p4)*zab2(p5,p1,p6,p4))/(za(p1,p6)**2*zab2(p1,p3,p4,p2)*t(p2,p3,p4))
     &  *Lnrat(-s(p1,p5),-s(p5,p6))
     & -(za(p1,p5)*za(p3,p4)*za(p5,p6)*zb(p2,p4)*zb(p4,p6))/(za(p1,p6)*zab2(p1,p3,p4,p2)*zab2(p6,p1,p5,p2))
     &  *Lnrat(-t(p1,p5,p6),-s(p3,p4))
     & -(za(p3,p6)*zb(p1,p5)*zb(p3,p4)*zab2(p3,p1,p2,p6))/(zab2(p6,p1,p2,p5)**2*zb(p2,p5))
     &  *Lnrat(-t(p1,p2,p5),-s(p1,p2))
     & +(zb(p3,p4))/(zb(p2,p5))*(-(za(p3,p6)*zb(p1,p2)**2*zab2(p3,p1,p5,p6))/(zb(p1,p5)*zab2(p6,p1,p5,p2)**2)
     & -(zb(p1,p2)*zb(p2,p4)*zab2(p3,p2,p4,p1)*za(p3,p4))/(zb(p1,p5)*zab2(p6,p1,p5,p2)**2)
     & +(zb(p4,p5)*zab2(p3,p2,p6,p1)*za(p3,p4))/(zab2(p6,p1,p2,p5)**2)
     &  )*Lnrat(-s(p3,p4),-t(p1,p2,p5))
      amp=amp
     & -((zb(p4,p5)*t(p3,p4,p5)*zab2(p2,p1,p6,p5)*zab2(p5,p1,p2,p4))/(zb(p3,p4)*zab2(p1,p2,p6,p5)*zab2(p6,p1,p2,p5)**2)
     & +(za(p2,p5)*zb(p4,p6)*zb(p4,p5)*t(p3,p4,p5))/(zb(p3,p4)*zab2(p1,p2,p6,p5)*zab2(p6,p1,p2,p5))
     &  )*L0(-t(p1,p2,p6),-s(p3,p4))
     & -((zb(p2,p4)*zab2(p2,p1,p5,p4)*zab2(p5,p1,p6,p2)**2)/(zb(p3,p4)*zab2(p1,p3,p4,p2)*zab2(p6,p1,p5,p2)**2)
     & -(za(p2,p5)*zb(p2,p4)*zb(p4,p6)*zab2(p5,p1,p6,p2))/(zb(p3,p4)*zab2(p1,p3,p4,p2)*zab2(p6,p1,p5,p2))
     &  )*L0(-t(p1,p5,p6),-s(p3,p4))
     & -(za(p1,p5)**2*za(p3,p4)*zb(p1,p4)*zab2(p1,p2,p3,p4))/(zb(p5,p6)*za(p1,p6)**2*za(p5,p6)*zab2(p1,p3,p4,p2))
     &  *L0(-t(p1,p5,p6),-s(p5,p6))
     & +s(p3,p4)*((za(p3,p6)**2*zb(p1,p6)*zb(p3,p4))/(zab2(p6,p1,p2,p5)**2*zab2(p6,p1,p5,p2))
     & -(za(p3,p6)*zb(p1,p6)*zb(p2,p4)*zab2(p6,p2,p5,p1))/(zb(p1,p5)*zab2(p6,p1,p2,p5)*zab2(p6,p1,p5,p2)**2)
     &  )*L0(-t(p1,p2,p5),-s(p3,p4))
      amp=amp
     & -((za(p1,p3)*za(p2,p6)*zb(p1,p6)*t(p3,p4,p5)*zab2(p6,p1,p2,p4))
     & /(zb(p1,p2)*zb(p1,p5)*za(p1,p2)*za(p1,p6)**2*zab2(p6,p1,p2,p5))
     & -(za(p2,p6)*za(p3,p6)*zb(p4,p6)*t(p3,p4,p5))/(zab2(p6,p1,p2,p5)**2*za(p1,p6))
     &  )*L0(-t(p1,p2,p6),-s(p1,p2))
     & -(za(p3,p5)*zb(p1,p5)*zb(p3,p4)*zab2(p3,p1,p2,p5))/(zab2(p6,p1,p2,p5)**2*zb(p2,p5))
     &  *L0(-t(p1,p2,p5),-s(p1,p2))
     & -(za(p2,p3)*zb(p1,p2)**2*zb(p3,p4)*zab2(p3,p1,p5,p2))/(zb(p1,p5)*zb(p2,p5)*zab2(p6,p1,p5,p2)**2)
     &  *L0(-t(p1,p2,p5),-s(p1,p5))
     & -((za(p3,p4)*za(p5,p6)*zb(p1,p2)*zb(p4,p6)*zab2(p6,p1,p5,p4))/(zb(p1,p5)*za(p1,p6)*zab2(p6,p1,p5,p2)**2)
     & -(za(p3,p4)*za(p5,p6)*zb(p4,p6)*zab2(p6,p1,p5,p4))/(zb(p1,p5)*za(p1,p6)**2*zab2(p6,p1,p5,p2))
     &  )*L0(-t(p1,p5,p6),-s(p1,p5))
      amp=amp
     & +(za(p3,p6)*zb(p1,p6)*zb(p4,p6)*zab2(p6,p2,p5,p1))/(zb(p1,p5)*zab2(p6,p1,p2,p5)*zab2(p6,p1,p5,p2))
     &  *L1(-t(p1,p2,p5),-s(p3,p4))
     & -(za(p3,p4)*za(p5,p6)**2*zb(p4,p6)**2)/(zb(p1,p5)*za(p1,p6)*za(p1,p5)*zab2(p6,p1,p5,p2))
     &  *L1(-t(p1,p5,p6),-s(p1,p5))
     & +(za(p2,p6)*za(p3,p6)*zb(p1,p6)*zb(p4,p6)*t(p3,p4,p5))/(zb(p1,p2)*zb(p1,p5)*za(p1,p2)*za(p1,p6)*zab2(p6,p1,p2,p5))
     &  *L1(-t(p1,p2,p6),-s(p1,p2))
     & -(za(p3,p4)*za(p5,p6)*zb(p1,p4)*zb(p4,p6))/(zb(p1,p5)*za(p1,p6)*zab2(p6,p1,p5,p2))
     & -Wpamp_Lmp(p1,p2,p3,p4,p5,p6,za,zb)

      Alo=zb(p1,p6)*za(p2,p3)*zab2(p5,p2,p3,p4)
     & /(zb(p1,p5)*za(p1,p6)*t(p2,p3,p4))

      Vpole=
     &  -0.5_dp*lnrat(musq,-s(p1,p2))**2
     &  -epinv*epinv2
     &  -epinv*lnrat(musq,-s(p1,p2))
     &  -1.5_dp*(epinv+lnrat(musq,-t(p1,p2,p6))+2._dp)

      amp=amp+Alo*Vpole

      amp=amp*prp34

c      write(6,*) 'SL Qd mp bubrat',
c     &+(b256mp_sl_Wpgamg(p1,p2,p3,p4,p5,p6,za,zb)
c     & +b256pm_sl_Wpgamg(p1,p2,p3,p4,p6,p5,za,zb))*lnrat(-s(p3,p4),-t(p2,p5,p6))
c     &+b125pm_sl_Wpgamg(p2,p1,p4,p3,p5,p6,zb,za)*lnrat(-s(p3,p4),-t(p1,p2,p5))
c     &+b56mp_sl_Qd_Wpgamg(p1,p2,p3,p4,p5,p6,za,zb)*lnrat(-s(p3,p4),-s(p5,p6))
c     &+(b126mp_sl_Wpgamg(p1,p2,p3,p4,p5,p6,za,zb)
c     & +b125pm_sl_Wpgamg(p1,p2,p3,p4,p6,p5,za,zb)
c     & -b126mp_sl_Qu_Wpgamg(p1,p2,p3,p4,p5,p6,za,zb))*lnrat(-s(p3,p4),-t(p1,p2,p6))
c     &+b12mp_sl_Qd_Wpgamg(p1,p2,p3,p4,p5,p6,za,zb)*lnrat(-s(p3,p4),-s(p1,p2))
c     &+b26mp_sl_Qd_Wpgamg(p1,p2,p3,p4,p5,p6,za,zb)*lnrat(-s(p3,p4),-s(p2,p6))
c     &+ratpm_sl_Qu_Wpgamg(p2,p1,p4,p3,p5,p6,zb,za)
c     & -1.5_dp*(lnrat(musq,-s(p3,p4))+2._dp)*Alo

      return
      end
