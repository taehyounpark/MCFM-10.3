!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine WpampLC_uppe_mp_paper(p1,p2,p3,p4,p5,p6,za,zb,amp)
c     arXiv:2105.00954, Eq.(C.10)
      implicit none
      include 'types.f'
      include 'cplx.h'
      include 'constants.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      integer p1,p2,p3,p4,p5,p6
      complex(dp):: amp,zab2,prp34,Alo,
     & Lsm1,Lsm1_2mh,I3m,Lnrat,Vpole,Wpamp_Tsum,Wpamp_TSL,Wpamp_Lmp
      real(dp):: t

c Statement functions
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      t(p1,p2,p3)=s(p1,p2)+s(p1,p3)+s(p2,p3)

      prp34=cone/cplx2(s(p3,p4)-wmass**2,wmass*wwidth)

      amp=
     & 1._dp/(t(p1,p5,p6))*((za(p3,p4)*zab2(p6,p1,p5,p4)**2*zab2(p5,p1,p6,p2)**2)/(za(p1,p6)*zab2(p6,p1,p5,p2)**3)
     & +(za(p2,p3)**2*zb(p1,p6)**2*zb(p3,p4))/(zab2(p2,p1,p6,p5)*zb(p1,p5))
     & )*Lsm1_2mh(s(p2,p6),t(p1,p5,p6),s(p1,p5),s(p3,p4))
     & -(za(p1,p2)**3*za(p3,p4)*zb(p4,p5)**2*t(p1,p2,p6)**2)/(za(p1,p6)*za(p6,p2)*zab2(p1,p2,p6,p5)**3*zab2(p2,p1,p6,p5))
     & *Lsm1_2mh(s(p1,p5),t(p1,p2,p6),s(p2,p6),s(p3,p4))
     & +(zb(p1,p6)**2*za(p2,p3)**2*zb(p3,p4))/(zb(p1,p5)*t(p1,p5,p6)*zab2(p2,p1,p6,p5))
     & *Lsm1(-s(p1,p5),-t(p1,p5,p6),-s(p1,p6),-t(p1,p5,p6))
     & -(zb(p1,p6)*zab2(p2,p1,p6,p4)*za(p2,p3))/(zb(p1,p5)*za(p1,p6)*zab2(p2,p1,p6,p5))
     & *Lsm1(-s(p1,p6),-t(p1,p2,p6),-s(p2,p6),-t(p1,p2,p6))
     & -(Wpamp_Tsum(p1,p2,p3,p4,p5,p6,za,zb)-Wpamp_TSL(p1,p2,p3,p4,p5,p6,za,zb)-Wpamp_TSL(p2,p1,p4,p3,p6,p5,zb,za))
     &  *I3m(s(p1,p5),s(p2,p6),s(p3,p4))
     & +Wpamp_Lmp(p1,p2,p3,p4,p5,p6,za,zb)

      Alo=-zb(p1,p6)*za(p2,p3)*zab2(p5,p2,p3,p4)/(zb(p1,p5)*za(p1,p6)*t(p2,p3,p4))

      Vpole=
     &  -1.5_dp*(lnrat(musq,-t(p1,p2,p6))+2._dp)
     &  -0.5_dp*(lnrat(musq,-s(p1,p6))**2+lnrat(musq,-s(p2,p6))**2)
     &  -2._dp*epinv*epinv2
     &  -epinv*(lnrat(musq,-s(p1,p6))+lnrat(musq,-s(p2,p6)))
     &  -1.5_dp*epinv

      amp=amp+Alo*Vpole

      amp = amp*prp34

      return
      end
