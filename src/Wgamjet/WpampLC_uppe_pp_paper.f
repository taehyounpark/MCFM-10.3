!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine WpampLC_uppe_pp_paper(p1,p2,p3,p4,p5,p6,za,zb,amp)
c--- Version from paper arXiv:2105.00954, Eq.(C.2)
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
      complex(dp):: amp,zab2,prp34,Alo,L0,L1,Lsm1,Lsm1_2me,Lnrat,Vpole
      real(dp):: t

c Statement functions
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      t(p1,p2,p3)=s(p1,p2)+s(p1,p3)+s(p2,p3)

      prp34=cone/cplx2(s(p3,p4)-wmass**2,wmass*wwidth)

      amp =
     & (za(p2,p3)**2*zb(p4,p3))/(za(p1,p5)*za(p2,p6)*za(p5,p6))*Lsm1(-s(p1,p5),-t(p1,p5,p6),-s(p1,p6),-t(p1,p5,p6))
     & +(za(p1,p2)*za(p2,p3)**2*zb(p3,p4))/(za(p1,p5)*za(p1,p6)*za(p2,p5)*za(p2,p6))
     & *(Lsm1_2me(t(p1,p2,p6),t(p1,p5,p6),s(p1,p6),s(p3,p4))
     & +Lsm1(-s(p1,p6),-t(p1,p2,p6),-s(p2,p6),-t(p1,p2,p6)))
     & +(za(p1,p2)*zb(p4,p5))/(2*za(p1,p5)*za(p1,p6)*za(p2,p6))*
     & (3*za(p2,p3)*L0(-t(p1,p2,p6),-s(p3,p4))
     & +za(p3,p5)*zab2(p2,p1,p6,p5)*(L1(-t(p1,p2,p6),-s(p3,p4)))/(s(p3,p4)))
     & +1._dp/2._dp*(za(p1,p2)*zb(p1,p6))/(za(p1,p5)*za(p1,p6)*za(p2,p5))*(za(p2,p3)*(3*za(p1,p3)*zb(p3,p4)-2*za(p1,p5)*zb(p4,p5))
     &  -za(p1,p3)*zab2(p2,p1,p5,p4))*(L0(-t(p1,p2,p6),-s(p2,p6)))/(s(p2,p6))
     & -(za(p1,p2)**2*za(p1,p3)*zb(p1,p6)*zb(p1,p4))/(2*za(p1,p5)*za(p1,p6)*za(p2,p5))*(L1(-t(p1,p2,p6),-s(p2,p6)))/(s(p2,p6))
     & +(za(p1,p2)**2*za(p2,p3)*zb(p1,p4))/(2*za(p1,p5)*za(p1,p6)*za(p2,p5)*za(p2,p6))
     & +(za(p1,p2)*za(p2,p3)*zb(p4,p5))/(za(p1,p5)*za(p1,p6)*za(p2,p6))
     & +(za(p2,p3)*zb(p5,p6)*zab2(p1,p2,p3,p4))/(2*za(p1,p5)*za(p1,p6)*t(p2,p3,p4))

      Alo=za(p1,p2)*za(p2,p3)**2*zb(p4,p3)
     & /(za(p1,p5)*za(p1,p6)*za(p2,p5)*za(p2,p6))

      Vpole=
     &  -1.5_dp*(lnrat(musq,-t(p1,p2,p6))+2._dp)
     &  -0.5_dp*(lnrat(musq,-s(p1,p6))**2+lnrat(musq,-s(p2,p6))**2)
     &  -2._dp*epinv*epinv2
     &  -epinv*(lnrat(musq,-s(p1,p6))+lnrat(musq,-s(p2,p6)))
     &  -1.5_dp*epinv

      amp=amp+Alo*Vpole

      amp=amp*prp34

      return
      end

