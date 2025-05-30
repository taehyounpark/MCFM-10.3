!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function A52(j1,j2,j3,j4,j5,za,zb)
c  Amplitudes taken from Appendix IV of
c  Z.~Bern, L.~J.~Dixon and D.~A.~Kosower,
c  %``One loop amplitudes for e+ e- to four partons,''
c  Nucl.\ Phys.\ B {\bf 513}, 3 (1998)
c  [hep-ph/9708239].
c  Modified to remove momentum conservation relations
      implicit none
      include 'types.f'
      complex(dp):: A52
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      include 'epinv2.f'
      integer:: j1,j2,j3,j4,j5
      complex(dp):: Vcc,Fcc,Vsc,Fsc,l12,l45,L0,L1,Lsm1,A5lom
      complex(dp):: lnrat,zab2
      real(dp):: s123

      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      s123=s(j1,j2)+s(j2,j3)+s(j3,j1)

      l12=lnrat(musq,-s(j1,j2))
      l45=lnrat(musq,-s(j4,j5))
c    -i * A5tree
c     Eq(IV.6)
      A5lom=za(j2,j4)*zab2(j2,j1,j3,j5)/(za(j2,j3)*za(j3,j1)*s123)
c     Eq(IV.7)
      Vcc=-(epinv*epinv2+epinv*l12+half*l12**2)
     & -two*(epinv+l45)-four
c     Eq(IV.8)
      Fcc=-za(j2,j4)*zab2(j2,j1,j3,j5)/(za(j2,j3)*za(j3,j1)*s123)
     & *Lsm1(-s(j1,j2),-s123,-s(j1,j3), -s123)
     & +zab2(j2,j1,j3,j5)*(za(j1,j2)*za(j3,j4)-za(j1,j4)*za(j2,j3))
     &  /(za(j2,j3)*za(j1,j3)**2*s123)
     & *Lsm1(-s(j1,j2),-s123,-s(j2,j3),-s123)
     & +two*zb(j1,j3)*za(j1,j4)*zab2(j2,j1,j3,j5)/(za(j1,j3)*s123)
     & *L0(-s(j2,j3),-s123)/s123

c--subleading N

c     Eq(IV.9)
      Vsc=half*(epinv+l45)+half
c     Eq(IV.10)
      Fsc=za(j1,j4)*za(j2,j3)*zab2(j1,j2,j3,j5)/(za(j1,j3)**3*s123)
     & *Lsm1(-s(j1,j2),-s123,-s(j2,j3),-s123)
     & +half*za(j4,j1)*zb(j1,j3)**2*za(j2,j3)*zab2(j1,j2,j3,j5)
     & /(za(j1,j3)*s123)*L1(-s123,-s(j2,j3))/s(j2,j3)**2
     & +za(j1,j4)*za(j2,j3)*zb(j3,j1)*zab2(j1,j2,j3,j5)
     & /(za(j1,j3)**2*s123)*L0(-s123,-s(j2,j3))/s(j2,j3)
     & -za(j2,j1)*zb(j1,j3)*za(j4,j3)*zb(j3,j5)/za(j1,j3)
     & *L1(-s123,-s(j1,j2))/s(j1,j2)**2
     & -za(j2,j1)*zb(j1,j3)*za(j3,j4)*zab2(j1,j2,j3,j5)
     & /(za(j1,j3)**2*s123)*L0(-s123,-s(j1,j2))/s(j1,j2)
     & +half*zab2(j4,j1,j2,j3)*(zb(j1,j3)*zb(j2,j5)+zb(j2,j3)*zb(j1,j5))
     & /(zb(j1,j2)*zb(j2,j3)*za(j1,j3)*s123)

      A52=(Vcc+Vsc)*A5lom+Fcc+Fsc

      return
      end


