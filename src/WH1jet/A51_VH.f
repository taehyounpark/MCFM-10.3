!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function A51_VH(j1,j2,j3,j4,j5,za,zb)

c  Amplitudes taken from Appendix IV of
c  Z.~Bern, L.~J.~Dixon and D.~A.~Kosower,
c  %``One loop amplitudes for e+ e- to four partons,''
c  Nucl.\ Phys.\ B {\bf 513}, 3 (1998)
c  [hep-ph/9708239].
c  Modified to remove momentum conservation relations
c===== C.Williams added missing terms which go like <3|P_125|4] in
c===== W + 1 jet and vanish for that case but are non zero for us

      implicit none
      include 'types.f'
      complex(dp):: A51_VH
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      include 'epinv2.f'
      integer:: j1,j2,j3,j4,j5
      complex(dp):: Vcc,Fcc,Vsc,Fsc,l12,l23,L0,L1,Lsm1,A5lom
      complex(dp):: lnrat,zab2,zab3
      real(dp):: s123
      complex(dp):: MPC
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab3(j1,j2,j3,j4,j5)=za(j1,j2)*zb(j2,j5)+za(j1,j3)*zb(j3,j5)
     &     +za(j1,j4)*zb(j4,j5)

c    -i * A5tree
c==== debug
c      epinv=0_dp
c      musq=1_dp

      s123=s(j1,j2)+s(j2,j3)+s(j3,j1)
      A5lom=-za(j3,j4)*zab2(j3,j1,j2,j5)/(za(j1,j2)*za(j2,j3)*s123)
      l12=lnrat(musq,-s(j1,j2))
      l23=lnrat(musq,-s(j2,j3))


c--leading N
      Vcc=
     & -(epinv*epinv2+epinv*l12+half*l12**2)
     & -(epinv*epinv2+epinv*l23+half*l23**2)
     & -two*(epinv+l23)-four

      Fcc=zab2(j3,j1,j2,j5)/(za(j1,j2)*za(j2,j3)*s123)
     & *(za(j3,j4)*Lsm1(-s(j1,j2),-s123,-s(j2,j3),-s123)
     & +two*za(j3,j1)*zab2(j4,j2,j3,j1)
     &   *L0(-s(j2,j3),-s123)/s123)

      Vsc =half*(epinv+l23)+one
      Fsc =-za(j3,j4)*za(j3,j1)*zb(j1,j5)
     & /(za(j1,j2)*za(j2,j3))*L0(-s(j2,j3),-s123)/s123
     & +half*za(j3,j1)**2*zb(j1,j5)*zab2(j4,j2,j3,j1)
     & /(za(j1,j2)*za(j2,j3))*L1(-s(j2,j3),-s123)/s123**2

c===== correction for piece
      MPC=-(3._dp*L0(-s(j2,j3),-s123)*za(j1,j3)*
     &   zab3(j4,j1,j2,j3,j5)*zb(j2,j1))/(2.*s123**2*za(j1,j2))

      A51_VH=(Vcc+Vsc)*A5Lom+Fcc+Fsc+MPC

c      write(6,*) 'LO LC ',A5LOm
c      write(6,*) 'NLO, epinv (2) ',-2_dp*A5LOm
c      write(6,*) 'NLO, epinv(1) ', (-l12-l23-2_dp+half)*A5LOm
c      write(6,*) 'NLO   ',A51_VH

c      write(6,*) 'Box coefficient ',zab2(j3,j1,j2,j5)/(za(j1,j2)*za(j2,j3)*s123!)
c     & *(za(j3,j4))
c      write(6,*) 'rational coefficient', +half*za(j3,j1)**2*zb(j1,j5)*zab2(j4,j2,j3,j1)
c     & /(za(j1,j2)*za(j2,j3))/s123**2*(cone/cplx1(one-s(j2,j3)/s123))
c      pause


      return
      end
