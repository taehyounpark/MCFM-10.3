!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function a6vQLslc(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a6vQLslc
c***************************************************
c virtual amplitude for
c 0->q(p1)+qb(p2)+glu(p3)+gam(p4)+lb*p5)+l(p6)
c where one photon is coming from quark line
c and another photon coming from lepton line
c***************************************************
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: st
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: a6treeQLslc,vQLslc,fQLslc
c-----
      a6vQLslc = a6treeQLslc(st,j1,j2,j3,j4,j5,j6,za,zb)
     &           *vQLslc(st,j1,j2,j3,j4,j5,j6)
c-----
      a6vQLslc = a6vQLslc+fQLslc(st,j1,j2,j3,j4,j5,j6,za,zb)
c-----
      return
      end


      function vQLslc(st,j1,j2,j3,j4,j5,j6)
      implicit none
      include 'types.f'
      complex(dp):: vQLslc
c-----divergent part
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f'
      include 'epinv2.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: Lnrat,xl12,xl456
      integer:: st
      real(dp):: t
c-----
      xl12=Lnrat(musq,-s(j1,j2))
      xl456=Lnrat(musq,-t(j4,j5,j6))
c-----
      vQLslc = -(7._dp/2._dp)
     &         -(epinv*epinv2+epinv*xl12+half*xl12**2)
     &         -(3._dp/2._dp)*(epinv+xl456)
c-----
      return
      end


      function fQLslc(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: fQLslc
c-----finite part
      integer:: st
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: Lsm1,L0,L1
      real(dp):: t

c---- helicity stamps
c     'q+qb-g+g+lb-l+'=4
c     'q+qb-g+g-lb-l+'=5

      if(st==4) then
      fQLslc =
     &      + za(j2,j5)**2/(za(j1,j3)*za(j2,j3)*za(j4,j5)*za(j4,j6))
     &        *Lsm1(-s(j1,j2),-t(j4,j5,j6),-s(j1,j3),-t(j4,j5,j6))
     &      + za(j1,j2)**2*za(j3,j5)**2
     &        /(za(j1,j3)**3*za(j2,j3)*za(j4,j5)*za(j4,j6))
     &        *Lsm1(-s(j1,j2),-t(j4,j5,j6),-s(j2,j3),-t(j4,j5,j6))
     &      - ( 2._dp*s(j1,j3)*za(j1,j5)*za(j2,j5)
     &         -za(j2,j3)*zb(j3,j1)*za(j1,j5)**2 )
     &        /(za(j1,j3)**2*za(j4,j5)*za(j4,j6))
     &        *L0(-t(j4,j5,j6),-s(j2,j3))/s(j2,j3)
     &      - zb(j1,j3)**2*za(j2,j3)*za(j1,j5)**2
     &        /(2._dp*za(j1,j3)*za(j4,j5)*za(j4,j6))
     &        *L1(-t(j4,j5,j6),-s(j2,j3))/s(j2,j3)**2
     &      - za(j2,j1)*zb(j1,j3)*za(j1,j5)*za(j3,j5)
     &        /(za(j1,j3)**2*za(j4,j5)*za(j4,j6))
     &        *L0(-t(j4,j5,j6),-s(j1,j2))/s(j1,j2)
     &      + za(j2,j1)*zb(j1,j3)*za(j5,j3)
     &        *( za(j5,j4)*zb(j4,j3) + za(j5,j6)*zb(j6,j3) )
     &        /(za(j1,j3)*za(j4,j5)*za(j4,j6))
     &        *L1(-t(j4,j5,j6),-s(j1,j2))/s(j1,j2)**2
     &      -  ( za(j5,j1)*zb(j1,j3)+za(j5,j2)*zb(j2,j3) )
     &        *( zb(j1,j3)*( za(j5,j4)*zb(j4,j2) + za(j5,j6)*zb(j6,j2) )
     &          +zb(j2,j3)*( za(j5,j4)*zb(j4,j1) + za(j5,j6)*zb(j6,j1) )
     &         )/(2._dp*t(j4,j5,j6)
     &            *zb(j1,j2)*zb(j2,j3)*za(j1,j3)*za(j4,j5)*za(j4,j6))
      elseif(st==5) then
      fQLslc =
     &      -  (za(j2,j1)*zb(j1,j6)+ za(j2,j3)*zb(j3,j6))
     &        *(za(j2,j4)*zb(j4,j6)+ za(j2,j5)*zb(j5,j6))
     &        /(t(j4,j5,j6)*za(j1,j3)*za(j2,j3)*zb(j4,j5)*zb(j4,j6))
     &        *Lsm1(-s(j1,j2),-t(j4,j5,j6),-s(j1,j3),-t(j4,j5,j6))
     &      - za(j1,j2)**2*(za(j3,j1)*zb(j1,j6)+ za(j3,j2)*zb(j2,j6))
     &        *(za(j3,j4)*zb(j4,j6)+ za(j3,j5)*zb(j5,j6))
     &        /(t(j4,j5,j6)*za(j1,j3)**3*za(j2,j3)*zb(j4,j5)*zb(j4,j6))
     &        *Lsm1(-s(j1,j2),-t(j4,j5,j6),-s(j2,j3),-t(j4,j5,j6))
     &      + (za(j1,j2)*zb(j2,j6)+ za(j1,j3)*zb(j3,j6))
     &        *( 2._dp*s(j1,j3)*(za(j2,j4)*zb(j4,j6)+ za(j2,j5)*zb(j5,j6))
     &          -za(j2,j3)*zb(j3,j1)
     &           *(za(j1,j4)*zb(j4,j6)+ za(j1,j5)*zb(j5,j6)) )
     &        /(t(j4,j5,j6)*za(j1,j3)**2*zb(j4,j5)*zb(j4,j6))
     &        *L0(-t(j4,j5,j6),-s(j2,j3))/s(j2,j3)
     &      + zb(j1,j3)**2*za(j2,j3)
     &        *(za(j1,j2)*zb(j2,j6)+ za(j1,j3)*zb(j3,j6))
     &        *(za(j1,j4)*zb(j4,j6)+ za(j1,j5)*zb(j5,j6))
     &        /(2._dp*t(j4,j5,j6)*za(j1,j3)*zb(j4,j5)*zb(j4,j6))
     &        *L1(-t(j4,j5,j6),-s(j2,j3))/s(j2,j3)**2
     &      + za(j2,j1)*zb(j1,j3)
     &        *(za(j3,j1)*zb(j1,j6)+ za(j3,j2)*zb(j2,j6))
     &        *(za(j1,j4)*zb(j4,j6)+ za(j1,j5)*zb(j5,j6))
     &        /(t(j4,j5,j6)*za(j1,j3)**2*zb(j4,j5)*zb(j4,j6))
     &        *L0(-t(j4,j5,j6),-s(j1,j2))/s(j1,j2)
     &      + za(j2,j1)*zb(j1,j3)*zb(j6,j3)
     &        *( za(j3,j4)*zb(j4,j6) + za(j3,j5)*zb(j5,j6) )
     &        /(za(j1,j3)*zb(j4,j5)*zb(j4,j6))
     &        *L1(-t(j4,j5,j6),-s(j1,j2))/s(j1,j2)**2
     &      - zb(j3,j6)*( zb(j1,j3)*zb(j6,j2) + zb(j2,j3)*zb(j6,j1) )
     &        /(2._dp*zb(j1,j2)*zb(j2,j3)*za(j1,j3)*zb(j4,j5)*zb(j4,j6))
      else
      write(6,*) 'unimplemented st'
      stop
      endif
c-----
      return
      end

