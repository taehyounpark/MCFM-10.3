!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function a6virtLL(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a6virtLL
c***************************************************
c virtual amplitude for
c 0->q(p1)+qb(p2)+gam(p3)+gam(p4)+lb*p5)+l(p6)
c where both photon are coming from lepton line
c***************************************************
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer st
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: a6treeg,vll
c-----
      a6virtLL = a6treeg(st,j1,j2,j3,j4,j5,j6,za,zb)
     &              *vll(st,j5,j6,j3,j4,j1,j2)
c-----
      return
      end


      function vll(st,j1,j2,j3,j4,j5,j6)
      implicit none
      include 'types.f'
      complex(dp):: vll
c-----divergent part
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f'
      include 'epinv2.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: Lnrat,xl12
      integer st
c-----
      xl12=Lnrat(musq,-s(j1,j2))
c-----
      vll = -4._dp
     &      -(epinv*epinv2+epinv*xl12+half*xl12**2)
     &      -(3._dp/2._dp)*(epinv+xl12)
c-----
      return
      end
