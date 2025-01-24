!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c----- fills amplitudes for gg=>HZ process
c===== these are the pieces which come from triangle topologies
      subroutine gg_HZ_tri(p,amp,mt2)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'cplx.h'
      include 'scalarselect.f'
      real(dp):: p(mxpart,4),mt2
      complex(dp):: amp(2,2,2)
      complex(dp):: ggHZ_pp_tri
      external ggHZ_pp_tri
      complex(dp):: prop_34,prop_12

      call spinoru(4,p,za,zb)
c====== z propagtors
      prop_34=s(3,4)/cplx2(s(3,4)-zmass**2,zmass*zwidth)
      prop_12=s(1,2)/cplx2(s(1,2)-zmass**2,zmass*zwidth)
c---- for checking total xs use the formula below
c      prop_12=s(1,2)/cplx2(s(1,2)-zmass**2,zip)

 !------ left lepton amplitudes
      amp(2,2,1)=ggHZ_pp_tri(1,2,3,4,za,zb,mt2)
      amp(2,1,1)=czip
      amp(1,2,1)=czip
      amp(1,1,1)=-ggHZ_pp_tri(1,2,4,3,zb,za,mt2)

c------ right lepton amplitudes
      amp(2,2,2)=ggHZ_pp_tri(1,2,4,3,za,zb,mt2)
      amp(2,1,2)=czip
      amp(1,2,2)=czip
      amp(1,1,2)=-ggHZ_pp_tri(1,2,3,4,zb,za,mt2)

      amp(:,:,:)=prop_12*prop_34*amp(:,:,:)
      return
      end


      function ggHZ_pp_tri(i1,i2,i3,i4,za,zb,mt2)
      implicit none
      include 'types.f'
      complex(dp)::ggHZ_pp_tri
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      integer:: i1,i2,i3,i4
      real(dp):: mt2
      integer,parameter :: Nbox=3,Ntri=6
      integer,parameter :: c12=3
      complex(dp):: ci(Ntri),rat,cmzsq

      complex(dp):: D0(Nbox),C0(Ntri)
      common/ggZH_basisint/D0,C0
!$omp threadprivate(/ggZH_basisint/)

c----- note I'm using the complex mass scheme to ensure Gauge invariance,
      cmzsq=zmass**2-im*zmass*zwidth



      ci(:)=czip

c      C0(c12)=loopI3(s12,zip,zip,mt2,mt2,mt2,musq,0)

      ci(c12)=  (2*(-((mt2*za(i1,i3)*zb(i2,i1)*zb(i4,i1))/za(i1,i2)) -
     &      (mt2*za(i2,i3)*zb(i2,i1)*zb(i4,i2))/za(i1,i2)))/
     &  (za(i1,i2)*za(i3,i4)*zb(i2,i1)*zb(i4,i3))

c------ longitudinal pieces
      ci(c12)=ci(c12)+
     &(2*mt2*zb(i2,i1)**2*(za(i1,i3)*zb(i4,i1) + za(i2,i3)*zb(i4,i2)))/
     &  (cmzsq*s(i1,i2)*s(i3,i4))

      rat= (-2*(-((za(i1,i3)*zb(i2,i1)*zb(i4,i1))/za(i1,i2)) -
     &      (za(i2,i3)*zb(i2,i1)*zb(i4,i2))/za(i1,i2)))/
     &  (za(i1,i2)*za(i3,i4)*zb(i2,i1)*zb(i4,i3))

c------ longitudinal pieces

      rat=rat+
     &  (-2*zb(i2,i1)**2*
     &    (za(i1,i3)*zb(i4,i1) + za(i2,i3)*zb(i4,i2)))/
     &  (cmzsq*s(i1,i2)*s(i3,i4))

      ggHZ_pp_tri=C0(c12)*ci(c12)-rat/two

      return
      end
