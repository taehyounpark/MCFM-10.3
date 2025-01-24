!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine hfill(j,lnmuonpt,h)
!     Generalization from 1307.0025 Eq 10
!     for j=0 returns powers of as/4/pi of h for gluons
!     for j=0 returns powers of as/4/pi of h for quarks
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'zeta.f'
      include 'nfl.f'
      include 'first.f'
      integer::j,k
      real(dp)::lnmuonpt,h(3),LBNR
      real(dp),save:: GammaA0(0:1),GammaA1(0:1),GammaA2(0:1),CR(0:1),
     & beta0,beta1,gamma0(0:1),gamma1(0:1)
!$omp threadprivate(GammaA0,GammaA1,GammaA2,CR,beta0,beta1,gamma0,gamma1)

      if (first) then
        CR(0)=CA
        CR(1)=CF
        beta0=(11*CA-4*TF*nfl)/3d0
        beta1=34/3d0*CA**2-20/3d0*CA*TF*nfl-4*CF*TF*nfl
* 0903.1126 Eq(A5)
        gamma0(0)=-beta0;
        gamma1(0)=CA**2*(-692/27._dp+11*pisq/18._dp+2*zeta3)
     &   +CA*TF*nfl*(256/27._dp-2*pisq/9._dp)+4*CF*TF*nfl
        gamma0(1)=-3*CF
        gamma1(1)=CF**2*(-3/2._dp+2*pisq-24*zeta3)
     &   +CF*CA*(-961/54._dp-11*pisq/6._dp+26*zeta3)
     &   +CF*TF*nfl*(130/27._dp+2*pisq/3_dp);
        do k=0,1
        GammaA0(k)=4*CR(k)
        GammaA1(k)=16*CR(k)*(67/12d0-pisq/4d0-nfl*5/18d0)
        GammaA2(k)=64*CR(k)
     &   *((735/32d0+33/8d0*zeta3-67/24d0*pisq+11/80d0*pisq**2)
     &   +nfl*(-319/144d0-13/12d0*zeta3+5/36d0*pisq)
     &   -nfl**2/108d0)
        enddo
        first=.false.
      endif

      LBNR=2._dp*lnmuonpt
      
      h(1)=GammaA0(j)*LBNR**2/4._dp-gamma0(j)*LBNR
      h(2)=GammaA0(j)*beta0*LBNR**3/12._dp
     & +(GammaA1(j)-2*gamma0(j)*beta0)*LBNR**2/4._dp-gamma1(j)*LBNR
      h(3)=GammaA0(j)*beta0**2/24._dp*LBNR**4
     & +(1/12._dp*GammaA0(j)*beta1+1/6._dp*GammaA1(j)*beta0
     & -1/3._dp*gamma0(j)*beta0**2)*LBNR**3
     & +(-0.5*gamma0(j)*beta1-gamma1(j)*beta0)*LBNR**2
      return
      end

