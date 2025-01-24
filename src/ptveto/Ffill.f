!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine Ffill(j,lnmuonpt,R,F)
!     *R dependence of BNR F function
!     from 1307.0025 Eq.60 and A20 and A21
!     for j=0 returns powers of as/4/pi in gluon exponent 
!     for j=1 returns powers of as/4/pi in quark exponent 
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'zeta.f'
      include 'nfl.f'
      include 'first.f'
      integer::j,k
      real(dp)::lnmuonpt,R,F(3),LBNR,d2veto,d3veto
      real(dp),save:: beta0,beta1,beta2,beta3,beta4,dRA,DRf,CAA,CAF,CFA,CFF,
     & GammaA0(0:1),GammaA1(0:1),GammaA2(0:1),GammaA3(0:1),GammaA4(0:1),CR(0:1)
      real(dp) :: gammagq0(0:1),gammagq1(0:1),gammagq2(0:1),gammagq3(0:1)
      real(dp) :: gamma0S,gamma1S,gamma2S,gamma3S
!$omp threadprivate(GammaA0,GammaA1,GammaA2,GammaA3,GammaA4,CR,beta0,beta1,beta2,beta3,beta4)

      if (first) then
!--- !     Taken from 9701390, Eq11
!---         beta0=11._dp-2._dp*nfl/3._dp
!---         beta1=102._dp-38*nfl/3._dp
!---         beta2=2857/2._dp-5033/18._dp*nfl+325/54._dp*nfl**2
!---         beta3=149753._dp/6._dp+3564*zeta3
!---      &   -(1078361/162._dp+6508/27._dp*zeta3)*nfl
!---      &   +(50065/162._dp+6472/81._dp*zeta3)*nfl**2
!---      &   +1093/729._dp*nfl**3
!---         CR(0)=CA
!---         CR(1)=CF
!---
!---         do k=0,1
!--- !        GammaA0(k)=4*CR(k)
!--- !        GammaA1(k)=16*CR(k)
!--- !     &   *(67/36d0*CA-5*TF/9*nfl-CA*pisq/12d0)
!--- !        GammaA2(k)=64*CR(k)
!--- !     &  *(CA**2*(11/24d0*zeta3+245/96d0-67/216d0*pisq+11/720d0*pisq**2)
!--- !     &   +CF*TF*nfl*(zeta3-55/48d0)
!--- !     &   +CA*TF*nfl*(-7*zeta3/6d0-209/216d0+5*pisq/54d0)
!--- !     &  -(TF*nfl)**2/27d0)
!--- !         GammaA3(k)=256*CR(k)
!--- !     &   *(CA**3*(1309/432d0*zeta3-11/144d0*pisq*zeta3
!--- !     &   -zeta3**2/16d0-451*zeta5/288d0+42139/10368d0
!--- !     &   -5525*pisq/7776d0+451*pisq**2/5760d0-313*pisq**3/90720d0)
!--- !     &   +CA**2*TF*nfl*(-361*zeta3/54d0+7*pisq*zeta3/36d0+131*zeta5/72d0
!--- !     &   -24137/10368d0+635*pisq/1944d0-11*pisq**2/2160d0)
!--- !     &   +CA*CF*TF*nfl*(29/9d0*zeta3-pisq*zeta3/6d0+5/4d0*zeta5
!--- !     &   -17033/5184d0+55/288d0*pisq-11/720d0*pisq**2)
!--- !     &   +CF**2*TF*nfl*(37/24d0*zeta3-5/2d0*zeta5+143/288d0)
!--- !     &   +CA*(TF*nfl)**2*(35/27d0*zeta3-7/1080d0*pisq**2-19/972d0*pisq
!--- !     &   +923/5184d0)
!--- !     &  +CF*(TF*nfl)**2*(-10/9d0*zeta3+1/180d0*pisq**2+299/648d0)
!--- !     &  +(TF*nfl)**3*(-1/81d0+2*zeta3/27d0))
!---
!---         GammaA0(k)=4*CR(k)
!---         GammaA1(k)=16*CR(k)*(67/12d0-pisq/4d0-nfl*5/18d0)
!---         GammaA2(k)=64*CR(k)
!---      &   *((735/32d0+33/8d0*zeta3-67/24d0*pisq+11/80d0*pisq**2)
!---      &   +nfl*(-319/144d0-13/12d0*zeta3+5/36d0*pisq)
!---      &   -nfl**2/108d0)
!---        GammaA3(k)=256*CR(k)
!---      &   *(42139/384d0+1309/16d0*zeta3-1353/32d0*zeta5-27/16d0*zeta3**2
!---      &   -5525/288d0*pisq-33/16d0*pisq*zeta3+1353/640d0*pisq**2
!---      &   -313/3360d0*pisq**3
!---      &   +nfl*(-344345/20736d0+1219/144d0*zeta5-2405/108d0*zeta3
!---      &   +50/27d0*pisq+13/24d0*pisq*zeta3-77/1440d0*pisq**2)
!---      &   +nfl**2*(+17875/62208d0+65/108d0*zeta3-19/1296d0*pisq
!---      &   -13/4320d0*pisq**2)
!---      &   +nfl**3*(-1/648d0+zeta3/108d0))
!---
!---         enddo
!---
!--- !     Taken from 1911.10174, Eq 6.5
!---         CAA=135/8d0        !CAA=N^2*(N^2+36)/24
!---         CAF=15/16d0        !CAF=N*(N^2+6)/48;
!---         CFA=15/16d0*8d0/3d0        !CFA=N*(N^2+6)/48;
!---         CFF=5/36d0         !CFF=(N^4-6*N^2+18)/96/N^2;
!---         dRA=55/12._dp*zeta5+1/6._dp*zeta3-3/2._dp*zeta3**2
!---      &  -pisq/12._dp-31/7560._dp*pisq**3
!---         dRf=pisq/6._dp-zeta3/3._dp-5/3._dp*zeta5
!--- !     Taken from 1911.10174, Eq 6.4
!---         GammaA3(0)=GammaA3(0)+256*(CAA*dRA+nfl*CAF*dRF)
!---         GammaA3(1)=GammaA3(1)+256*(CFA*dRA+nfl*CFF*dRF)
!---
!---         write(6,*) 'old'
!---         write(6,*) 'CR(0)',CR(0)
!---         write(6,*) 'CR(1)',CR(1)
!---         write(6,*) 'beta0',beta0
!---         write(6,*) 'beta1',beta1
!---         write(6,*) 'beta2',beta2
!---         write(6,*) 'beta3',beta3
!---         do k=0,1
!---         write(6,*) 'k,GammaA0',k,GammaA0(k)
!---         write(6,*) 'k,GammaA1',k,GammaA1(k)
!---         write(6,*) 'k,GammaA2',k,GammaA2(k)
!---         write(6,*) 'k,GammaA3',k,GammaA3(k)
!---         enddo

         CR(0)=CA
         CR(1)=CF
         call Anomfill(nfl,
     &    beta0,beta1,beta2,beta3,beta4,
     &    GammaA0,GammaA1,GammaA2,GammaA3,GammaA4,
     &    gammagq0,gammagq1,gammagq2,gammagq3,
     &    gamma0S,gamma1S,gamma2S,gamma3S)

!        write(6,*) 'new'
!        write(6,*) 'CR(0)',CR(0)
!        write(6,*) 'CR(1)',CR(1)
!        write(6,*) 'beta0',beta0
!        write(6,*) 'beta1',beta1
!        write(6,*) 'beta2',beta2
!        write(6,*) 'beta3',beta3
!        do k=0,1
!        write(6,*) 'k,GammaA0',k,GammaA0(k)
!        write(6,*) 'k,GammaA1',k,GammaA1(k)
!        write(6,*) 'k,GammaA2',k,GammaA2(k)
!        write(6,*) 'k,GammaA3',k,GammaA3(k)
!        enddo
!        pause

         first=.false.

      endif

      LBNR=2*lnmuonpt
! d1veto=0 (MS-bar scheme)
      F(1)=GammaA0(j)*LBNR
      F(2)=GammaA0(j)*beta0*LBNR**2/2d0
     & +GammaA1(j)*LBNR+d2veto(R,CR(j))
      F(3)=GammaA0(j)*beta0**2*LBNR**3/3d0
     & +(GammaA0(j)*beta1+2*GammaA1(j)*beta0)*LBNR**2/2d0
     & +LBNR*(GammaA2(j)+2*beta0*d2veto(R,CR(j)))+d3veto(R,CR(j))
      return
      end

      function d2veto(R,CB)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'zeta.f'
      include 'nfl.f'
      real(dp)::d2veto,R,d2b,CB,fBNR
      d2b=CB*((808/27d0-28*zeta3)*CA-224/27d0*TF*nfl)
      d2veto=d2b-32*CB*fBNR(R,CB)
      return
      end

      function d2veto0(R,CB)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'zeta.f'
      include 'nfl.f'
      real(dp)::d2veto0,R,d2b,CB,fBNR0
      d2b=CB*((808/27d0-28*zeta3)*CA-224/27d0*TF*nfl)
      d2veto0=d2b-32*CB*fBNR0(R,CB)
      return
      end

      function d2veto2(R,CB)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'zeta.f'
      include 'nfl.f'
      real(dp)::d2veto2,R,d2b,CB,fBNR2
      d2b=CB*((808/27d0-28*zeta3)*CA-224/27d0*TF*nfl)
      d2veto2=d2b-32*CB*fBNR2(R,CB)
      return
      end

      function d2veto4(R,CB)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'zeta.f'
      include 'nfl.f'
      real(dp)::d2veto4,R,d2b,CB,fBNR4
      d2b=CB*((808/27d0-28*zeta3)*CA-224/27d0*TF*nfl)
      d2veto4=d2b-32*CB*fBNR4(R,CB)
      return
      end

      function d3veto(R,CB)
          use ptveto
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nfl.f'
      include 'first.f'
      include 'taucut.f'
      include 'mpicommon.f'
      real(dp)::d3veto,R,CB,lnR

      if (useBanfid3veto) then
! leading log value taken from 1511.02886, Eq.(C2)
! (reasonable variation of R0 is between 0.5 and 2)
        lnR=log(R/d3vetoR0)
        d3veto=-64*CB*lnR**2*(1.803136_dp*CA**2-0.589237_dp*CA*nfl
     &                      +0.36982_dp*CF*nfl-0.05893_dp*nfl**2)
        if (first) then
!$omp master
          if (rank == 0) then
          write(6,*) 'Using Banfi et al. approximation for d3veto, R0 = ',d3vetoR0
          endif
!$omp end master
          first = .false.
        endif
      else
! Estimate from Eq. (66) of 1307.0025, with additional overall minus
! sign that I believe should have been present to reproduce results
! therein (reasonable variation of d3vetokappa is between -4 and 4)
        d3veto=-d3vetokappa*64*CB*CA**2*log(2._dp/R)**2
        if (first) then
!$omp master
          if (rank == 0) then
          write(6,*) 'Using BNR estimate of d3veto, kappa = ',d3vetokappa
          endif
!$omp end master
          first = .false.
        endif
      endif

      return
      end
      
