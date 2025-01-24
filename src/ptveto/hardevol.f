!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function hardevol(ii,order,qsq,mu,muh)
      use evolution
      use LHAPDF, only: getalphas
      use ptveto, only: gghsinglestep, timelikemusq
! hard evolution from Appendix B of 1109.6027 (Becher, Neubert, Wilhelm)
! i=0 (gg), i=1 (qqbar)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nfl.f'
      include 'zeta.f'
      include 'first.f'
      integer ii,order,k
      real(dp) :: hardevol,qsq,mu,muh
      real(dp) :: aGamma,aGamma0,aGamma1,aGamma2,aGamma3
      real(dp) :: agammai,agammai0,agammai1,agammai2,agammai3
      real(dp) :: S,S0,S1,S2,S3
      real(dp) :: alphasMu,alphasMuH,r
      real(dp),save :: beta0,beta1,beta2,beta3,beta4
      real(dp),save :: CR(0:1),GammaA0(0:1),GammaA1(0:1),GammaA2(0:1),GammaA3(0:1),GammaA4(0:1)
      real(dp),save :: gammagq0(0:1),gammagq1(0:1),gammagq2(0:1),gammagq3(0:1)
      real(dp),save :: gamma0S,gamma1S,gamma2S,gamma3S
      real(dp) :: Gamma0,Gamma1,Gamma2,Gamma3,Gamma4
      real(dp) :: gammai0,gammai1,gammai2,gammai3
      complex(dp) :: alphasMu_c,alphasMuH_c,S_c,aGamma_c,agammai_c
      complex(dp) :: getason2pic
!$omp threadprivate(beta0,beta1,beta2,beta3,beta4)
!$omp threadprivate(CR,GammaA0,GammaA1,GammaA2,GammaA3,GammaA4)
!$omp threadprivate(gammagq0,gammagq1,gammagq2,gammagq3)
!$omp threadprivate(gamma0S,gamma1S,gamma2S,gamma3S)

      if (first) then
         call Anomfill(nfl,
     &   beta0,beta1,beta2,beta3,beta4,
     &   GammaA0,GammaA1,GammaA2,GammaA3,GammaA4,
     &   gammagq0,gammagq1,gammagq2,gammagq3,
     &   gamma0S,gamma1S,gamma2S,gamma3S)
         if (gghsinglestep .eqv. .false.) then
! In this implementation use gammaS instead of gammag, which comes
! with an evolution factor 4*S-2*aGamma*log-2*agammaS rather than
! 4*S-2*aGamma*log-4*agammai that is coded below, hence factor of 1/2
! (see for example Eq. (A3) of arXiv:1205.3806)
           gammagq0(0)=gamma0S/2._dp
           gammagq1(0)=gamma1S/2._dp
           gammagq2(0)=gamma2S/2._dp
           gammagq3(0)=gamma3S/2._dp
         endif
         first=.false.
      endif

!     Fill expansion coefficients with flavor-specific values
      Gamma0=GammaA0(ii)
      Gamma1=GammaA1(ii)
      Gamma2=GammaA2(ii)
      Gamma3=GammaA3(ii)
      Gamma4=GammaA4(ii)
      gammai0=gammagq0(ii)
      gammai1=gammagq1(ii)
      gammai2=gammagq2(ii)
      gammai3=gammagq3(ii)

! hard evolution from Eq. (B1) in 1109.6027
! except we note that it extends to complex parameters via the use of absolute value
! and the power should be (-qsq/muh**2), c.f. Eq. (22) of 0809.4283
      alphasMu=getalphas(mu)
      alphasMuH=getalphas(muh)

      if (timelikemusq) then
!------ Section for treating complex alphas (musq < 0)

! factorization scale is space-like
        alphasMu_c=cmplx(alphasMu,0._dp,kind=dp)
! matching scale is time-like, translation to complex alphas from 0809.4283, Eq. (25)
        alphasMuH_c=getason2pic(muh**2)*twopi

        call evolfill(alphasMu_c,alphasMuH_c,order,
     &   Gamma0,Gamma1,Gamma2,Gamma3,Gamma4,
     &   gammai0,gammai1,gammai2,gammai3,
     &   beta0,beta1,beta2,beta3,beta4,
     &   S_c,aGamma_c,agammai_c)

! Note that (-qsq/muh**2) -> (qsq/muh**2) since matching scale is time-like
        hardevol=abs(exp(4*S_c-4*agammai_c)*(qsq/muh**2)**(-2*aGamma_c))

      else
!------ Section for treating real alphas (musq > 0)

        call evolfill(alphasMu,alphasMuH,order,
     &   Gamma0,Gamma1,Gamma2,Gamma3,Gamma4,
     &   gammai0,gammai1,gammai2,gammai3,
     &   beta0,beta1,beta2,beta3,beta4,
     &   S,aGamma,agammai)

        hardevol=exp(4*S-4*agammai)*(qsq/muh**2)**(-2*aGamma)
! equivalent form
!        hardevol=abs(exp(4*S-4*agammai-2*aGamma*log(cmplx(-qsq/muh**2,0._dp,dp))))

      endif

      return
      end
      
