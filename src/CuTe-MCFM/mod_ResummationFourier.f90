
!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module qtResummationFourier
      use types
      use constants
      use iso_fortran_env
      implicit none

      public :: fourierM

      private

      real(dp), save :: nf
!$omp threadprivate(nf)

      real(dp), save :: beta0,beta1,beta2,beta3
!$omp threadprivate(beta0,beta1,beta2,beta3)

      real(dp), save :: Gamma0, Gamma1, Gamma2, Gamma3
!$omp threadprivate(Gamma0, Gamma1,Gamma2,Gamma3)

      real(dp), save :: gammaq0,gammaq1, gammaq2, gammag0, gammag1, gammag2
!$omp threadprivate(gammaq0,gammaq1,gammaq2,gammag0,gammag1,gammag2)

      real(dp), save :: d2,d3,d4adj,d4fund
!$omp threadprivate(d2,d3,d4adj,d4fund)

      integer, save :: order, nn
!$omp threadprivate(order,nn)

      logical, save :: initQuark
!$omp threadprivate(initQuark)

      real(dp), save :: qt, q2, mu, alphasMu
!$omp threadprivate(qt,q2,mu,alphasMu)

    contains

    ! wrapper to amos zbesk, see src/Lib/amos/
    function kv(v, z)
        implicit none
        real(dp), intent(in) :: v
        complex(dp), intent(in) :: z
        complex(dp) :: kv

        real(dp) :: cyR, cyI
        integer :: nz, ierr

        call zbesk(real(z), aimag(z), v, 1, 1, cyR, cyI, nz, ierr)
        if (ierr /= 0) then
            !write (*,*) "WARNING: ierr = ", ierr, "in zbesk"
        endif

        kv = cmplx(cyR, cyI, dp)
    end function

    function besselJ0(xx)
        implicit none
        real(dp), intent(in) :: xx
        complex(dp) :: besselJ0

        real(dp) :: jR, jI
        integer :: nz, ierr

        ! second argument, imaginary part set to 0d0
        call zbesj(xx, 0d0, 0d0, 1, 1, jR, jI, nz, ierr)

        besselJ0  = cmplx(jR, jI, dp)
    end function

    function fourierIntegrand_xt(xt)
        use qtResummation_params, only: scalevar_rapidity, &
            scalevar_rapidity_i, scalevar_rapidity_mult
        implicit none
        real(dp), intent(in) :: xt
        real(dp) :: fourierIntegrand_xt

        real(dp) :: aa, etaf
        real(dp) :: Lperp

        real(dp) :: as, Ci, eta
        complex(dp) :: gi
        real(dp) :: gamma1rx, gamma2rx, gamma3rx
        real(dp) :: gammai0, gammai1, gammai2

        real(dp) :: logR

        if (initQuark) then
            Ci = CF
            gammai0 = gammaq0
            gammai1 = gammaq1
            gammai2 = gammaq2
            gamma1rx = -d2*CF
            gamma2rx = -d3*CF
            gamma3rx = -d4fund
        else
            Ci = CA
            gammai0 = gammag0
            gammai1 = gammag1
            gammai2 = gammag2
            gamma1rx = -d2*CA
            gamma2rx = -d3*CA
            gamma3rx = -d4adj
        endif

        aa = Ci*alphasMu/pi
        etaf = aa*log(Q2/mu**2)
        eta = etaf

        if (scalevar_rapidity_i > 0) then
            logR = log(scalevar_rapidity_mult(scalevar_rapidity_i))      
        else
            logR = 0._dp
        endif

        Lperp = log(exp(2._dp*EulerGamma)*xt**2*mu**2/4._dp)

        as = alphasMu/4._dp/pi

        gi = -(eta*Lperp) - as*Ci*Gamma0*logR*Lperp -  &
        (as*beta0*eta*Lperp**2)/2._dp - (as*Ci*Gamma0*Lperp**2)/2._dp

        if (order >= 4) then
            gi = gi + (as*eta*gamma1rX)/(Ci*Gamma0) - (as*eta*Gamma1*Lperp)/Gamma0 - &
                2*as*gammai0*Lperp - as**2*Ci*Gamma1*logR*Lperp - &
                (as**2*beta1*eta*Lperp**2)/2._dp - (as**2*Ci*Gamma1*Lperp**2)/2._dp - &
                (as**2*beta0*eta*Gamma1*Lperp**2)/Gamma0 - &
                as**2*beta0*gammai0*Lperp**2 - &
                (as**2*beta0*Ci*Gamma0*logR*Lperp**2)/2._dp - &
                (as**2*beta0**2*eta*Lperp**3)/3._dp - &
                (as**2*beta0*Ci*Gamma0*Lperp**3)/3._dp - &
                (as**3*beta0**2*Ci*Gamma0*logR*Lperp**3)/3._dp - &
                (as**3*beta0**3*eta*Lperp**4)/4._dp - &
                (as**3*beta0**2*Ci*Gamma0*Lperp**4)/4._dp
        endif

        if (order >= 6) then
            gi = gi + (as**2*eta*gamma2rX)/(Ci*Gamma0) + as**2*gamma1rX*logR + &
                    as**2*gamma1rX*Lperp + &
                    (2*as**2*beta0*eta*gamma1rX*Lperp)/(Ci*Gamma0) - &
                    (as**2*eta*Gamma2*Lperp)/Gamma0 - 2*as**2*gammai1*Lperp + &
                    2*as**3*beta0*gamma1rX*logR*Lperp - &
                    as**3*Ci*Gamma2*logR*Lperp - (as**3*beta2*eta*Lperp**2)/2._dp - &
                    (as**3*beta1*eta*Gamma1*Lperp**2)/Gamma0 + &
                    2*as**3*beta0*gamma1rX*Lperp**2 + &
                    (3*as**3*beta0**2*eta*gamma1rX*Lperp**2)/(Ci*Gamma0) - &
                    (as**3*Ci*Gamma2*Lperp**2)/2._dp - &
                    (3*as**3*beta0*eta*Gamma2*Lperp**2)/(2._dp*Gamma0) - &
                    as**3*beta1*gammai0*Lperp**2 - 2*as**3*beta0*gammai1*Lperp**2 - &
                    (as**3*beta1*Ci*Gamma0*logR*Lperp**2)/2._dp - &
                    as**3*beta0*Ci*Gamma1*logR*Lperp**2 - &
                    (5*as**3*beta0*beta1*eta*Lperp**3)/6._dp - &
                    (as**3*beta1*Ci*Gamma0*Lperp**3)/3._dp - &
                    (2*as**3*beta0*Ci*Gamma1*Lperp**3)/3._dp - &
                    (as**3*beta0**2*eta*Gamma1*Lperp**3)/Gamma0 - &
                    (2*as**3*beta0**2*gammai0*Lperp**3)/3._dp - &
                    (5*as**4*beta0*beta1*Ci*Gamma0*logR*Lperp**3)/6._dp - &
                    as**4*beta0**2*Ci*Gamma1*logR*Lperp**3 - &
                    (13*as**4*beta0**2*beta1*eta*Lperp**4)/12._dp - &
                    (5*as**4*beta0*beta1*Ci*Gamma0*Lperp**4)/8._dp - &
                    (3*as**4*beta0**2*Ci*Gamma1*Lperp**4)/4._dp - &
                    (as**4*beta0**3*eta*Gamma1*Lperp**4)/Gamma0 - &
                    (as**4*beta0**3*gammai0*Lperp**4)/2._dp - &
                    (as**4*beta0**3*Ci*Gamma0*logR*Lperp**4)/4._dp - &
                    (as**4*beta0**4*eta*Lperp**5)/5._dp - &
                    (as**4*beta0**3*Ci*Gamma0*Lperp**5)/5._dp - &
                    (as**5*beta0**4*Ci*Gamma0*logR*Lperp**5)/5._dp - &
                    (as**5*beta0**5*eta*Lperp**6)/6._dp - &
                    (as**5*beta0**4*Ci*Gamma0*Lperp**6)/6._dp
        endif

         if (order >= 8) then
            gi = gi + (as**3*eta*gamma3rX)/(Ci*Gamma0) + as**3*gamma2rX*logR + &
                    (2*as**3*beta1*eta*gamma1rX*Lperp)/(Ci*Gamma0) + &
                    as**3*gamma2rX*Lperp + &
                    (3*as**3*beta0*eta*gamma2rX*Lperp)/(Ci*Gamma0) - &
                    (as**3*eta*Gamma3*Lperp)/Gamma0 - 2*as**3*gammai2*Lperp + &
                    2*as**4*beta1*gamma1rX*logR*Lperp + &
                    3*as**4*beta0*gamma2rX*logR*Lperp - &
                    as**4*Ci*Gamma3*logR*Lperp - (as**4*beta3*eta*Lperp**2)/2._dp - &
                    (as**4*beta2*eta*Gamma1*Lperp**2)/Gamma0 + &
                    2*as**4*beta1*gamma1rX*Lperp**2 + &
                    (7*as**4*beta0*beta1*eta*gamma1rX*Lperp**2)/(Ci*Gamma0) - &
                    (3*as**4*beta1*eta*Gamma2*Lperp**2)/(2._dp*Gamma0) + &
                    3*as**4*beta0*gamma2rX*Lperp**2 + &
                    (6*as**4*beta0**2*eta*gamma2rX*Lperp**2)/(Ci*Gamma0) - &
                    (as**4*Ci*Gamma3*Lperp**2)/2._dp - &
                    (2*as**4*beta0*eta*Gamma3*Lperp**2)/Gamma0 - &
                    as**4*beta2*gammai0*Lperp**2 - 2*as**4*beta1*gammai1*Lperp**2 - &
                    3*as**4*beta0*gammai2*Lperp**2 - &
                    (as**4*beta2*Ci*Gamma0*logR*Lperp**2)/2._dp - &
                    as**4*beta1*Ci*Gamma1*logR*Lperp**2 + &
                    3*as**4*beta0**2*gamma1rX*logR*Lperp**2 - &
                    (3*as**4*beta0*Ci*Gamma2*logR*Lperp**2)/2._dp - &
                    (as**4*beta1**2*eta*Lperp**3)/2._dp - &
                    as**4*beta0*beta2*eta*Lperp**3 - &
                    (as**4*beta2*Ci*Gamma0*Lperp**3)/3._dp - &
                    (2*as**4*beta1*Ci*Gamma1*Lperp**3)/3._dp - &
                    (7*as**4*beta0*beta1*eta*Gamma1*Lperp**3)/(3._dp*Gamma0) + &
                    3*as**4*beta0**2*gamma1rX*Lperp**3 + &
                    (4*as**4*beta0**3*eta*gamma1rX*Lperp**3)/(Ci*Gamma0) - &
                    as**4*beta0*Ci*Gamma2*Lperp**3 - &
                    (2*as**4*beta0**2*eta*Gamma2*Lperp**3)/Gamma0 - &
                    (5*as**4*beta0*beta1*gammai0*Lperp**3)/3._dp - &
                    2*as**4*beta0**2*gammai1*Lperp**3 - &
                    (as**5*beta1**2*Ci*Gamma0*logR*Lperp**3)/2._dp - &
                    as**5*beta0*beta2*Ci*Gamma0*logR*Lperp**3 - &
                    (7*as**5*beta0*beta1*Ci*Gamma1*logR*Lperp**3)/3._dp + &
                    4*as**5*beta0**3*gamma1rX*logR*Lperp**3 - &
                    2*as**5*beta0**2*Ci*Gamma2*logR*Lperp**3 - &
                    (35*as**5*beta0*beta1**2*eta*Lperp**4)/24._dp - &
                    (3*as**5*beta0**2*beta2*eta*Lperp**4)/2._dp - &
                    (3*as**5*beta1**2*Ci*Gamma0*Lperp**4)/8._dp - &
                    (3*as**5*beta0*beta2*Ci*Gamma0*Lperp**4)/4._dp - &
                    (7*as**5*beta0*beta1*Ci*Gamma1*Lperp**4)/4._dp - &
                    (47*as**5*beta0**2*beta1*eta*Gamma1*Lperp**4)/(12._dp*Gamma0) + &
                    4*as**5*beta0**3*gamma1rX*Lperp**4 + &
                    (5*as**5*beta0**4*eta*gamma1rX*Lperp**4)/(Ci*Gamma0) - &
                    (3*as**5*beta0**2*Ci*Gamma2*Lperp**4)/2._dp - &
                    (5*as**5*beta0**3*eta*Gamma2*Lperp**4)/(2._dp*Gamma0) - &
                    (13*as**5*beta0**2*beta1*gammai0*Lperp**4)/6._dp - &
                    2*as**5*beta0**3*gammai1*Lperp**4 - &
                    (13*as**5*beta0**2*beta1*Ci*Gamma0*logR*Lperp**4)/12._dp - &
                    as**5*beta0**3*Ci*Gamma1*logR*Lperp**4 - &
                    (77*as**5*beta0**3*beta1*eta*Lperp**5)/60._dp - &
                    (13*as**5*beta0**2*beta1*Ci*Gamma0*Lperp**5)/15._dp - &
                    (4*as**5*beta0**3*Ci*Gamma1*Lperp**5)/5._dp - &
                    (as**5*beta0**4*eta*Gamma1*Lperp**5)/Gamma0 - &
                    (2*as**5*beta0**4*gammai0*Lperp**5)/5._dp - &
                    (77*as**6*beta0**3*beta1*Ci*Gamma0*logR*Lperp**5)/60._dp - &
                    as**6*beta0**4*Ci*Gamma1*logR*Lperp**5 - &
                    (29*as**6*beta0**4*beta1*eta*Lperp**6)/20._dp - &
                    (77*as**6*beta0**3*beta1*Ci*Gamma0*Lperp**6)/72._dp - &
                    (5*as**6*beta0**4*Ci*Gamma1*Lperp**6)/6._dp - &
                    (as**6*beta0**5*eta*Gamma1*Lperp**6)/Gamma0 - &
                    (as**6*beta0**5*gammai0*Lperp**6)/3._dp - &
                    (as**6*beta0**5*Ci*Gamma0*logR*Lperp**6)/6._dp - &
                    (as**6*beta0**6*eta*Lperp**7)/7._dp - &
                    (as**6*beta0**5*Ci*Gamma0*Lperp**7)/7._dp - &
                    (as**7*beta0**6*Ci*Gamma0*logR*Lperp**7)/7._dp - &
                    (as**7*beta0**7*eta*Lperp**8)/8._dp - &
                    (as**7*beta0**6*Ci*Gamma0*Lperp**8)/8._dp
        endif

        fourierIntegrand_xt = 0.5_dp*xt*besselJ0(xt*qT)*exp(gi)*Lperp**real(nn,dp)

        !fourierIntegrand = aimag(-kv(0._dp, cmplx(xi*qT,0._dp,dp))*xi/pi*exp(gi)*Lperp**real(nn,dp))

    end function

    function fourierIntegrand(xi)
        use qtResummation_params, only: scalevar_rapidity, &
            scalevar_rapidity_i, scalevar_rapidity_mult
        implicit none
        real(dp), intent(in) :: xi
        real(dp) :: fourierIntegrand

        real(dp) :: aa, etaf
        complex(dp) :: Lperp

        real(dp) :: as, Ci, eta
        complex(dp) :: gi
        real(dp) :: gamma1rx, gamma2rx, gamma3rx
        real(dp) :: gammai0, gammai1, gammai2

        real(dp) :: logR

        if (initQuark) then
            Ci = CF
            gammai0 = gammaq0
            gammai1 = gammaq1
            gammai2 = gammaq2
            gamma1rx = -d2*CF
            gamma2rx = -d3*CF
            gamma3rx = -d4adj
        else
            Ci = CA
            gammai0 = gammag0
            gammai1 = gammag1
            gammai2 = gammag2
            gamma1rx = -d2*CA
            gamma2rx = -d3*CA
            gamma3rx = -d4fund
        endif

        aa = Ci*alphasMu/pi
        etaf = aa*log(Q2/mu**2)
        eta = etaf

        if (scalevar_rapidity_i > 0) then
            logR = log(scalevar_rapidity_mult(scalevar_rapidity_i))      
        else
            logR = 0._dp
        endif

        ! checked 9/26
        ! after replacement xT -> I*xT
        Lperp = cmplx(log(exp(2._dp*EulerGamma)*xi**2*mu**2/4._dp), pi, dp)

        as = alphasMu/4._dp/pi

        gi = -(eta*Lperp) - as*Ci*Gamma0*logR*Lperp -  &
        (as*beta0*eta*Lperp**2)/2._dp - (as*Ci*Gamma0*Lperp**2)/2._dp

        if (order >= 4) then
            gi = gi + (as*eta*gamma1rX)/(Ci*Gamma0) - (as*eta*Gamma1*Lperp)/Gamma0 - &
                2*as*gammai0*Lperp - as**2*Ci*Gamma1*logR*Lperp - &
                (as**2*beta1*eta*Lperp**2)/2._dp - (as**2*Ci*Gamma1*Lperp**2)/2._dp - &
                (as**2*beta0*eta*Gamma1*Lperp**2)/Gamma0 - &
                as**2*beta0*gammai0*Lperp**2 - &
                (as**2*beta0*Ci*Gamma0*logR*Lperp**2)/2._dp - &
                (as**2*beta0**2*eta*Lperp**3)/3._dp - &
                (as**2*beta0*Ci*Gamma0*Lperp**3)/3._dp - &
                (as**3*beta0**2*Ci*Gamma0*logR*Lperp**3)/3._dp - &
                (as**3*beta0**3*eta*Lperp**4)/4._dp - &
                (as**3*beta0**2*Ci*Gamma0*Lperp**4)/4._dp
        endif

        if (order >= 6) then
            gi = gi + (as**2*eta*gamma2rX)/(Ci*Gamma0) + as**2*gamma1rX*logR + &
                    as**2*gamma1rX*Lperp + &
                    (2*as**2*beta0*eta*gamma1rX*Lperp)/(Ci*Gamma0) - &
                    (as**2*eta*Gamma2*Lperp)/Gamma0 - 2*as**2*gammai1*Lperp + &
                    2*as**3*beta0*gamma1rX*logR*Lperp - &
                    as**3*Ci*Gamma2*logR*Lperp - (as**3*beta2*eta*Lperp**2)/2._dp - &
                    (as**3*beta1*eta*Gamma1*Lperp**2)/Gamma0 + &
                    2*as**3*beta0*gamma1rX*Lperp**2 + &
                    (3*as**3*beta0**2*eta*gamma1rX*Lperp**2)/(Ci*Gamma0) - &
                    (as**3*Ci*Gamma2*Lperp**2)/2._dp - &
                    (3*as**3*beta0*eta*Gamma2*Lperp**2)/(2._dp*Gamma0) - &
                    as**3*beta1*gammai0*Lperp**2 - 2*as**3*beta0*gammai1*Lperp**2 - &
                    (as**3*beta1*Ci*Gamma0*logR*Lperp**2)/2._dp - &
                    as**3*beta0*Ci*Gamma1*logR*Lperp**2 - &
                    (5*as**3*beta0*beta1*eta*Lperp**3)/6._dp - &
                    (as**3*beta1*Ci*Gamma0*Lperp**3)/3._dp - &
                    (2*as**3*beta0*Ci*Gamma1*Lperp**3)/3._dp - &
                    (as**3*beta0**2*eta*Gamma1*Lperp**3)/Gamma0 - &
                    (2*as**3*beta0**2*gammai0*Lperp**3)/3._dp - &
                    (5*as**4*beta0*beta1*Ci*Gamma0*logR*Lperp**3)/6._dp - &
                    as**4*beta0**2*Ci*Gamma1*logR*Lperp**3 - &
                    (13*as**4*beta0**2*beta1*eta*Lperp**4)/12._dp - &
                    (5*as**4*beta0*beta1*Ci*Gamma0*Lperp**4)/8._dp - &
                    (3*as**4*beta0**2*Ci*Gamma1*Lperp**4)/4._dp - &
                    (as**4*beta0**3*eta*Gamma1*Lperp**4)/Gamma0 - &
                    (as**4*beta0**3*gammai0*Lperp**4)/2._dp - &
                    (as**4*beta0**3*Ci*Gamma0*logR*Lperp**4)/4._dp - &
                    (as**4*beta0**4*eta*Lperp**5)/5._dp - &
                    (as**4*beta0**3*Ci*Gamma0*Lperp**5)/5._dp - &
                    (as**5*beta0**4*Ci*Gamma0*logR*Lperp**5)/5._dp - &
                    (as**5*beta0**5*eta*Lperp**6)/6._dp - &
                    (as**5*beta0**4*Ci*Gamma0*Lperp**6)/6._dp
        endif

         if (order >= 8) then
            gi = gi + (as**3*eta*gamma3rX)/(Ci*Gamma0) + as**3*gamma2rX*logR + &
                    (2*as**3*beta1*eta*gamma1rX*Lperp)/(Ci*Gamma0) + &
                    as**3*gamma2rX*Lperp + &
                    (3*as**3*beta0*eta*gamma2rX*Lperp)/(Ci*Gamma0) - &
                    (as**3*eta*Gamma3*Lperp)/Gamma0 - 2*as**3*gammai2*Lperp + &
                    2*as**4*beta1*gamma1rX*logR*Lperp + &
                    3*as**4*beta0*gamma2rX*logR*Lperp - &
                    as**4*Ci*Gamma3*logR*Lperp - (as**4*beta3*eta*Lperp**2)/2._dp - &
                    (as**4*beta2*eta*Gamma1*Lperp**2)/Gamma0 + &
                    2*as**4*beta1*gamma1rX*Lperp**2 + &
                    (7*as**4*beta0*beta1*eta*gamma1rX*Lperp**2)/(Ci*Gamma0) - &
                    (3*as**4*beta1*eta*Gamma2*Lperp**2)/(2._dp*Gamma0) + &
                    3*as**4*beta0*gamma2rX*Lperp**2 + &
                    (6*as**4*beta0**2*eta*gamma2rX*Lperp**2)/(Ci*Gamma0) - &
                    (as**4*Ci*Gamma3*Lperp**2)/2._dp - &
                    (2*as**4*beta0*eta*Gamma3*Lperp**2)/Gamma0 - &
                    as**4*beta2*gammai0*Lperp**2 - 2*as**4*beta1*gammai1*Lperp**2 - &
                    3*as**4*beta0*gammai2*Lperp**2 - &
                    (as**4*beta2*Ci*Gamma0*logR*Lperp**2)/2._dp - &
                    as**4*beta1*Ci*Gamma1*logR*Lperp**2 + &
                    3*as**4*beta0**2*gamma1rX*logR*Lperp**2 - &
                    (3*as**4*beta0*Ci*Gamma2*logR*Lperp**2)/2._dp - &
                    (as**4*beta1**2*eta*Lperp**3)/2._dp - &
                    as**4*beta0*beta2*eta*Lperp**3 - &
                    (as**4*beta2*Ci*Gamma0*Lperp**3)/3._dp - &
                    (2*as**4*beta1*Ci*Gamma1*Lperp**3)/3._dp - &
                    (7*as**4*beta0*beta1*eta*Gamma1*Lperp**3)/(3._dp*Gamma0) + &
                    3*as**4*beta0**2*gamma1rX*Lperp**3 + &
                    (4*as**4*beta0**3*eta*gamma1rX*Lperp**3)/(Ci*Gamma0) - &
                    as**4*beta0*Ci*Gamma2*Lperp**3 - &
                    (2*as**4*beta0**2*eta*Gamma2*Lperp**3)/Gamma0 - &
                    (5*as**4*beta0*beta1*gammai0*Lperp**3)/3._dp - &
                    2*as**4*beta0**2*gammai1*Lperp**3 - &
                    (as**5*beta1**2*Ci*Gamma0*logR*Lperp**3)/2._dp - &
                    as**5*beta0*beta2*Ci*Gamma0*logR*Lperp**3 - &
                    (7*as**5*beta0*beta1*Ci*Gamma1*logR*Lperp**3)/3._dp + &
                    4*as**5*beta0**3*gamma1rX*logR*Lperp**3 - &
                    2*as**5*beta0**2*Ci*Gamma2*logR*Lperp**3 - &
                    (35*as**5*beta0*beta1**2*eta*Lperp**4)/24._dp - &
                    (3*as**5*beta0**2*beta2*eta*Lperp**4)/2._dp - &
                    (3*as**5*beta1**2*Ci*Gamma0*Lperp**4)/8._dp - &
                    (3*as**5*beta0*beta2*Ci*Gamma0*Lperp**4)/4._dp - &
                    (7*as**5*beta0*beta1*Ci*Gamma1*Lperp**4)/4._dp - &
                    (47*as**5*beta0**2*beta1*eta*Gamma1*Lperp**4)/(12._dp*Gamma0) + &
                    4*as**5*beta0**3*gamma1rX*Lperp**4 + &
                    (5*as**5*beta0**4*eta*gamma1rX*Lperp**4)/(Ci*Gamma0) - &
                    (3*as**5*beta0**2*Ci*Gamma2*Lperp**4)/2._dp - &
                    (5*as**5*beta0**3*eta*Gamma2*Lperp**4)/(2._dp*Gamma0) - &
                    (13*as**5*beta0**2*beta1*gammai0*Lperp**4)/6._dp - &
                    2*as**5*beta0**3*gammai1*Lperp**4 - &
                    (13*as**5*beta0**2*beta1*Ci*Gamma0*logR*Lperp**4)/12._dp - &
                    as**5*beta0**3*Ci*Gamma1*logR*Lperp**4 - &
                    (77*as**5*beta0**3*beta1*eta*Lperp**5)/60._dp - &
                    (13*as**5*beta0**2*beta1*Ci*Gamma0*Lperp**5)/15._dp - &
                    (4*as**5*beta0**3*Ci*Gamma1*Lperp**5)/5._dp - &
                    (as**5*beta0**4*eta*Gamma1*Lperp**5)/Gamma0 - &
                    (2*as**5*beta0**4*gammai0*Lperp**5)/5._dp - &
                    (77*as**6*beta0**3*beta1*Ci*Gamma0*logR*Lperp**5)/60._dp - &
                    as**6*beta0**4*Ci*Gamma1*logR*Lperp**5 - &
                    (29*as**6*beta0**4*beta1*eta*Lperp**6)/20._dp - &
                    (77*as**6*beta0**3*beta1*Ci*Gamma0*Lperp**6)/72._dp - &
                    (5*as**6*beta0**4*Ci*Gamma1*Lperp**6)/6._dp - &
                    (as**6*beta0**5*eta*Gamma1*Lperp**6)/Gamma0 - &
                    (as**6*beta0**5*gammai0*Lperp**6)/3._dp - &
                    (as**6*beta0**5*Ci*Gamma0*logR*Lperp**6)/6._dp - &
                    (as**6*beta0**6*eta*Lperp**7)/7._dp - &
                    (as**6*beta0**5*Ci*Gamma0*Lperp**7)/7._dp - &
                    (as**7*beta0**6*Ci*Gamma0*logR*Lperp**7)/7._dp - &
                    (as**7*beta0**7*eta*Lperp**8)/8._dp - &
                    (as**7*beta0**6*Ci*Gamma0*Lperp**8)/8._dp
        endif

        fourierIntegrand = aimag(-kv(0._dp, cmplx(xi*qT,0._dp,dp))*xi/pi*exp(gi)*Lperp**real(nn,dp))

    end function

    ! transformed to interval 0,1
    function fourierIntegrand01(t)
        implicit none
        real(dp), intent(in) :: t
        real(dp) :: fourierIntegrand01

        if (t == 0._dp) then
            fourierIntegrand01 = 0._dp
            return
        endif

        if (t == 1._dp) then
            fourierIntegrand01 = 0._dp
            return
        endif

        fourierIntegrand01 = fourierIntegrand((1._dp-t)/t)/t**2


    end function

    function fourierIntegrandCos(theta)
        implicit none
        real(dp), intent(in) :: theta
        real(dp) :: fourierIntegrandCos

        fourierIntegrandCos = 2._dp*fourierIntegrand((cos(theta/2._dp)/sin(theta/2._dp))**2)*sin(theta)/(1._dp-cos(theta))**2

    end function

    function fourierIntegrandCos_xt(theta)
        implicit none
        real(dp), intent(in) :: theta
        real(dp) :: fourierIntegrandCos_xt

        fourierIntegrandCos_xt = 2._dp*fourierIntegrand_xt((cos(theta/2._dp)/sin(theta/2._dp))**2)*sin(theta)/(1._dp-cos(theta))**2

    end function

    ! checked 12/31/19 including order=1 against MM notebook
    function fourierM(nn_in, order_in, nf_in, initQuark_in, qt_in, q2_in, mu_in, alphasMu_in)
        use Quadpack, only : qagi, qag
        implicit none
        integer, intent(in) :: nn_in, order_in, nf_in
        logical, intent(in) :: initQuark_in
        real(dp), intent(in) :: qt_in, q2_in, mu_in, alphasMu_in
        real(dp) :: fourierM
        
        real(dp) :: absErr
        integer(int64) :: neval, ierr, neval2, neval3

        real(dp), parameter :: precisionGoal = 1e-12_dp
        real(dp) :: fourierM_xt, fourierM_xt2

        nn = nn_in
        order = order_in
        call update_nf_parameters(nf_in)
        initQuark = initQuark_in
        qt = qt_in
        q2 = q2_in
        mu = mu_in
        alphasMu = alphasMu_in

        if (nn < 5 .and. order < 8) then
            call qag(fourierIntegrandCos, 0._dp, pi, precisionGoal, precisionGoal, 1_int64, &
                fourierM, absErr, neval, ierr)
        else
            ! for higher Lperp powers this seems to be more stable
            ! for higher orders (8) this is actually necessary for stability
             call qagi(fourierIntegrand_xt, 0._dp, 1_int64, precisionGoal, precisionGoal, &
                 fourierM, absErr, neval, ierr)
        endif


    end function

    subroutine update_nf_parameters(nf_in)
        implicit none
        integer, intent(in) :: nf_in

        real(dp), parameter :: dAANA = 135._dp/8._dp
        real(dp), parameter :: dRANA = 15._dp/16._dp
        real(dp), parameter :: dRRNA = 5._dp/96._dp

        if (nf /= nf_in) then
            ! update all nf dependent constants
            nf = nf_in

            d2 = CA*(808._dp/27._dp - 28._dp*zeta3) - 224._dp/27._dp*tf*nf
            d3 = &
                (-10*CA*NF*(31313 - 618*Pi**2 + 27*Pi**4 - 12204*zeta3) + &
                  NF*(160*NF*(58 + 81*zeta3) + &
                     27*CF*(-8555 + 24*Pi**4 + 4560*zeta3)) + &
                  CA**2*(1485145 - 2079*Pi**4 - 1664280*zeta3 + &
                     60*Pi**2*(-799 + 594*zeta3) + 699840*zeta5))/3645._dp

            ! 2205.02249 eq. (6.10) and (6.11), adjusted by factor of (-2)
            d4adj = -2*(333.8_dp + 5506._dp*nf - 851.2_dp*nf**2 + 18.16_dp*nf**3)
            d4fund = -2*(-350.8_dp + 2428._dp*nf - 378.3_dp*nf**2 + 8.072_dp*nf**3)


            beta0 = 11._dp/3._dp*CA - 4._dp/3._dp*TF*nf
            beta1 = 34._dp/3._dp*CA**2 - 20._dp/3._dp*CA*TF*nf - 4._dp*CF*TF*nf
            beta2 = 2857._dp/54._dp*CA**3 + nf*(CF**2 - 205._dp/18._dp*CF*CA - &
                    1415._dp/54._dp*CA**2) + nf**2*(11._dp/9._dp*CF + 79._dp/54._dp*CA)
            beta3 = (CA**4*(150653 - 2376*zeta3) + 864*dAANA*(-5 + 132*zeta3) +  &
                  6*CA**3*nf*TF*(-39143 + 3672*zeta3) + &
                  4*nf*(CF*TF*(5589*CF**2 + 616*nf**2*TF**2 + &
                        36*CF*nf*TF*(169 - 264*zeta3)) + &
                     864*dRRNA*nf*(-11 + 24*zeta3) - 1728*dRANA*(-4 + 39*zeta3)) + &
                  8*CA*nf*TF*(106*nf**2*TF**2 + 16*CF*nf*TF*(268 + 189*zeta3) + &
                     9*CF**2*(-1051 + 264*zeta3)) + &
                  2*CA**2*nf*TF*(CF*(7073 - 17712*zeta3) + 6*nf*TF*(3965 + 1008*zeta3))) &
                 /486._dp

            Gamma0 = 4._dp
            Gamma1 = 4._dp*((67._dp/9._dp - pi**2/3._dp)*CA - 20._dp/9._dp*TF*nf)
            Gamma2 = 16._dp*((245._dp/24._dp - 67._dp/54._dp*pi**2 + 11._dp*pi**4/180._dp &
                    + 11._dp/6._dp*zeta3)*CA**2 - (209._dp/108._dp - 5._dp*pi**2/27._dp &
                    + 7._dp/3._dp*zeta3)*CA*nf - (55._dp/24._dp - 2._dp*zeta3)*CF*nf - nf**2/27._dp)
             ! 1911.10174 eq (6.4)
            Gamma3 = 15526.512384780493_dp - 3879.1186236243348_dp*nf + &
            146.68291933718706_dp*nf**2 + 2.454258338353606_dp*nf**3

            gammaq0 = -3._dp*CF
            gammaq1 = CF*nf*(130._dp/27._dp + 2._dp*pi**2/3._dp)*TF  &
                    + CF**2*(-3._dp/2._dp + 2._dp*pi**2 - 24._dp*zeta3) &
                    + CA*CF*(-961._dp/54._dp - 11._dp*pi**2/6._dp + 26._dp*zeta3)
            gammaq2 = (-139345*CA**2*CF)/2916._dp - (151*CA*CF**2)/4._dp - (29*CF**3)/2._dp - &
                (7163*CA**2*CF*Pi**2)/486._dp + (205*CA*CF**2*Pi**2)/9._dp - 3*CF**3*Pi**2 - &
                (83*CA**2*CF*Pi**4)/90._dp + (247*CA*CF**2*Pi**4)/135._dp - &
                (8*CF**3*Pi**4)/5._dp - (17318*CA*CF*nf*TF)/729._dp + &
                (2953*CF**2*nf*TF)/27._dp + (2594*CA*CF*nf*Pi**2*TF)/243._dp - &
                (26*CF**2*nf*Pi**2*TF)/9._dp + (22*CA*CF*nf*Pi**4*TF)/45._dp - &
                (28*CF**2*nf*Pi**4*TF)/27._dp + (9668*CF*nf**2*TF**2)/729._dp - &
                (40*CF*nf**2*Pi**2*TF**2)/27._dp + (3526*CA**2*CF*zeta3)/9._dp - &
                (844*CA*CF**2*zeta3)/3._dp - 68*CF**3*zeta3 - &
                (44*CA**2*CF*Pi**2*zeta3)/9._dp - (8*CA*CF**2*Pi**2*zeta3)/3._dp + &
                (16*CF**3*Pi**2*zeta3)/3._dp - (1928*CA*CF*nf*TF*zeta3)/27._dp + &
                (512*CF**2*nf*TF*zeta3)/9._dp - (32*CF*nf**2*TF**2*zeta3)/27._dp - &
                136*CA**2*CF*zeta5 - 120*CA*CF**2*zeta5 + 240*CF**3*zeta5

            gammag0 = (-11*CA)/3._dp + (4*nf*TF)/3._dp
            gammag1 = (-692*CA**2)/27._dp + (11*CA**2*Pi**2)/18._dp + (256*CA*nf*TF)/27._dp + &
                4*CF*nf*TF - (2*CA*nf*Pi**2*TF)/9._dp + 2*CA**2*zeta3
            gammag2 = (-97186*CA**3)/729._dp + (6109*CA**3*Pi**2)/486._dp - &
                (319*CA**3*Pi**4)/270._dp + (30715*CA**2*nf*TF)/729._dp + &
                (2434*CA*CF*nf*TF)/27._dp - 2*CF**2*nf*TF - &
                (1198*CA**2*nf*Pi**2*TF)/243._dp - (2*CA*CF*nf*Pi**2*TF)/3._dp + &
                (82*CA**2*nf*Pi**4*TF)/135._dp - (8*CA*CF*nf*Pi**4*TF)/45._dp - &
                (538*CA*nf**2*TF**2)/729._dp - (44*CF*nf**2*TF**2)/9._dp + &
                (40*CA*nf**2*Pi**2*TF**2)/81._dp + (122*CA**3*zeta3)/3._dp - &
                (20*CA**3*Pi**2*zeta3)/9._dp + (712*CA**2*nf*TF*zeta3)/27._dp - &
                (304*CA*CF*nf*TF*zeta3)/9._dp - (224*CA*nf**2*TF**2*zeta3)/27._dp - &
                16*CA**3*zeta5
        endif

    end subroutine


end module
