
!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module qtResummationHardEvolution
      use types
      use constants
      use iso_fortran_env
      implicit none

      public :: hardEvolutionThres, hardEvolution

      private

      real(dp), save :: nf = 0
!$omp threadprivate(nf)

      real(dp), save :: beta0,beta1,beta2,beta3,beta4
!$omp threadprivate(beta0,beta1,beta2,beta3,beta4)

      real(dp), save :: Gamma0, Gamma1, Gamma2, Gamma3, Gamma4
!$omp threadprivate(Gamma0, Gamma1,Gamma2,Gamma3,Gamma4)

      real(dp), save :: gammaq0,gammaq1, gammaq2, gammaq3, gammag0, gammag1, gammag2, gammag3
!$omp threadprivate(gammaq0,gammaq1,gammaq2,gammaq3,gammag0,gammag1,gammag2,gammag3)

      contains

    ! takes into account thresholds
    function hardEvolutionThres(q2,muH2,alphasMuH,alphasMu,order,initQuark, &
          nf_in, nfmax, mb, alphas4p, alphas4m, mc, alphas3p, alphas3m)
        implicit none
        real(dp), intent(in) :: q2, muH2, alphasMuH, alphasMu
        integer, intent(in) :: order
        logical, intent(in) :: initQuark

        integer, intent(in) :: nf_in, nfmax
        real(dp), intent(in), optional :: mb, alphas4p, alphas4m
        real(dp), intent(in), optional :: mc, alphas3p, alphas3m

        real(dp) :: hardEvolutionThres

        ! nf_in tells us which in which nf-flavor region mu and therefore alphasMu is
        ! nfmax tells us nf(muH), i.e. the point from which to evolve downwards

        hardEvolutionThres = 0._dp

        if (nf_in > nfmax) then
            write (*,*) __FILE__//": nf > nfmax, this should not happen"
            error stop
        elseif (nfmax <= 2) then
            write (*,*) __FILE__//": nfmax <= 2, please fix PDF"
            error stop
        endif

        if (nf_in == 5 .and. nfmax > 4) then
            ! just evolve from muH to mu
            call update_nf_parameters(5)
            hardEvolutionThres = hardEvolution(q2,muH2,alphasMuH,alphasMu,order,initQuark)
            return
        elseif (nf_in == 4 .and. nfmax > 3) then
            if ((.not. present(mb)) .or. (.not. present(alphas4p)) .or. &
                (.not. present(alphas4m))) then
                write (*,*) __FILE__//": mb,alphas4p,alphas4m must be present for nf=4 evolution"
                error stop
            endif

            ! run down to above mb threshold (alphas4p)
            if (nfmax > 4) then
                call update_nf_parameters(5)
                hardEvolutionThres = hardEvolution(q2,muH2,alphasMuH,alphas4p,order,initQuark)
            endif
            ! run from below mb threshold as hard scale with alphas4m to mu (alphasMu)
            if (nfmax > 3) then
                call update_nf_parameters(4)
                hardEvolutionThres = hardEvolutionThres * &
                    hardEvolution(q2,mb**2,alphas4m,alphasMu,order,initQuark)
            endif
            return
        elseif (nf_in == 3 .and. nfmax > 2) then
            if ((.not. present(mb)) .or. (.not. present(alphas4p)) .or. &
                (.not. present(alphas4m)) .or. (.not. present(mc)) .or. &
                (.not. present(alphas3p)) .or. (.not. present(alphas3m))) then
                write (*,*) __FILE__//": mb,alphas4p,alpha4m must be present for nf=3 evolution"
                write (*,*) __FILE__//": mc,alphas3p,alpha3m must be present for nf=3 evolution"
                error stop
            endif

            ! run down to above mb threshold (alphas4p)
            if (nfmax > 4) then
                call update_nf_parameters(5)
                hardEvolutionThres = hardEvolution(q2,muH2,alphasMuH,alphas4p,order,initQuark)
            endif
            ! run from below mb threshold as hard scale with alphas4m to alphas3p
            if (nfmax > 3) then
                call update_nf_parameters(4)
                hardEvolutionThres = hardEvolutionThres * &
                    hardEvolution(q2,mb**2,alphas4m,alphas3p,order,initQuark)
            endif
            ! run from below mc threshold as hard scale with alphas3m to alphasMu
            if (nfmax > 2) then
                call update_nf_parameters(3)
                hardEvolutionThres = hardEvolutionThres * &
                    hardEvolution(q2,mc**2,alphas3m,alphasMu,order,initquark)
            endif
            return
        else
            write (*,*) __FILE__//": Nf = ", nf, " not supported in hard evolution"
            write (*,*) __FILE__//": nfmax = ", nfmax
            error stop
        endif

    end function hardEvolutionThres

    function hardEvolution(q2,muH2,alphasMuH,alphasMu,order,initQuark)
        implicit none
        real(dp), intent(in) :: q2, muH2, alphasMuH, alphasMu
        integer, intent(in) :: order
        logical, intent(in) :: initQuark
        real(dp) :: hardEvolution

        real(dp) :: a_anom, a_cusp, S, r

        real(dp) :: gammai0, gammai1, gammai2, gammai3

        a_cusp = (Gamma0*Log(alphasMu/alphasMuH))/(2._dp*beta0)
        if (order >= 4) then
            a_cusp = a_cusp + ((alphasMu - alphasMuH)*(-(beta1*Gamma0) + beta0*Gamma1))/(8._dp*beta0**2*Pi)
        endif
        if (order >= 6) then
            a_cusp = a_cusp + ((alphasMu**2 - alphasMuH**2)* &
          (beta1**2*Gamma0 - beta0*beta2*Gamma0 - beta0*beta1*Gamma1 + &
            beta0**2*Gamma2))/(64._dp*beta0**3*Pi**2)
        endif
        if (order >= 8) then
            a_cusp = a_cusp &
             -((alphasMu**3 - alphasMuH**3)*(beta1**3*Gamma0 - beta0*beta1**2*Gamma1 + &
              beta0*beta1*(-2*beta2*Gamma0 + beta0*Gamma2) + &
              beta0**2*(beta3*Gamma0 + beta2*Gamma1 - beta0*Gamma3)))/ &
              (384.*beta0**4*Pi**3)
        endif

        if (initQuark) then
            gammai0 = gammaq0
            gammai1 = gammaq1
            gammai2 = gammaq2
            gammai3 = gammaq3
        else
            gammai0 = gammag0
            gammai1 = gammag1
            gammai2 = gammag2
            gammai3 = gammag3
        endif

        a_anom = (gammai0*Log(alphasMu/alphasMuH))/(2._dp*beta0)
        if (order >= 4) then
            a_anom = a_anom + ((alphasMu - alphasMuH)*(-(beta1*gammai0) + beta0*gammai1))/(8._dp*beta0**2*Pi)
        endif
        if (order >= 6) then
            a_anom = a_anom + ((alphasMu**2 - alphasMuH**2)* &
          (beta1**2*gammai0 - beta0*beta2*gammai0 - beta0*beta1*gammai1 + &
            beta0**2*gammai2))/(64._dp*beta0**3*Pi**2)
        endif
        if (order >= 8) then
            a_anom = a_anom &
             -((alphasMu**3 - alphasMuH**3)*(beta1**3*gammai0 - beta0*beta1**2*gammai1 + &
              beta0*beta1*(-2*beta2*gammai0 + beta0*gammai2) + &
              beta0**2*(beta3*gammai0 + beta2*gammai1 - beta0*gammai3)))/ &
              (384.*beta0**4*Pi**3)
        endif

        S = ((-1 + alphasMu/alphasMuH)*(beta1*Gamma0 - beta0*Gamma1))/ &
         (4._dp*beta0**3) + ((-(beta1*Gamma0) + beta0*Gamma1)* &
           Log(alphasMu/alphasMuH))/(4._dp*beta0**3) + &
        (beta1*Gamma0*Log(alphasMu/alphasMuH)**2)/(8._dp*beta0**3) + &
        (Gamma0*Pi*(-1 + alphasMu/alphasMuH - &
             (alphasMu*Log(alphasMu/alphasMuH))/alphasMuH))/(alphasMu*beta0**2)

         if (order >= 4) then
             S = S + (-((alphasMu - alphasMuH)* &
             (-(alphasMuH*(beta1**2*Gamma0 + beta0*beta2*Gamma0 - &
                    3*beta0*beta1*Gamma1 + beta0**2*Gamma2)) + &
               alphasMu*(beta1**2*Gamma0 - beta0*beta2*Gamma0 - &
                  beta0*beta1*Gamma1 + beta0**2*Gamma2))) + &
          2*alphasMuH*(alphasMuH*(beta1**2 - beta0*beta2)*Gamma0 + &
             alphasMu*beta1*(-(beta1*Gamma0) + beta0*Gamma1))*&
           Log(alphasMu/alphasMuH))/(32._dp*alphasMuH*beta0**4*Pi)
         endif

         if (order >= 6) then
            S = S + ((alphasMu - alphasMuH)*(4*alphasMu**2* &
              (beta1**3*Gamma0 - beta0*beta1**2*Gamma1 + &
                beta0*beta1*(-2*beta2*Gamma0 + beta0*Gamma2) + &
                beta0**2*(beta3*Gamma0 + beta2*Gamma1 - beta0*Gamma3)) + &
             alphasMuH**2*(-8*beta1**3*Gamma0 + 11*beta0*beta1**2*Gamma1 + &
                beta0*beta1*(7*beta2*Gamma0 - 5*beta0*Gamma2) + &
                beta0**2*(beta3*Gamma0 - 8*beta2*Gamma1 + 2*beta0*Gamma3)) + &
             alphasMu*alphasMuH*(4*beta1**3*Gamma0 - beta0*beta1**2*Gamma1 - &
                5*beta0*beta1*(beta2*Gamma0 + beta0*Gamma2) + &
                beta0**2*(beta3*Gamma0 + 4*beta2*Gamma1 + 2*beta0*Gamma3))) - &
                6*alphasMuH*(alphasMuH**2* &
              (beta1**3 - 2*beta0*beta1*beta2 + beta0**2*beta3)*Gamma0 + &
             alphasMu**2*beta1*(-(beta1**2*Gamma0) + beta0*beta2*Gamma0 + &
                beta0*beta1*Gamma1 - beta0**2*Gamma2))*Log(alphasMu/alphasMuH))/ &
                (768._dp*alphasMuH*beta0**5*Pi**2)
         endif
         
         if (order >= 8) then
            r = alphasMu/alphasMuH
            S = S + alphasMuH**3*(-3*beta1**4*Gamma0*(-1 + r)**2*(7 + r*(8 + 3*r)) + &
             beta0*beta1**3*Gamma1*(-1 + r)*(-25 + r*(-7 + r*(11 + 9*r))) - &
             beta0**2*beta1*(-1 + r)* &
              (2*beta3*Gamma0*(-7 + 9*r)*(1 + r + r**2) + &
                2*beta2*Gamma1*(-20 + r*(1 + r)*(-2 + 9*r)) + &
                beta0*Gamma3*(7 + r*(7 + (7 - 9*r)*r))) + &
             beta0*beta1**2*(-1 + r)* &
              (beta0*Gamma2*(13 + r*(13 - r*(5 + 9*r))) + &
                beta2*Gamma0*(-41 + r*(-5 + r*(31 + 27*r)))) + &
             beta0**2*(9*beta0*beta2*Gamma2*(-1 + r**2)**2 - &
                beta2**2*Gamma0*(5 + r**2*(-18 + r*(4 + 9*r))) + &
                beta0*(beta4*Gamma0*(-1 + r**3*(-8 + 9*r)) - &
                   3*(beta3*Gamma1*(-5 + 6*r + 2*r**3 - 3*r**4) + &
                      beta0*Gamma4*(1 + r**3*(-4 + 3*r))))) - &
             12*(-((beta1**4 - 3*beta0*beta1**2*beta2 + &
                     2*beta0**2*beta1*beta3 + &
                     beta0**2*(beta2**2 - beta0*beta4))*Gamma0) + &
                beta1*(beta1**3*Gamma0 - beta0*beta1**2*Gamma1 + &
                   beta0*beta1*(-2*beta2*Gamma0 + beta0*Gamma2) + &
                   beta0**2*(beta3*Gamma0 + beta2*Gamma1 - beta0*Gamma3))*r**3 &
                )*Log(r))/(9216.*beta0**6*Pi**3)
         endif


        if (initQuark) then
            hardEvolution = exp(4._dp*CF*S - 4._dp*a_anom)*(Q2/muh2)**(-2._dp*CF*a_cusp)
        else
            hardEvolution = exp(4._dp*CA*S - 4._dp*a_anom)*(Q2/muh2)**(-2._dp*CA*a_cusp)
        endif

    end function

    subroutine update_nf_parameters(nf_in)
        implicit none
        integer, intent(in) :: nf_in

        real(dp), parameter :: dAANA = 135._dp/8._dp
        real(dp), parameter :: dRANA = 15._dp/16._dp
        real(dp), parameter :: dRRNA = 5._dp/96._dp

        if (nf /= real(nf_in,dp)) then
            ! update all nf dependent constants
            nf = nf_in

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

! Expressions beta3alt and beta4 taken from arXiv:1606.08659, Eqs. (4) and (5)
! Note that we want the powers of 4 in the denominator (as/4/pi) so they are removed
!	    beta3alt = (149753._dp/6._dp+3564*zeta3-(1078361._dp/162._dp+6508._dp/27._dp*zeta3)*nf &
!	               +(50065._dp/162._dp+6472._dp/81._dp*zeta3)*nf**2+1093._dp/729._dp*nf**3)
		       
	    beta4 = (8157455._dp/16._dp+621885._dp/2._dp*zeta3-88209._dp/2._dp*zeta4-288090._dp*zeta5 &
	            +nf*(-336460813._dp/1944._dp-4811164._dp/81._dp*zeta3 &
		         +33935._dp/6._dp*zeta4+1358995._dp/27._dp*zeta5) &
		    +nf**2*(25960913._dp/1944._dp+698531._dp/81._dp*zeta3-10526._dp/9._dp*zeta4-381760._dp/81._dp*zeta5) &
		    +nf**3*(-630559._dp/5832._dp-48722._dp/243._dp*zeta3+1618._dp/27._dp*zeta4+460._dp/9._dp*zeta5) &
		    +nf**4*(1205._dp/2916._dp-152._dp/81._dp*zeta3))

            Gamma0 = 4._dp
            Gamma1 = 4._dp*((67._dp/9._dp - pi**2/3._dp)*CA - 20._dp/9._dp*TF*nf)
            Gamma2 = 16._dp*((245._dp/24._dp - 67._dp/54._dp*pi**2 + 11._dp*pi**4/180._dp &
                    + 11._dp/6._dp*zeta3)*CA**2 - (209._dp/108._dp - 5._dp*pi**2/27._dp &
                    + 7._dp/3._dp*zeta3)*CA*nf - (55._dp/24._dp - 2._dp*zeta3)*CF*nf - nf**2/27._dp)
             ! 1911.10174 eq (6.4)
            Gamma3 = 15526.512384780493_dp - 3879.1186236243348_dp*nf + &
            146.68291933718706_dp*nf**2 + 2.454258338353606_dp*nf**3

! Numerical values taken from arXiv:1812.11818, Eq. (13)
!            Gamma0=0.42441*fourpi/CF
            if (nf_in == 3) then
!              Gamma1=0.42441*0.7266*fourpi**2/CF
!              Gamma2=0.42441*0.7341*fourpi**3/CF
!              Gamma3=0.42441*0.665*fourpi**4/CF
              Gamma4=0.42441*1.3*fourpi**5/CF
            elseif (nf_in == 4) then
!              Gamma1=0.42441*0.6382*fourpi**2/CF
!              Gamma2=0.42441*0.5100*fourpi**3/CF
!              Gamma3=0.42441*0.317*fourpi**4/CF
              Gamma4=0.42441*0.8*fourpi**5/CF
            elseif (nf_in == 5) then
!              Gamma1=0.42441*0.5497*fourpi**2/CF
!              Gamma2=0.42441*0.2840*fourpi**3/CF
!              Gamma3=0.42441*0.013*fourpi**4/CF
              Gamma4=0.42441*0.5*fourpi**5/CF
             else
              write(6,*) 'Gamma4 not defined for nf = ',nf_in
              call exit(1)
             endif

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
                
! Numerical value taken from 2102.09725, Eq. (22)
! Note that we need an overall factor of (-1/2), c.f. comment in 2002.04617 before Eq. (8)
            gammaq3 = -0.5_dp*( &
               + dRANA * ( - 2195.670535473) &
               + CA**3*CF * ( - 13.809312037) &
               + CA**2*CF**2 * ( 2438.569338812) &
               + CA*CF**3 * (- 1373.764650948) &
               +CF**4  * ( 392.899478384) &
               +nf*dRRNA * ( - 425.019550390) &
               +nf*CA**2*CF * ( - 274.147360589) &
               +nf*CA*CF**2 * ( - 912.844845636) &
               +nf*CF**3  * ( 151.933788877) &
               +nf**2*CA*CF  * ( 109.081415293) &
               +nf**2*CF**2 * ( - 12.5342425083) &
               +nf**3*CF  * ( 4.88682798281) &
               )
 
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

! Numerical value taken from 2102.09725, Eq. (23)
! Note that we need an overall factor of (-1/2), c.f. comment in 2002.04617 before Eq. (8)
            gammag3 = -0.5_dp*( &
               + dAANA * ( - 2451.040712450) &
               +CA**4 * ( 1557.4287417889) &
               +nf*dRANA * ( - 41.2080190194) &
               +nf*CA**3 * ( - 1033.98729659) &
               +nf*CA**2*CF * ( - 57.9377499658) &
               +nf*CA*CF**2 * ( - 100.315097910) &
               +nf*CF**3  * ( 46) &
               +nf**2*dRRNA * ( 253.857645167) &
               +nf**2*CA**2  * ( 70.7744401902) &
               +nf**2*CA*CF  * ( 73.9372035966) &
               +nf**2*CF**2 * ( - 21.9767440643) &
               +nf**3*CA  * ( 0.405507202650) &
               +nf**3*CF  * ( 1.26748971193) &
               )

        endif

    end subroutine

end module

