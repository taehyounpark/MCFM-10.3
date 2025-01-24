!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

!
! "The Diphoton qT spectrum at N3LL'+NNLO", Tobias Neumann, 2107.12478
!

module GammaGamma3L
    use types
    use GammaGammaABCD
    implicit none

    public :: gamgamhard

    private

    contains

    subroutine gamgamhard(p,tree,hard,twoloop_nfv2,threeloop_nfv,threeloop_nfv2)
        use types
        use constants
        implicit none
        include 'nf.f'
        include 'mxpart.f'
        include 'zprods_decl.f'
        include 'sprods_com.f'
        include 'scale.f'

        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(out) :: tree, hard(3)
        real(dp), intent(out) :: twoloop_nfv2
        real(dp), intent(out) :: threeloop_nfv, threeloop_nfv2

        complex(dp) :: Lmm(0:3), Lmp(0:3), Lpm(0:3), Lpp(0:3)
        complex(dp) :: Lmm_nfv, Lmp_nfv, Lpm_nfv, Lpp_nfv
        complex(dp) :: Lmm_nfv2(0:3), Lmp_nfv2(0:3), Lpm_nfv2(0:3), Lpp_nfv2(0:3)

        complex(dp) :: c_alpha(0:3), c_beta(0:3), c_gamma(0:3), c_delta(0:3)
        complex(dp) :: c_alpha_NfV, c_beta_NfV, c_gamma_NfV, c_delta_NfV
        complex(dp) :: c_alpha_NfV2(0:3), c_beta_NfV2(0:3), c_gamma_NfV2(0:3), c_delta_NfV2(0:3)
        real(dp) :: xx
        real(dp) :: Ls12mu

        c_alpha(:) = 0._dp
        c_beta(:) = 0._dp
        c_gamma(:) = 0._dp
        c_delta(:) = 0._dp

        c_alpha_NfV = 0._dp
        c_beta_NfV = 0._dp
        c_gamma_NfV = 0._dp
        c_delta_NfV = 0._dp

        c_alpha_NfV2(:) = 0._dp
        c_beta_NfV2(:) = 0._dp
        c_gamma_NfV2(:) = 0._dp
        c_delta_NfV2(:) = 0._dp


        call spinoru(4,p,za,zb)

        xx = -s(1,3)/s(1,2)

        ! tree level
        c_alpha(0) = 0._dp
        c_beta(0) = 1._dp
        c_gamma(0) = 1._dp

        ! one loop
        c_alpha(1) = alpha_1l(xx)
        c_beta(1) = beta_1l(xx)
        c_gamma(1) = gamma_1l(xx)

        call inithpl(xx)

        ! two loop
        c_alpha(2) = alphaHPL_2l_noNfV2(xx)
        c_beta(2) = betaHPL_2l_noNfV2(xx)
        c_gamma(2) = gammaHPL_2l_noNfV2(xx)

        ! three loop
        c_alpha(3) = alphaHPL_3l_no_nfv(xx)
        c_beta(3) = betaHPL_3l_no_nfv(xx)
        c_gamma(3) = gammaHPL_3l_no_nfv(xx)

        ! delta by symmetry
        c_delta(:) = -c_alpha(:)

        Lmm = 2._dp*za(3,4)**2/za(1,3)/zb(2,3) * c_alpha
        Lmp = 2._dp*za(2,4)*zb(1,3)/za(2,3)/zb(2,4) * c_beta
        Lpm = 2._dp*za(2,3)*zb(4,1)/za(2,4)/zb(3,2) * c_gamma
        Lpp = 2._dp*za(3,4)**2/za(3,1)/zb(2,3) * c_delta

        ! reinstante renormalization scale dependence
        Ls12mu = log(s(1,2)/musq)

        call rge_evolve(Ls12mu, Lmm)
        call rge_evolve(Ls12mu, Lmp)
        call rge_evolve(Ls12mu, Lpm)
        call rge_evolve(Ls12mu, Lpp)
 
        tree = abs(Lmm(0))**2 + abs(Lmp(0))**2 + abs(Lpm(0))**2 + abs(Lpp(0))**2

        hard(1) = 2._dp*real(Lmm(0)*conjg(Lmm(1))) + &
                  2._dp*real(Lmp(0)*conjg(Lmp(1))) + &
                  2._dp*real(Lpm(0)*conjg(Lpm(1))) + &
                  2._dp*real(Lpp(0)*conjg(Lpp(1)))

        hard(2) = 2._dp*real(Lmm(0)*conjg(Lmm(2))) + &
                  2._dp*real(Lmp(0)*conjg(Lmp(2))) + &
                  2._dp*real(Lpm(0)*conjg(Lpm(2))) + &
                  2._dp*real(Lpp(0)*conjg(Lpp(2))) + &
                        real(Lmm(1)*conjg(Lmm(1))) + &
                        real(Lmp(1)*conjg(Lmp(1))) + &
                        real(Lpm(1)*conjg(Lpm(1))) + &
                        real(Lpp(1)*conjg(Lpp(1)))

        hard(3) = 2._dp*real(Lmm(0)*conjg(Lmm(3))) + &
                  2._dp*real(Lmp(0)*conjg(Lmp(3))) + &
                  2._dp*real(Lpm(0)*conjg(Lpm(3))) + &
                  2._dp*real(Lpp(0)*conjg(Lpp(3))) + &
                  2._dp*real(Lmm(1)*conjg(Lmm(2))) + &
                  2._dp*real(Lmp(1)*conjg(Lmp(2))) + &
                  2._dp*real(Lpm(1)*conjg(Lpm(2))) + &
                  2._dp*real(Lpp(1)*conjg(Lpp(2)))

        ! two loop and three loop nfv2 piece
        c_alpha_nfv2(2) = alphaHPL_2l_NfV2(xx)
        c_beta_nfv2(2) = betaHPL_2l_NfV2(xx)
        c_gamma_nfv2(2) = gammaHPL_2l_NfV2(xx)

        c_alpha_nfv2(3) = alphaHPL_3l_nfv2(xx)
        c_beta_nfv2(3) = betaHPL_3l_nfv2(xx)
        c_gamma_nfv2(3) = gammaHPL_3l_nfv2(xx)

        ! delta by symmetry
        c_delta_nfv2(:) = -c_alpha_nfv2(:)

        Lmm_nfv2 = 2._dp*za(3,4)**2/za(1,3)/zb(2,3) * c_alpha_nfv2
        Lmp_nfv2 = 2._dp*za(2,4)*zb(1,3)/za(2,3)/zb(2,4) * c_beta_nfv2
        Lpm_nfv2 = 2._dp*za(2,3)*zb(4,1)/za(2,4)/zb(3,2) * c_gamma_nfv2
        Lpp_nfv2 = 2._dp*za(3,4)**2/za(3,1)/zb(2,3) * c_delta_nfv2

        call rge_evolve(Ls12mu, Lmm_nfv2)
        call rge_evolve(Ls12mu, Lmp_nfv2)
        call rge_evolve(Ls12mu, Lpm_nfv2)
        call rge_evolve(Ls12mu, Lpp_nfv2)

        twoloop_nfv2 = 2._dp*real(Lmm(0)*conjg(Lmm_nfv2(2))) + &
                       2._dp*real(Lmp(0)*conjg(Lmp_nfv2(2))) + &
                       2._dp*real(Lpm(0)*conjg(Lpm_nfv2(2))) + &
                       2._dp*real(Lpp(0)*conjg(Lpp_nfv2(2)))

        threeloop_nfv2 = 2._dp*real(Lmm(1)*conjg(Lmm_nfv2(2))) + &
                         2._dp*real(Lmp(1)*conjg(Lmp_nfv2(2))) + &
                         2._dp*real(Lpm(1)*conjg(Lpm_nfv2(2))) + &
                         2._dp*real(Lpp(1)*conjg(Lpp_nfv2(2))) + &
                         2._dp*real(Lmm(0)*conjg(Lmm_nfv2(3))) + &
                         2._dp*real(Lmp(0)*conjg(Lmp_nfv2(3))) + &
                         2._dp*real(Lpm(0)*conjg(Lpm_nfv2(3))) + &
                         2._dp*real(Lpp(0)*conjg(Lpp_nfv2(3)))

        ! three loop nfv piece
        c_alpha_nfv = alphaHPL_3l_nfv(xx)
        c_beta_nfv = betaHPL_3l_nfv(xx)
        c_gamma_nfv = gammaHPL_3l_nfv(xx)
        c_delta_nfv = -c_alpha_nfv

        Lmm_nfv  = 2._dp*za(3,4)**2/za(1,3)/zb(2,3) * c_alpha_nfv
        Lmp_nfv  = 2._dp*za(2,4)*zb(1,3)/za(2,3)/zb(2,4) * c_beta_nfv
        Lpm_nfv  = 2._dp*za(2,3)*zb(4,1)/za(2,4)/zb(3,2) * c_gamma_nfv
        Lpp_nfv  = 2._dp*za(3,4)**2/za(3,1)/zb(2,3) * c_delta_nfv

        ! no musq dependence for nfv piece since lower orders are zero

        threeloop_nfv = 2._dp*real(Lmm(0)*conjg(Lmm_nfv)) + &
                        2._dp*real(Lmp(0)*conjg(Lmp_nfv)) + &
                        2._dp*real(Lpm(0)*conjg(Lpm_nfv)) + &
                        2._dp*real(Lpp(0)*conjg(Lpp_nfv))

        ! since we only included left-handed quarks
        ! right-handed are only a factor of two
        tree = tree * 2._dp
        hard(:) = hard(:) * 2._dp
        twoloop_nfv2 = twoloop_nfv2 * 2._dp
        threeloop_nfv2 = threeloop_nfv2 * 2._dp
        threeloop_nfv = threeloop_nfv * 2._dp

    end subroutine gamgamhard

    subroutine rge_evolve(Ls12mu,amp)
        use types
        use constants
        implicit none
        real(dp), intent(in) :: Ls12mu
        complex(dp), intent(inout) :: amp(0:3)

        include 'nf.f'

        real(dp) :: beta0, beta1
        real(dp) :: gamma0q, gamma1q, gamma2q
        real(dp) :: gamma0cusp, gamma1cusp, gamma2cusp

        gamma0q = -3._dp*CF
        gamma1q = CF*nf*(130._dp/27._dp + 2._dp*pi**2/3._dp)*TF  &
                + CF**2*(-3._dp/2._dp + 2._dp*pi**2 - 24._dp*zeta3) &
                + CA*CF*(-961._dp/54._dp - 11._dp*pi**2/6._dp + 26._dp*zeta3)

        gamma0cusp = 4._dp
        gamma1cusp = 4._dp*((67._dp/9._dp - pi**2/3._dp)*CA - 20._dp/9._dp*TF*nf)
        gamma2cusp = 16._dp*((245._dp/24._dp - 67._dp/54._dp*pi**2 + 11._dp*pi**4/180._dp &
                    + 11._dp/6._dp*zeta3)*CA**2 - (209._dp/108._dp - 5._dp*pi**2/27._dp &
                    + 7._dp/3._dp*zeta3)*CA*nf - (55._dp/24._dp - 2._dp*zeta3)*CF*nf - nf**2/27._dp)

        beta0 = 11._dp/3._dp*CA - 4._dp/3._dp*TF*nf
        beta1 = 34._dp/3._dp*CA**2 - 20._dp/3._dp*CA*TF*nf - 4._dp*CF*TF*nf

        if (amp(0) == 0._dp .and. amp(1) == 0._dp) then
            amp(3) = (-8*beta0*Ls12mu*amp(2) - 4*gamma0q*Ls12mu*amp(2) - &
                CF*gamma0cusp*Ls12mu**2*amp(2) + &
                2*CF*gamma0cusp*im*Ls12mu*Pi*amp(2) + 8*amp(3))/8._dp
        else
            amp(3) = (-384*gamma2q*Ls12mu*amp(0) + 384*beta0*gamma1q*Ls12mu**2*amp(0) + &
                384*gamma0q*gamma1q*Ls12mu**2*amp(0) - &
                96*CF*gamma2cusp*Ls12mu**2*amp(0) - &
                128*beta0**2*gamma0q*Ls12mu**3*amp(0) - &
                192*beta0*gamma0q**2*Ls12mu**3*amp(0) - &
                64*gamma0q**3*Ls12mu**3*amp(0) + &
                64*beta0*CF*gamma1cusp*Ls12mu**3*amp(0) + &
                96*CF*gamma0q*gamma1cusp*Ls12mu**3*amp(0) + &
                96*CF*gamma0cusp*gamma1q*Ls12mu**3*amp(0) - &
                16*beta0**2*CF*gamma0cusp*Ls12mu**4*amp(0) - &
                80*beta0*CF*gamma0cusp*gamma0q*Ls12mu**4*amp(0) - &
                48*CF*gamma0cusp*gamma0q**2*Ls12mu**4*amp(0) + &
                24*CF**2*gamma0cusp*gamma1cusp*Ls12mu**4*amp(0) - &
                8*beta0*CF**2*gamma0cusp**2*Ls12mu**5*amp(0) - &
                12*CF**2*gamma0cusp**2*gamma0q*Ls12mu**5*amp(0) - &
                CF**3*gamma0cusp**3*Ls12mu**6*amp(0) + &
                192*CF*gamma2cusp*im*Ls12mu*Pi*amp(0) - &
                192*beta0*CF*gamma1cusp*im*Ls12mu**2*Pi*amp(0) - &
                192*CF*gamma0q*gamma1cusp*im*Ls12mu**2*Pi*amp(0) - &
                192*CF*gamma0cusp*gamma1q*im*Ls12mu**2*Pi*amp(0) + &
                64*beta0**2*CF*gamma0cusp*im*Ls12mu**3*Pi*amp(0) + &
                192*beta0*CF*gamma0cusp*gamma0q*im*Ls12mu**3*Pi*amp(0) + &
                96*CF*gamma0cusp*gamma0q**2*im*Ls12mu**3*Pi*amp(0) - &
                96*CF**2*gamma0cusp*gamma1cusp*im*Ls12mu**3*Pi*amp(0) + &
                40*beta0*CF**2*gamma0cusp**2*im*Ls12mu**4*Pi*amp(0) + &
                48*CF**2*gamma0cusp**2*gamma0q*im*Ls12mu**4*Pi*amp(0) + &
                6*CF**3*gamma0cusp**3*im*Ls12mu**5*Pi*amp(0) + &
                96*CF**2*gamma0cusp*gamma1cusp*im**2*Ls12mu**2*Pi**2*amp(0) - &
                48*beta0*CF**2*gamma0cusp**2*im**2*Ls12mu**3*Pi**2*amp(0) - &
                48*CF**2*gamma0cusp**2*gamma0q*im**2*Ls12mu**3*Pi**2*amp(0) - &
                12*CF**3*gamma0cusp**3*im**2*Ls12mu**4*Pi**2*amp(0) + &
                8*CF**3*gamma0cusp**3*im**3*Ls12mu**3*Pi**3*amp(0) + &
                32*beta1*Ls12mu*(6*gamma0q*Ls12mu*amp(0) + &
                   CF*gamma0cusp*Ls12mu*(Ls12mu - 3*im*Pi)*amp(0) - 24*amp(1)) - &
                768*gamma1q*Ls12mu*amp(1) + 768*beta0**2*Ls12mu**2*amp(1) + &
                1152*beta0*gamma0q*Ls12mu**2*amp(1) + &
                384*gamma0q**2*Ls12mu**2*amp(1) - &
                192*CF*gamma1cusp*Ls12mu**2*amp(1) + &
                256*beta0*CF*gamma0cusp*Ls12mu**3*amp(1) + &
                192*CF*gamma0cusp*gamma0q*Ls12mu**3*amp(1) + &
                24*CF**2*gamma0cusp**2*Ls12mu**4*amp(1) + &
                384*CF*gamma1cusp*im*Ls12mu*Pi*amp(1) - &
                576*beta0*CF*gamma0cusp*im*Ls12mu**2*Pi*amp(1) - &
                384*CF*gamma0cusp*gamma0q*im*Ls12mu**2*Pi*amp(1) - &
                96*CF**2*gamma0cusp**2*im*Ls12mu**3*Pi*amp(1) + &
                96*CF**2*gamma0cusp**2*im**2*Ls12mu**2*Pi**2*amp(1) - &
                3072*beta0*Ls12mu*amp(2) - 1536*gamma0q*Ls12mu*amp(2) - &
                384*CF*gamma0cusp*Ls12mu**2*amp(2) + &
                768*CF*gamma0cusp*im*Ls12mu*Pi*amp(2) + 3072*amp(3))/3072._dp
        endif

        amp(2) = (3*CF**2*gamma0cusp**2*Ls12mu**4*amp(0) + &
          4*CF*gamma0cusp*Ls12mu**3*&
           (2*beta0 + 6*gamma0q - 3*CF*gamma0cusp*im*Pi)*amp(0) + &
          12*Ls12mu**2*((4*beta0*gamma0q + 4*gamma0q**2 - 2*CF*gamma1cusp - &
                2*beta0*CF*gamma0cusp*im*Pi - 4*CF*gamma0cusp*gamma0q*im*Pi + &
                CF**2*gamma0cusp**2*im**2*Pi**2)*amp(0) - 4*CF*gamma0cusp*amp(1)&
             ) - 48*Ls12mu*(2*gamma1q*amp(0) - CF*gamma1cusp*im*Pi*amp(0) + &
             4*beta0*amp(1) + 4*gamma0q*amp(1) - 2*CF*gamma0cusp*im*Pi*amp(1)) &
           + 384*amp(2))/384._dp

        amp(1) = amp(1) + 0.5_dp*amp(0)*(-gamma0q + 0.5_dp*im*pi*CF*gamma0cusp)*Ls12mu &
                    -0.125_dp*amp(0)*CF*gamma0cusp*Ls12mu**2

    end subroutine rge_evolve

end module
