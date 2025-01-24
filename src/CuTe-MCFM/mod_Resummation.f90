!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
    module qtResummation
        use types
        use constants
        use iso_fortran_env
        use qtResummation_params
        implicit none

        public :: resummation
        public :: qtsubtraction
        public :: recoilBoost
        public :: gen_res

    private

    integer, save :: nf, nfmax
!$omp threadprivate(nf, nfmax)

    real(dp), save :: d2,d3
!$omp threadprivate(d2,d3)

    real(dp), save :: beta0,beta1,beta2
!$omp threadprivate(beta0,beta1,beta2)

    real(dp), save :: Gamma0, Gamma1,Gamma2
!$omp threadprivate(Gamma0,Gamma1,Gamma2)

    real(dp), save :: gammaq0, gammaq1, gammag0, gammag1, gammaq2, gammag2
!$omp threadprivate(gammaq0,gammaq1,gammag0,gammag1,gammaq2,gammag2)

    !! shared parameters with Fourier integrand
    integer, save :: order, nn
!$omp threadprivate(order,nn)
    real(dp), save :: qt, q2, mu, alphasMu
!$omp threadprivate(qt,q2,mu,alphasMu)

    ! parameter that defines whether this module is set up for
    ! a quark initiated process or gluon initiated
    logical, save :: initQuark = .true.
!$omp threadprivate(initQuark)

    contains

    subroutine update_nf_parameters(nf_in)
        implicit none
        integer, intent(in) :: nf_in

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

            beta0 = 11._dp/3._dp*CA - 4._dp/3._dp*TF*nf
            beta1 = 34._dp/3._dp*CA**2 - 20._dp/3._dp*CA*TF*nf - 4._dp*CF*TF*nf
            beta2 = 2857._dp/54._dp*CA**3 + nf*(CF**2 - 205._dp/18._dp*CF*CA - &
                    1415._dp/54._dp*CA**2) + nf**2*(11._dp/9._dp*CF + 79._dp/54._dp*CA)

            Gamma0 = 4._dp
            Gamma1 = 4._dp*((67._dp/9._dp - pi**2/3._dp)*CA - 20._dp/9._dp*TF*nf)
            Gamma2 = 16._dp*((245._dp/24._dp - 67._dp/54._dp*pi**2 + 11._dp*pi**4/180._dp &
                    + 11._dp/6._dp*zeta3)*CA**2 - (209._dp/108._dp - 5._dp*pi**2/27._dp &
                    + 7._dp/3._dp*zeta3)*CA*nf - (55._dp/24._dp - 2._dp*zeta3)*CF*nf - nf**2/27._dp)

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

    subroutine recoilBoost(pt,phi,p)
        implicit none
        include 'src/Inc/mxpart.f'
        include 'src/Inc/npart.f'
        real(dp), intent(in) :: pt,phi
        real(dp), intent(inout) :: p(mxpart,4)

        real(dp) :: qtx, qty
        real(dp) :: Lxx, Lxy, Lyy
        real(dp) :: energy, rt, btx, bty, gam
        integer :: j

        real(dp) :: boostxy(4,4)

        energy = - p(1,4) - p(2,4)
        qtx = pt*cos(phi)
        qty = pt*sin(phi)
        rt = sqrt(energy**2 + pt**2)
        btx = qtx/rt
        bty = qty/rt
        gam = 1._dp/sqrt(1._dp-btx**2-bty**2)

        Lxx = (gam-1._dp)*btx**2/(btx**2 + bty**2) + 1._dp
        Lyy = (gam-1._dp)*bty**2/(btx**2 + bty**2) + 1._dp
        Lxy = (gam-1._dp)*btx*bty/(btx**2 + bty**2)

        boostxy = reshape([ Lxx, Lxy, 0._dp, -btx*gam, &
                    Lxy, Lyy, 0._dp, -bty*gam, &
                    0._dp, 0._dp, 1._dp, 0._dp, &
                    -btx*gam, -bty*gam, 0._dp, gam ],[4,4])

        p(1,:) = matmul(boostxy, p(1,:))
        p(2,:) = matmul(boostxy, p(2,:))

        do j=3,npart+2
            p(j,:) = matmul(boostxy, p(j,:))
        enddo

    end subroutine

    function qstar(x)
        use LHAPDF, only: getalphas
        implicit none
        real(dp), intent(in) :: x
        real(dp) :: qstar

        real(dp) :: Ci

        if (initQuark) then
            Ci = CF
        else
            Ci = CA
        endif

        qstar = x**2 - q2 * exp(-pi/Ci/getalphas(x))

    end function

    subroutine qtsubtraction(p, q2, x1, x2, f1, f2, hard, muH, &
            order, muMult, resarray)
        use LHAPDF
        use Beamfunctions3L
        use SCET
        implicit none
        include 'src/Inc/mxpart.f'
        include 'src/Inc/kpart.f'
        include 'src/Inc/taucut.f'
        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(in) :: q2, x1, x2, muH
        integer, intent(in) :: f1, f2
        real(dp), intent(in) :: hard(0:3)
        integer, intent(in) :: order
        real(dp), intent(in), optional :: muMult
        real(dp), intent(out), allocatable :: resarray(:)

        include 'beamtype.f'

        logical :: initQuark
        real(dp) :: alphasMuH, beamcache(2,10)
        integer :: j

        real(dp) :: ci,gamma1rx,gamma2rx,gammai0,gammai1,gammai2
        real(dp) :: musq, logQsqDivMusq
        real(dp), allocatable :: logMusqDivCutsq(:), taucuts(:)
        real(dp), allocatable :: outpiece(:)

        allocate(resarray(size(tcutarray)+1))
        allocate(outpiece(size(tcutarray)+1))
        allocate(logMusqDivCutsq(size(tcutarray)+1))
        allocate(taucuts(size(tcutarray)+1))

        taucuts(1) = taucut
        taucuts(2:size(tcutarray)+1) = tcutarray(:)

        if (dynamictau) then
           do j=1,size(taucuts)
               taucuts(j) = getdynamictau(p, taucuts(j))
           enddo
        endif

        alphasMuH = getalphas(muH)

        beamcache(:,:) = 0._dp

        beamcache(1,1) = getbeam(ih1,f1,0,0,x1,muH,1)
        beamcache(2,1) = getbeam(ih2,f2,0,0,x2,muH,2)

        if (order >= 1) then
            beamcache(1,2) = getbeam(ih1,f1,1,0,x1,muH,1)
            beamcache(2,2) = getbeam(ih2,f2,1,0,x2,muH,2)

            beamcache(1,3) = getbeam(ih1,f1,1,1,x1,muH,1)
            beamcache(2,3) = getbeam(ih2,f2,1,1,x2,muH,2)

        endif

        if (order >= 2) then
            beamcache(1,4) = getbeam(ih1,f1,2,0,x1,muH,1)
            beamcache(2,4) = getbeam(ih2,f2,2,0,x2,muH,2)

            beamcache(1,5) = getbeam(ih1,f1,2,1,x1,muH,1)
            beamcache(2,5) = getbeam(ih2,f2,2,1,x2,muH,2)

            beamcache(1,6) = getbeam(ih1,f1,2,2,x1,muH,1)
            beamcache(2,6) = getbeam(ih2,f2,2,2,x2,muH,2)
        endif

        if (order >= 3) then
            beamcache(1,7) = getbeam(ih1,f1,3,0,x1,muH,1)
            beamcache(2,7) = getbeam(ih2,f2,3,0,x2,muH,2)

            beamcache(1,8) = getbeam(ih1,f1,3,1,x1,muH,1)
            beamcache(2,8) = getbeam(ih2,f2,3,1,x2,muH,2)

            beamcache(1,9) = getbeam(ih1,f1,3,2,x1,muH,1)
            beamcache(2,9) = getbeam(ih2,f2,3,2,x2,muH,2)

            beamcache(1,10) = getbeam(ih1,f1,3,3,x1,muH,1)
            beamcache(2,10) = getbeam(ih2,f2,3,3,x2,muH,2)
        endif

        if (f1 == 0 .and. f2 == 0) then
            initQuark = .false.
        else
            initQuark = .true.
        endif

        musq = muH**2
        call update_nf_parameters(getnumflavors(muH))
        logQsqDivMusq = log(q2/musq)
        logMusqDivCutsq = log(musq/taucuts**2)

        if (initQuark) then
            Ci = CF
            gammai0 = gammaq0
            gammai1 = gammaq1
            gammai2 = gammaq2
            gamma1rx = -d2*CF
            gamma2rx = -d3*CF
        else
            Ci = CA
            gammai0 = gammag0
            gammai1 = gammag1
            gammai2 = gammag2
            gamma1rx = -d2*CA
            gamma2rx = -d3*CA
        endif

         if (coeffonly) then
             resarray(:) = 0._dp
         else
             resarray(:) = beamcache(1,1)*beamcache(2,1)
         endif

        if ((coeffonly .and. order == 1) .or. (order >= 1 .and. (.not. coeffonly))) then
        !if (order == 1) then
            outpiece(:) =  -(Ci*Gamma0*LogMusqDivCutSq**2*beamcache(1,1)*beamcache(2,1))/2._dp + &
            beamcache(1,2)*beamcache(2,1) + beamcache(1,1)*beamcache(2,2) +  &
            LogMusqDivCutSq*(-2*gammai0*beamcache(1,1)*beamcache(2,1) -  &
               Ci*Gamma0*LogQsqDivMusq*beamcache(1,1)*beamcache(2,1) +  &
               beamcache(1,3)*beamcache(2,1) + beamcache(1,1)*beamcache(2,3)) + &
            beamcache(1,1)*beamcache(2,1)*hard(1)

            outpiece(:) = outpiece(:) * (alphasMuH/4._dp/pi)

            resarray = resarray + outpiece
        endif

        if ((coeffonly .and. order == 2) .or. ( order >= 2 .and. (.not. coeffonly))) then
        !if (order == 2) then
            outpiece(:) = (Ci**2*Gamma0**2*LogMusqDivCutSq**4*beamcache(1,1)*beamcache(2,1))/8._dp + &
            (4*beta0*Ci*Gamma0*Zeta3*beamcache(1,1)*beamcache(2,1))/3._dp - &
            4*Ci*Gamma0*gammai0*Zeta3*beamcache(1,1)*beamcache(2,1) + &
            2*Ci*Gamma0*Zeta3*beamcache(1,3)*beamcache(2,1) + &
            beamcache(1,4)*beamcache(2,1) + &
            LogQsqDivMusq*(gamma1rX*beamcache(1,1)*beamcache(2,1) - &
               2*Ci**2*Gamma0**2*Zeta3*beamcache(1,1)*beamcache(2,1)) + &
            beamcache(1,2)*beamcache(2,2) + &
            2*Ci*Gamma0*Zeta3*beamcache(1,1)*beamcache(2,3) + &
            LogMusqDivCutSq**3*(-(beta0*Ci*Gamma0*beamcache(1,1)*beamcache(2,1))/&
                3._dp + Ci*Gamma0*gammai0*beamcache(1,1)*beamcache(2,1) + &
               (Ci**2*Gamma0**2*LogQsqDivMusq*beamcache(1,1)*beamcache(2,1))/2._dp - &
               (Ci*Gamma0*beamcache(1,3)*beamcache(2,1))/2._dp - &
               (Ci*Gamma0*beamcache(1,1)*beamcache(2,3))/2._dp) + &
            beamcache(1,1)*beamcache(2,4) + beamcache(1,2)*beamcache(2,1)*hard(1) + &
            beamcache(1,1)*beamcache(2,2)*hard(1) + &
            LogMusqDivCutSq**2*(-(Ci*Gamma1*beamcache(1,1)*beamcache(2,1))/2._dp - &
               beta0*gammai0*beamcache(1,1)*beamcache(2,1) + &
               2*gammai0**2*beamcache(1,1)*beamcache(2,1) + &
               (Ci**2*Gamma0**2*LogQsqDivMusq**2*beamcache(1,1)*beamcache(2,1))/&
                2._dp - (Ci*Gamma0*beamcache(1,2)*beamcache(2,1))/2._dp - &
               2*gammai0*beamcache(1,3)*beamcache(2,1) + &
               beamcache(1,6)*beamcache(2,1) - &
               (Ci*Gamma0*beamcache(1,1)*beamcache(2,2))/2._dp - &
               2*gammai0*beamcache(1,1)*beamcache(2,3) + &
               beamcache(1,3)*beamcache(2,3) + &
               LogQsqDivMusq*(-(beta0*Ci*Gamma0*beamcache(1,1)*beamcache(2,1))/2._dp + &
                  2*Ci*Gamma0*gammai0*beamcache(1,1)*beamcache(2,1) - &
                  Ci*Gamma0*beamcache(1,3)*beamcache(2,1) - &
                  Ci*Gamma0*beamcache(1,1)*beamcache(2,3)) + &
               beamcache(1,1)*beamcache(2,6) - &
               (Ci*Gamma0*beamcache(1,1)*beamcache(2,1)*hard(1))/2._dp) + &
            LogMusqDivCutSq*(gamma1rX*beamcache(1,1)*beamcache(2,1) - &
               2*gammai1*beamcache(1,1)*beamcache(2,1) - &
               2*Ci**2*Gamma0**2*Zeta3*beamcache(1,1)*beamcache(2,1) - &
               2*gammai0*beamcache(1,2)*beamcache(2,1) + &
               beamcache(1,5)*beamcache(2,1) - &
               2*gammai0*beamcache(1,1)*beamcache(2,2) + &
               beamcache(1,3)*beamcache(2,2) + beamcache(1,2)*beamcache(2,3) + &
               beamcache(1,1)*beamcache(2,5) - &
               2*gammai0*beamcache(1,1)*beamcache(2,1)*hard(1) + &
               beamcache(1,3)*beamcache(2,1)*hard(1) + &
               beamcache(1,1)*beamcache(2,3)*hard(1) + &
               LogQsqDivMusq*(-(Ci*Gamma1*beamcache(1,1)*beamcache(2,1)) - &
                  Ci*Gamma0*beamcache(1,2)*beamcache(2,1) - &
                  Ci*Gamma0*beamcache(1,1)*beamcache(2,2) - &
                  Ci*Gamma0*beamcache(1,1)*beamcache(2,1)*hard(1))) + &
            beamcache(1,1)*beamcache(2,1)*hard(2)

            outpiece(:) = outpiece(:) * (alphasMuH/4._dp/pi)**2

            resarray = resarray + outpiece
        endif

        if ((coeffonly .and. order == 3) .or. (order >= 3 .and. (.not. coeffonly))) then
        !if (order == 3) then
            outpiece(:) = -(Ci**3*Gamma0**3*LogMusqDivCutSq**6*beamcache(1,1)*beamcache(2,1))/48._dp + &
                (8*beta0*Ci*Gamma1*Zeta3*beamcache(1,1)*beamcache(2,1))/3._dp + &
                2*Ci*Gamma0*gamma1rX*Zeta3*beamcache(1,1)*beamcache(2,1) + &
                (8*beta0**2*gammai0*Zeta3*beamcache(1,1)*beamcache(2,1))/3._dp - &
                4*Ci*Gamma1*gammai0*Zeta3*beamcache(1,1)*beamcache(2,1) - &
                8*beta0*gammai0**2*Zeta3*beamcache(1,1)*beamcache(2,1) + &
                (16*gammai0**3*Zeta3*beamcache(1,1)*beamcache(2,1))/3._dp - &
                4*Ci*Gamma0*gammai1*Zeta3*beamcache(1,1)*beamcache(2,1) + &
                (2*Ci**3*Gamma0**3*LogQsqDivMusq**3*Zeta3*beamcache(1,1)*beamcache(2,1))/3._dp - &
                (10*Ci**3*Gamma0**3*Zeta3**2*beamcache(1,1)*beamcache(2,1))/3._dp - &
                8*beta0*Ci**2*Gamma0**2*Zeta5*beamcache(1,1)*beamcache(2,1) + &
                12*Ci**2*Gamma0**2*gammai0*Zeta5*beamcache(1,1)*beamcache(2,1) + &
                (4*Ci*Gamma0*Zeta3*beta1*beamcache(1,1)*beamcache(2,1))/3._dp + &
                (4*beta0*Ci*Gamma0*Zeta3*beamcache(1,2)*beamcache(2,1))/3._dp - &
                4*Ci*Gamma0*gammai0*Zeta3*beamcache(1,2)*beamcache(2,1) + &
                2*Ci*Gamma1*Zeta3*beamcache(1,3)*beamcache(2,1) + &
                4*beta0*gammai0*Zeta3*beamcache(1,3)*beamcache(2,1) - &
                8*gammai0**2*Zeta3*beamcache(1,3)*beamcache(2,1) - &
                6*Ci**2*Gamma0**2*Zeta5*beamcache(1,3)*beamcache(2,1) + &
                2*Ci*Gamma0*Zeta3*beamcache(1,5)*beamcache(2,1) + &
                8*gammai0*Zeta3*beamcache(1,6)*beamcache(2,1) + beamcache(1,7)*beamcache(2,1) - &
                4*Zeta3*beamcache(1,10)*beamcache(2,1) + &
                (4*beta0*Ci*Gamma0*Zeta3*beamcache(1,1)*beamcache(2,2))/3._dp - &
                4*Ci*Gamma0*gammai0*Zeta3*beamcache(1,1)*beamcache(2,2) + &
                2*Ci*Gamma0*Zeta3*beamcache(1,3)*beamcache(2,2) + beamcache(1,4)*beamcache(2,2) + &
                2*Ci*Gamma1*Zeta3*beamcache(1,1)*beamcache(2,3) + &
                4*beta0*gammai0*Zeta3*beamcache(1,1)*beamcache(2,3) - &
                8*gammai0**2*Zeta3*beamcache(1,1)*beamcache(2,3) - &
                6*Ci**2*Gamma0**2*Zeta5*beamcache(1,1)*beamcache(2,3) + &
                2*Ci*Gamma0*Zeta3*beamcache(1,2)*beamcache(2,3) + &
                8*gammai0*Zeta3*beamcache(1,3)*beamcache(2,3) - 4*Zeta3*beamcache(1,6)*beamcache(2,3) + &
                LogMusqDivCutSq**5*((beta0*Ci**2*Gamma0**2*beamcache(1,1)*beamcache(2,1))/6._dp - &
                   (Ci**2*Gamma0**2*gammai0*beamcache(1,1)*beamcache(2,1))/4._dp - &
                   (Ci**3*Gamma0**3*LogQsqDivMusq*beamcache(1,1)*beamcache(2,1))/8._dp + &
                   (Ci**2*Gamma0**2*beamcache(1,3)*beamcache(2,1))/8._dp + &
                   (Ci**2*Gamma0**2*beamcache(1,1)*beamcache(2,3))/8._dp) + &
                LogQsqDivMusq**2*(-2*beta0*Ci**2*Gamma0**2*Zeta3*beamcache(1,1)*beamcache(2,1) + &
                   4*Ci**2*Gamma0**2*gammai0*Zeta3*beamcache(1,1)*beamcache(2,1) - &
                   2*Ci**2*Gamma0**2*Zeta3*beamcache(1,3)*beamcache(2,1) - &
                   2*Ci**2*Gamma0**2*Zeta3*beamcache(1,1)*beamcache(2,3)) + &
                beamcache(1,2)*beamcache(2,4) + 2*Ci*Gamma0*Zeta3*beamcache(1,1)*beamcache(2,5) + &
                8*gammai0*Zeta3*beamcache(1,1)*beamcache(2,6) - 4*Zeta3*beamcache(1,3)*beamcache(2,6) + &
                beamcache(1,1)*beamcache(2,7) - 4*Zeta3*beamcache(1,1)*beamcache(2,10) + &
                (4*beta0*Ci*Gamma0*Zeta3*beamcache(1,1)*beamcache(2,1)*hard(1))/3._dp - &
                4*Ci*Gamma0*gammai0*Zeta3*beamcache(1,1)*beamcache(2,1)*hard(1) + &
                2*Ci*Gamma0*Zeta3*beamcache(1,3)*beamcache(2,1)*hard(1) + &
                beamcache(1,4)*beamcache(2,1)*hard(1) + beamcache(1,2)*beamcache(2,2)*hard(1) + &
                2*Ci*Gamma0*Zeta3*beamcache(1,1)*beamcache(2,3)*hard(1) + &
                beamcache(1,1)*beamcache(2,4)*hard(1) + &
                LogMusqDivCutSq**4*(-(beta0**2*Ci*Gamma0*beamcache(1,1)*beamcache(2,1))/4._dp + &
                   (Ci**2*Gamma0*Gamma1*beamcache(1,1)*beamcache(2,1))/4._dp + &
                   (7*beta0*Ci*Gamma0*gammai0*beamcache(1,1)*beamcache(2,1))/6._dp - &
                   Ci*Gamma0*gammai0**2*beamcache(1,1)*beamcache(2,1) - &
                   (Ci**3*Gamma0**3*LogQsqDivMusq**2*beamcache(1,1)*beamcache(2,1))/4._dp + &
                   (Ci**2*Gamma0**2*beamcache(1,2)*beamcache(2,1))/8._dp - &
                   (beta0*Ci*Gamma0*beamcache(1,3)*beamcache(2,1))/3._dp + &
                   Ci*Gamma0*gammai0*beamcache(1,3)*beamcache(2,1) - &
                   (Ci*Gamma0*beamcache(1,6)*beamcache(2,1))/2._dp + &
                   (Ci**2*Gamma0**2*beamcache(1,1)*beamcache(2,2))/8._dp - &
                   (beta0*Ci*Gamma0*beamcache(1,1)*beamcache(2,3))/3._dp + &
                   Ci*Gamma0*gammai0*beamcache(1,1)*beamcache(2,3) - &
                   (Ci*Gamma0*beamcache(1,3)*beamcache(2,3))/2._dp + &
                   LogQsqDivMusq*((7*beta0*Ci**2*Gamma0**2*beamcache(1,1)*beamcache(2,1))/12._dp - &
                      Ci**2*Gamma0**2*gammai0*beamcache(1,1)*beamcache(2,1) + &
                      (Ci**2*Gamma0**2*beamcache(1,3)*beamcache(2,1))/2._dp + &
                      (Ci**2*Gamma0**2*beamcache(1,1)*beamcache(2,3))/2._dp) - &
                   (Ci*Gamma0*beamcache(1,1)*beamcache(2,6))/2._dp + &
                   (Ci**2*Gamma0**2*beamcache(1,1)*beamcache(2,1)*hard(1))/8._dp) + &
                LogQsqDivMusq*(gamma2rX*beamcache(1,1)*beamcache(2,1) + &
                   (4*beta0**2*Ci*Gamma0*Zeta3*beamcache(1,1)*beamcache(2,1))/3._dp - &
                   4*Ci**2*Gamma0*Gamma1*Zeta3*beamcache(1,1)*beamcache(2,1) - &
                   8*beta0*Ci*Gamma0*gammai0*Zeta3*beamcache(1,1)*beamcache(2,1) + &
                   8*Ci*Gamma0*gammai0**2*Zeta3*beamcache(1,1)*beamcache(2,1) + &
                   6*Ci**3*Gamma0**3*Zeta5*beamcache(1,1)*beamcache(2,1) + &
                   gamma1rX*beamcache(1,2)*beamcache(2,1) - &
                   2*Ci**2*Gamma0**2*Zeta3*beamcache(1,2)*beamcache(2,1) + &
                   2*beta0*Ci*Gamma0*Zeta3*beamcache(1,3)*beamcache(2,1) - &
                   8*Ci*Gamma0*gammai0*Zeta3*beamcache(1,3)*beamcache(2,1) + &
                   4*Ci*Gamma0*Zeta3*beamcache(1,6)*beamcache(2,1) + &
                   gamma1rX*beamcache(1,1)*beamcache(2,2) - &
                   2*Ci**2*Gamma0**2*Zeta3*beamcache(1,1)*beamcache(2,2) + &
                   2*beta0*Ci*Gamma0*Zeta3*beamcache(1,1)*beamcache(2,3) - &
                   8*Ci*Gamma0*gammai0*Zeta3*beamcache(1,1)*beamcache(2,3) + &
                   4*Ci*Gamma0*Zeta3*beamcache(1,3)*beamcache(2,3) + &
                   4*Ci*Gamma0*Zeta3*beamcache(1,1)*beamcache(2,6) + &
                   gamma1rX*beamcache(1,1)*beamcache(2,1)*hard(1) - &
                   2*Ci**2*Gamma0**2*Zeta3*beamcache(1,1)*beamcache(2,1)*hard(1)) + &
                LogMusqDivCutSq**3*((-2*beta0*Ci*Gamma1*beamcache(1,1)*beamcache(2,1))/3._dp - &
                   (Ci*Gamma0*gamma1rX*beamcache(1,1)*beamcache(2,1))/2._dp - &
                   (2*beta0**2*gammai0*beamcache(1,1)*beamcache(2,1))/3._dp + &
                   Ci*Gamma1*gammai0*beamcache(1,1)*beamcache(2,1) + &
                   2*beta0*gammai0**2*beamcache(1,1)*beamcache(2,1) - &
                   (4*gammai0**3*beamcache(1,1)*beamcache(2,1))/3._dp + &
                   Ci*Gamma0*gammai1*beamcache(1,1)*beamcache(2,1) - &
                   (Ci**3*Gamma0**3*LogQsqDivMusq**3*beamcache(1,1)*beamcache(2,1))/6._dp + &
                   (5*Ci**3*Gamma0**3*Zeta3*beamcache(1,1)*beamcache(2,1))/3._dp - &
                   (Ci*Gamma0*Beta1*beamcache(1,1)*beamcache(2,1))/3._dp - &
                   (beta0*Ci*Gamma0*beamcache(1,2)*beamcache(2,1))/3._dp + &
                   Ci*Gamma0*gammai0*beamcache(1,2)*beamcache(2,1) - &
                   (Ci*Gamma1*beamcache(1,3)*beamcache(2,1))/2._dp - &
                   beta0*gammai0*beamcache(1,3)*beamcache(2,1) + &
                   2*gammai0**2*beamcache(1,3)*beamcache(2,1) - &
                   (Ci*Gamma0*beamcache(1,5)*beamcache(2,1))/2._dp - &
                   2*gammai0*beamcache(1,6)*beamcache(2,1) + beamcache(1,10)*beamcache(2,1) - &
                   (beta0*Ci*Gamma0*beamcache(1,1)*beamcache(2,2))/3._dp + &
                   Ci*Gamma0*gammai0*beamcache(1,1)*beamcache(2,2) - &
                   (Ci*Gamma0*beamcache(1,3)*beamcache(2,2))/2._dp - &
                   (Ci*Gamma1*beamcache(1,1)*beamcache(2,3))/2._dp - &
                   beta0*gammai0*beamcache(1,1)*beamcache(2,3) + &
                   2*gammai0**2*beamcache(1,1)*beamcache(2,3) - &
                   (Ci*Gamma0*beamcache(1,2)*beamcache(2,3))/2._dp - &
                   2*gammai0*beamcache(1,3)*beamcache(2,3) + beamcache(1,6)*beamcache(2,3) + &
                   LogQsqDivMusq**2*((beta0*Ci**2*Gamma0**2*beamcache(1,1)*beamcache(2,1))/2._dp - &
                      Ci**2*Gamma0**2*gammai0*beamcache(1,1)*beamcache(2,1) + &
                      (Ci**2*Gamma0**2*beamcache(1,3)*beamcache(2,1))/2._dp + &
                      (Ci**2*Gamma0**2*beamcache(1,1)*beamcache(2,3))/2._dp) - &
                   (Ci*Gamma0*beamcache(1,1)*beamcache(2,5))/2._dp - &
                   2*gammai0*beamcache(1,1)*beamcache(2,6) + beamcache(1,3)*beamcache(2,6) + &
                   beamcache(1,1)*beamcache(2,10) - &
                   (beta0*Ci*Gamma0*beamcache(1,1)*beamcache(2,1)*hard(1))/3._dp + &
                   Ci*Gamma0*gammai0*beamcache(1,1)*beamcache(2,1)*hard(1) - &
                   (Ci*Gamma0*beamcache(1,3)*beamcache(2,1)*hard(1))/2._dp - &
                   (Ci*Gamma0*beamcache(1,1)*beamcache(2,3)*hard(1))/2._dp + &
                   LogQsqDivMusq*(-(beta0**2*Ci*Gamma0*beamcache(1,1)*beamcache(2,1))/3._dp + &
                      Ci**2*Gamma0*Gamma1*beamcache(1,1)*beamcache(2,1) + &
                      2*beta0*Ci*Gamma0*gammai0*beamcache(1,1)*beamcache(2,1) - &
                      2*Ci*Gamma0*gammai0**2*beamcache(1,1)*beamcache(2,1) + &
                      (Ci**2*Gamma0**2*beamcache(1,2)*beamcache(2,1))/2._dp - &
                      (beta0*Ci*Gamma0*beamcache(1,3)*beamcache(2,1))/2._dp + &
                      2*Ci*Gamma0*gammai0*beamcache(1,3)*beamcache(2,1) - &
                      Ci*Gamma0*beamcache(1,6)*beamcache(2,1) + &
                      (Ci**2*Gamma0**2*beamcache(1,1)*beamcache(2,2))/2._dp - &
                      (beta0*Ci*Gamma0*beamcache(1,1)*beamcache(2,3))/2._dp + &
                      2*Ci*Gamma0*gammai0*beamcache(1,1)*beamcache(2,3) - &
                      Ci*Gamma0*beamcache(1,3)*beamcache(2,3) - &
                      Ci*Gamma0*beamcache(1,1)*beamcache(2,6) + &
                      (Ci**2*Gamma0**2*beamcache(1,1)*beamcache(2,1)*hard(1))/2._dp)) + &
                beamcache(1,2)*beamcache(2,1)*hard(2) + beamcache(1,1)*beamcache(2,2)*hard(2) + &
                LogMusqDivCutSq**2*(2*beta0*gamma1rX*beamcache(1,1)*beamcache(2,1) - &
                   (Gamma2*Ci*beamcache(1,1)*beamcache(2,1))/2._dp - &
                   2*gamma1rX*gammai0*beamcache(1,1)*beamcache(2,1) - &
                   2*beta0*gammai1*beamcache(1,1)*beamcache(2,1) + &
                   4*gammai0*gammai1*beamcache(1,1)*beamcache(2,1) - &
                   (20*beta0*Ci**2*Gamma0**2*Zeta3*beamcache(1,1)*beamcache(2,1))/3._dp + &
                   10*Ci**2*Gamma0**2*gammai0*Zeta3*beamcache(1,1)*beamcache(2,1) - &
                   gammai0*Beta1*beamcache(1,1)*beamcache(2,1) - &
                   (Ci*Gamma1*beamcache(1,2)*beamcache(2,1))/2._dp - &
                   beta0*gammai0*beamcache(1,2)*beamcache(2,1) + &
                   2*gammai0**2*beamcache(1,2)*beamcache(2,1) + gamma1rX*beamcache(1,3)*beamcache(2,1) - &
                   2*gammai1*beamcache(1,3)*beamcache(2,1) - &
                   5*Ci**2*Gamma0**2*Zeta3*beamcache(1,3)*beamcache(2,1) - &
                   (Ci*Gamma0*beamcache(1,4)*beamcache(2,1))/2._dp - &
                   2*gammai0*beamcache(1,5)*beamcache(2,1) + beamcache(1,9)*beamcache(2,1) - &
                   (Ci*Gamma1*beamcache(1,1)*beamcache(2,2))/2._dp - &
                   beta0*gammai0*beamcache(1,1)*beamcache(2,2) + &
                   2*gammai0**2*beamcache(1,1)*beamcache(2,2) - &
                   (Ci*Gamma0*beamcache(1,2)*beamcache(2,2))/2._dp - &
                   2*gammai0*beamcache(1,3)*beamcache(2,2) + beamcache(1,6)*beamcache(2,2) + &
                   gamma1rX*beamcache(1,1)*beamcache(2,3) - 2*gammai1*beamcache(1,1)*beamcache(2,3) - &
                   5*Ci**2*Gamma0**2*Zeta3*beamcache(1,1)*beamcache(2,3) - &
                   2*gammai0*beamcache(1,2)*beamcache(2,3) + beamcache(1,5)*beamcache(2,3) - &
                   (Ci*Gamma0*beamcache(1,1)*beamcache(2,4))/2._dp - &
                   2*gammai0*beamcache(1,1)*beamcache(2,5) + beamcache(1,3)*beamcache(2,5) + &
                   beamcache(1,2)*beamcache(2,6) + beamcache(1,1)*beamcache(2,9) - &
                   (Ci*Gamma1*beamcache(1,1)*beamcache(2,1)*hard(1))/2._dp - &
                   beta0*gammai0*beamcache(1,1)*beamcache(2,1)*hard(1) + &
                   2*gammai0**2*beamcache(1,1)*beamcache(2,1)*hard(1) - &
                   (Ci*Gamma0*beamcache(1,2)*beamcache(2,1)*hard(1))/2._dp - &
                   2*gammai0*beamcache(1,3)*beamcache(2,1)*hard(1) + &
                   beamcache(1,6)*beamcache(2,1)*hard(1) - &
                   (Ci*Gamma0*beamcache(1,1)*beamcache(2,2)*hard(1))/2._dp - &
                   2*gammai0*beamcache(1,1)*beamcache(2,3)*hard(1) + &
                   beamcache(1,3)*beamcache(2,3)*hard(1) + beamcache(1,1)*beamcache(2,6)*hard(1) + &
                   LogQsqDivMusq**2*(Ci**2*Gamma0*Gamma1*beamcache(1,1)*beamcache(2,1) + &
                      (Ci**2*Gamma0**2*beamcache(1,2)*beamcache(2,1))/2._dp + &
                      (Ci**2*Gamma0**2*beamcache(1,1)*beamcache(2,2))/2._dp + &
                      (Ci**2*Gamma0**2*beamcache(1,1)*beamcache(2,1)*hard(1))/2._dp) + &
                   LogQsqDivMusq*(-(beta0*Ci*Gamma1*beamcache(1,1)*beamcache(2,1)) - &
                      (3*Ci*Gamma0*gamma1rX*beamcache(1,1)*beamcache(2,1))/2._dp + &
                      2*Ci*Gamma1*gammai0*beamcache(1,1)*beamcache(2,1) + &
                      2*Ci*Gamma0*gammai1*beamcache(1,1)*beamcache(2,1) + &
                      5*Ci**3*Gamma0**3*Zeta3*beamcache(1,1)*beamcache(2,1) - &
                      (Ci*Gamma0*Beta1*beamcache(1,1)*beamcache(2,1))/2._dp - &
                      (beta0*Ci*Gamma0*beamcache(1,2)*beamcache(2,1))/2._dp + &
                      2*Ci*Gamma0*gammai0*beamcache(1,2)*beamcache(2,1) - &
                      Ci*Gamma1*beamcache(1,3)*beamcache(2,1) - &
                      Ci*Gamma0*beamcache(1,5)*beamcache(2,1) - &
                      (beta0*Ci*Gamma0*beamcache(1,1)*beamcache(2,2))/2._dp + &
                      2*Ci*Gamma0*gammai0*beamcache(1,1)*beamcache(2,2) - &
                      Ci*Gamma0*beamcache(1,3)*beamcache(2,2) - &
                      Ci*Gamma1*beamcache(1,1)*beamcache(2,3) - &
                      Ci*Gamma0*beamcache(1,2)*beamcache(2,3) - &
                      Ci*Gamma0*beamcache(1,1)*beamcache(2,5) - &
                      (beta0*Ci*Gamma0*beamcache(1,1)*beamcache(2,1)*hard(1))/2._dp + &
                      2*Ci*Gamma0*gammai0*beamcache(1,1)*beamcache(2,1)*hard(1) - &
                      Ci*Gamma0*beamcache(1,3)*beamcache(2,1)*hard(1) - &
                      Ci*Gamma0*beamcache(1,1)*beamcache(2,3)*hard(1)) - &
                   (Ci*Gamma0*beamcache(1,1)*beamcache(2,1)*hard(2))/2._dp) + &
                LogMusqDivCutSq*(gamma2rX*beamcache(1,1)*beamcache(2,1) - &
                   2*gammai2*beamcache(1,1)*beamcache(2,1) + &
                   4*beta0**2*Ci*Gamma0*Zeta3*beamcache(1,1)*beamcache(2,1) - &
                   4*Ci**2*Gamma0*Gamma1*Zeta3*beamcache(1,1)*beamcache(2,1) - &
                   (56*beta0*Ci*Gamma0*gammai0*Zeta3*beamcache(1,1)*beamcache(2,1))/3._dp + &
                   16*Ci*Gamma0*gammai0**2*Zeta3*beamcache(1,1)*beamcache(2,1) + &
                   6*Ci**3*Gamma0**3*Zeta5*beamcache(1,1)*beamcache(2,1) + &
                   gamma1rX*beamcache(1,2)*beamcache(2,1) - 2*gammai1*beamcache(1,2)*beamcache(2,1) - &
                   2*Ci**2*Gamma0**2*Zeta3*beamcache(1,2)*beamcache(2,1) + &
                   (16*beta0*Ci*Gamma0*Zeta3*beamcache(1,3)*beamcache(2,1))/3._dp - &
                   16*Ci*Gamma0*gammai0*Zeta3*beamcache(1,3)*beamcache(2,1) - &
                   2*gammai0*beamcache(1,4)*beamcache(2,1) + &
                   8*Ci*Gamma0*Zeta3*beamcache(1,6)*beamcache(2,1) + beamcache(1,8)*beamcache(2,1) + &
                   LogQsqDivMusq**2*(-(Ci*Gamma0*gamma1rX*beamcache(1,1)*beamcache(2,1)) + &
                      4*Ci**3*Gamma0**3*Zeta3*beamcache(1,1)*beamcache(2,1)) + &
                   gamma1rX*beamcache(1,1)*beamcache(2,2) - 2*gammai1*beamcache(1,1)*beamcache(2,2) - &
                   2*Ci**2*Gamma0**2*Zeta3*beamcache(1,1)*beamcache(2,2) - &
                   2*gammai0*beamcache(1,2)*beamcache(2,2) + beamcache(1,5)*beamcache(2,2) + &
                   (16*beta0*Ci*Gamma0*Zeta3*beamcache(1,1)*beamcache(2,3))/3._dp - &
                   16*Ci*Gamma0*gammai0*Zeta3*beamcache(1,1)*beamcache(2,3) + &
                   8*Ci*Gamma0*Zeta3*beamcache(1,3)*beamcache(2,3) + beamcache(1,4)*beamcache(2,3) - &
                   2*gammai0*beamcache(1,1)*beamcache(2,4) + beamcache(1,3)*beamcache(2,4) + &
                   beamcache(1,2)*beamcache(2,5) + 8*Ci*Gamma0*Zeta3*beamcache(1,1)*beamcache(2,6) + &
                   beamcache(1,1)*beamcache(2,8) + gamma1rX*beamcache(1,1)*beamcache(2,1)*hard(1) - &
                   2*gammai1*beamcache(1,1)*beamcache(2,1)*hard(1) - &
                   2*Ci**2*Gamma0**2*Zeta3*beamcache(1,1)*beamcache(2,1)*hard(1) - &
                   2*gammai0*beamcache(1,2)*beamcache(2,1)*hard(1) + &
                   beamcache(1,5)*beamcache(2,1)*hard(1) - &
                   2*gammai0*beamcache(1,1)*beamcache(2,2)*hard(1) + &
                   beamcache(1,3)*beamcache(2,2)*hard(1) + beamcache(1,2)*beamcache(2,3)*hard(1) + &
                   beamcache(1,1)*beamcache(2,5)*hard(1) - &
                   2*gammai0*beamcache(1,1)*beamcache(2,1)*hard(2) + &
                   beamcache(1,3)*beamcache(2,1)*hard(2) + beamcache(1,1)*beamcache(2,3)*hard(2) + &
                   LogQsqDivMusq*(2*beta0*gamma1rX*beamcache(1,1)*beamcache(2,1) - &
                      Gamma2*Ci*beamcache(1,1)*beamcache(2,1) - &
                      2*gamma1rX*gammai0*beamcache(1,1)*beamcache(2,1) - &
                      (28*beta0*Ci**2*Gamma0**2*Zeta3*beamcache(1,1)*beamcache(2,1))/3._dp + &
                      16*Ci**2*Gamma0**2*gammai0*Zeta3*beamcache(1,1)*beamcache(2,1) - &
                      Ci*Gamma1*beamcache(1,2)*beamcache(2,1) + gamma1rX*beamcache(1,3)*beamcache(2,1) - &
                      8*Ci**2*Gamma0**2*Zeta3*beamcache(1,3)*beamcache(2,1) - &
                      Ci*Gamma0*beamcache(1,4)*beamcache(2,1) - &
                      Ci*Gamma1*beamcache(1,1)*beamcache(2,2) - &
                      Ci*Gamma0*beamcache(1,2)*beamcache(2,2) + gamma1rX*beamcache(1,1)*beamcache(2,3) - &
                      8*Ci**2*Gamma0**2*Zeta3*beamcache(1,1)*beamcache(2,3) - &
                      Ci*Gamma0*beamcache(1,1)*beamcache(2,4) - &
                      Ci*Gamma1*beamcache(1,1)*beamcache(2,1)*hard(1) - &
                      Ci*Gamma0*beamcache(1,2)*beamcache(2,1)*hard(1) - &
                      Ci*Gamma0*beamcache(1,1)*beamcache(2,2)*hard(1) - &
                      Ci*Gamma0*beamcache(1,1)*beamcache(2,1)*hard(2))) + &
                beamcache(1,1)*beamcache(2,1)*hard(3)

            outpiece(:) = outpiece(:) * (alphasMuH/4._dp/pi)**3

            resarray = resarray + outpiece
        endif

        resarray(:) = resarray(:) * hard(0)

    end subroutine

    subroutine resummation(q2_in, qt_in, x1, x2, f1, f2, hard1, muH, order_in, muMult, res, resexp, hard2)
        use LHAPDF
        use Beamfunctions3L
        use newton_raphson, only: find_root
        use qtResummationHardEvolution
        use qtResummationFourier
        use qtResummation_params, only : scalevar_rapidity, scalevar_rapidity_mult, scalevar_rapidity_i
        use MCFMSettings, only: resexp_linPC
        implicit none
        real(dp), intent(in) :: q2_in, qt_in, x1, x2, muH
        integer, intent(in) :: f1, f2
        real(dp), intent(in) :: hard1
        integer, intent(in) :: order_in
        real(dp), intent(in), optional :: muMult
        real(dp), intent(out), optional :: res, resexp
        real(dp), intent(in), optional :: hard2

        include 'facscale.f'
        include 'beamtype.f'

        real(dp) :: logmuRmuF
        real(dp) :: beamFunf1(0:3), beamFunf2(0:3)
        real(dp) :: beamFunf1_grid(0:3), beamFunf2_grid(0:3)
        real(dp) :: beamFunf1_exp(0:2), beamFunf2_exp(0:2)
        real(dp) :: Mi(0:6)
        real(dp) :: alphasMuH, alphasMuMatch
        real(dp) :: mub, muc
        real(dp) :: alphas4p, alphas4m, alphas3p, alphas3m
        real(dp) :: Ci, gammai0, gammai1, gammai2, gamma1rx, gamma2rx, LM, LMU

        real(dp) :: root, muStar
        logical :: ierr

        integer :: j
        real(dp) :: resnew, resnew2
        real(dp) :: resexp1, resexp2, resexp3

        real(dp) :: chk

        real(dp) :: logR

        real(dp) :: ason4pi, beamcache(2,16)

        if (present(res) .and. present(resexp)) then
            ! because of setting nf based on mu or muMatch
            write (*,*) "res and resexp simultaneously currently unsupported"
            error stop
        endif

        !real(dp) :: test

        order = order_in
        q2 = q2_in
        qt = qt_in

        if (f1 == 0 .and. f2 == 0) then
            initQuark = .false.
        else
            initQuark = .true.
        endif

        call find_root(qstar, 1._dp, 1.0e-3_dp, 100, root, ierr)
        if (ierr .eqv. .false.) then
            write (*,*) "WARNING: could not determine qstar, using default of 1.88"
            root = 1.88_dp
        endif
        muStar = root

        ! Fourier integral is numerically unstable below mu~2GeV
        if (present(muMult)) then
            mu = max(muMult*(qt + muStar*exp(-qt/muStar)), 2.0d0)
        else
            mu = max(qt + muStar*exp(-qt/muStar), 2.0d0)
        endif

        ! alpha at resummation scale, hard scale and matching scale
        alphasMu = getalphas(mu)
        alphasMuH = getalphas(muH)


        ! update all nf constants in this module, including nf itself
        if (present(resexp)) then
            ! to do: we evaluate the PDFs at facscale, but everything else at muH
            ! logmuRmuF = log(muH**2/facscale**2)
            alphasMuMatch = getalphas(muH)
            call update_nf_parameters(getnumflavors(muH))
        else
            call update_nf_parameters(getnumflavors(mu))
        endif

        ! get number of flavors at hard scale, this is used in the hard function evolution
        ! to evolve downwards
        nfmax = getnumflavors(muH)

        ! update nf for beamfunctions

        ! PDFs in LHAPDF typically don't expose flavor thresholds,
        ! so we use the masses instead.
        ! Since the PDFs also use grid interpolation for alphas, the
        ! thresholds will be missing or washed out, but the code below
        ! will allow us to roll our own alphas running.
        if (nf < 5) then
            mub = getquarkMass(5)
            alphas4p = getalphas(mub + 1d-6)
            alphas4m = getalphas(mub - 1d-6)
        endif
        if (nf < 4) then
            muc = getquarkMass(4)
            alphas3p = getalphas(muc + 1d-6)
            alphas3m = getalphas(muc - 1d-6)
        endif

        if (initQuark) then
            Ci = CF
            gammai0 = gammaq0
            gammai1 = gammaq1
            gammai2 = gammaq2
            gamma1rx = -d2*CF
            gamma2rx = -d3*CF
        else
            Ci = CA
            gammai0 = gammag0
            gammai1 = gammag1
            gammai2 = gammag2
            gamma1rx = -d2*CA
            gamma2rx = -d3*CA
        endif

        if (present(res)) then

        Mi(0) = fourierM(0, order, nf, initQuark, qt, q2, mu, alphasMu)
        if (order >= 3) then
            Mi(1) = fourierM(1, order, nf, initQuark, qt, q2, mu, alphasMu)
            Mi(2) = fourierM(2, order, nf, initQuark, qt, q2, mu, alphasMu)
        endif
        if (order >= 5) then
            Mi(3) = fourierM(3, order, nf, initQuark, qt, q2, mu, alphasMu)
            Mi(4) = fourierM(4, order, nf, initQuark, qt, q2, mu, alphasMu)
        endif
        if (order >= 7) then
            Mi(5) = fourierM(5, order, nf, initQuark, qt, q2, mu, alphasMu)
            Mi(6) = fourierM(6, order, nf, initQuark, qt, q2, mu, alphasMu)
        endif

    
        ason4pi = alphasMu/4._dp/pi        

        beamcache(:,:) = 0._dp

        beamcache(1,1) = getbeam(ih1,f1,0,0,x1,mu,1)
        beamcache(2,1) = getbeam(ih2,f2,0,0,x2,mu,2)

        if (order >= 3) then
            beamcache(1,2) = getbeam(ih1,f1,1,0,x1,mu,1)
            beamcache(2,2) = getbeam(ih2,f2,1,0,x2,mu,2)

            beamcache(1,3) = getbeam(ih1,f1,1,1,x1,mu,1)
            beamcache(2,3) = getbeam(ih2,f2,1,1,x2,mu,2)

            beamcache(1,6) = getbeam(ih1,f1,2,2,x1,mu,1)
            beamcache(2,6) = getbeam(ih2,f2,2,2,x2,mu,2)

        endif

        if (order >= 5) then
            beamcache(1,4) = getbeam(ih1,f1,2,0,x1,mu,1)
            beamcache(2,4) = getbeam(ih2,f2,2,0,x2,mu,2)

            beamcache(1,5) = getbeam(ih1,f1,2,1,x1,mu,1)
            beamcache(2,5) = getbeam(ih2,f2,2,1,x2,mu,2)

            beamcache(1,9) = getbeam(ih1,f1,3,2,x1,mu,1)
            beamcache(2,9) = getbeam(ih2,f2,3,2,x2,mu,2)

            beamcache(1,10) = getbeam(ih1,f1,3,3,x1,mu,1)
            beamcache(2,10) = getbeam(ih2,f2,3,3,x2,mu,2)

            ! g^4
            beamcache(1,11) = getbeam(ih1,f1,4,4,x1,mu,1)
            beamcache(2,11) = getbeam(ih2,f2,4,4,x2,mu,2)
        endif

        if (order >= 7) then
            beamcache(1,7) = getbeam(ih1,f1,3,0,x1,mu,1)
            beamcache(2,7) = getbeam(ih2,f2,3,0,x2,mu,2)

            beamcache(1,8) = getbeam(ih1,f1,3,1,x1,mu,1)
            beamcache(2,8) = getbeam(ih2,f2,3,1,x2,mu,2)


            ! g^5
            beamcache(1,12) = getbeam(ih1,f1,4,3,x1,mu,1)
            beamcache(2,12) = getbeam(ih2,f2,4,3,x2,mu,2)

            beamcache(1,13) = getbeam(ih1,f1,5,5,x1,mu,1)
            beamcache(2,13) = getbeam(ih2,f2,5,5,x2,mu,2)


            ! g^6
            beamcache(1,14) = getbeam(ih1,f1,4,2,x1,mu,1)
            beamcache(2,14) = getbeam(ih2,f2,4,2,x2,mu,2)

            beamcache(1,15) = getbeam(ih1,f1,5,4,x1,mu,1)
            beamcache(2,15) = getbeam(ih2,f2,5,4,x2,mu,2)

            beamcache(1,16) = getbeam(ih1,f1,6,6,x1,mu,1)
            beamcache(2,16) = getbeam(ih2,f2,6,6,x2,mu,2)
        endif

        resnew = 0._dp
        resnew = Mi(0) * beamcache(1,1)*beamcache(2,1)

        if (scalevar_rapidity_i > 0) then
            logR = log(scalevar_rapidity_mult(scalevar_rapidity_i))      
        else
            logR = 0._dp
        endif

        if (order >= 3) then
            ! g
            resnew = resnew + ason4pi*(Ci*Gamma0*logR*beamcache(1,1)*beamcache(2,1) + &
                beamcache(1,3)*beamcache(2,1) + beamcache(1,1)*beamcache(2,3))*Mi(1)
            ! g^2
            resnew = resnew + ason4pi*(beamcache(1,2)*beamcache(2,1) + &
                   beamcache(1,1)*beamcache(2,2))*Mi(0) + &
                (ason4pi**2*(Ci**2*Gamma0**2*logR**2*beamcache(1,1)*beamcache(2,1) + &
                     Ci*Gamma0*logR*(beta0*beamcache(1,1)*beamcache(2,1) + &
                        2*beamcache(1,3)*beamcache(2,1) + &
                        2*beamcache(1,1)*beamcache(2,3)) + &
                     2*(beamcache(1,6)*beamcache(2,1) + beamcache(1,3)*beamcache(2,3) + &
                        beamcache(1,1)*beamcache(2,6)))*Mi(2))/2._dp
        endif

        if (order >= 5) then
            ! g^3
            resnew = resnew + ason4pi**2*(Ci*Gamma1*logR*beamcache(1,1)*beamcache(2,1) + &
                   beamcache(1,5)*beamcache(2,1) + beamcache(1,3)*beamcache(2,2) + &
                   Ci*Gamma0*logR*(beamcache(1,2)*beamcache(2,1) + &
                      beamcache(1,1)*beamcache(2,2)) + beamcache(1,2)*beamcache(2,3) + &
                   beamcache(1,1)*beamcache(2,5))*Mi(1) + &
                (ason4pi**3*(Ci**3*Gamma0**3*logR**3*beamcache(1,1)*beamcache(2,1) + &
                     3*Ci**2*Gamma0**2*logR**2*&
                      (beta0*beamcache(1,1)*beamcache(2,1) + &
                        beamcache(1,3)*beamcache(2,1) + beamcache(1,1)*beamcache(2,3)) &
                      + Ci*Gamma0*logR*(2*beta0**2*beamcache(1,1)*beamcache(2,1) + &
                        3*beta0*(beamcache(1,3)*beamcache(2,1) + &
                           beamcache(1,1)*beamcache(2,3)) + &
                        6*(beamcache(1,6)*beamcache(2,1) + &
                           beamcache(1,3)*beamcache(2,3) + beamcache(1,1)*beamcache(2,6)&
                           )) + 6*(beamcache(1,10)*beamcache(2,1) + &
                        beamcache(1,6)*beamcache(2,3) + beamcache(1,3)*beamcache(2,6) + &
                        beamcache(1,1)*beamcache(2,10)))*Mi(3))/6._dp
            ! g^4
        resnew = resnew + ason4pi**2*(-(gamma1rX*logR*beamcache(1,1)*beamcache(2,1)) +  &
           beamcache(1,4)*beamcache(2,1) + beamcache(1,2)*beamcache(2,2) +  &
           beamcache(1,1)*beamcache(2,4))*Mi(0) +  &
        (ason4pi**3*(Ci**2*Gamma0**2*logR**2* &
              (beamcache(1,2)*beamcache(2,1) + beamcache(1,1)*beamcache(2,2)) +  &
             Ci*Gamma0*logR*(beta1*beamcache(1,1)*beamcache(2,1) +  &
                2*Ci*Gamma1*logR*beamcache(1,1)*beamcache(2,1) +  &
                beta0*beamcache(1,2)*beamcache(2,1) +  &
                2*beamcache(1,5)*beamcache(2,1) +  &
                beta0*beamcache(1,1)*beamcache(2,2) +  &
                2*beamcache(1,3)*beamcache(2,2) +  &
                2*beamcache(1,2)*beamcache(2,3) +  &
                2*beamcache(1,1)*beamcache(2,5)) +  &
             2*(beamcache(1,9)*beamcache(2,1) + beamcache(1,6)*beamcache(2,2) +  &
                beamcache(1,5)*beamcache(2,3) +  &
                Ci*Gamma1*logR*(beta0*beamcache(1,1)*beamcache(2,1) +  &
                   beamcache(1,3)*beamcache(2,1) + beamcache(1,1)*beamcache(2,3) &
                   ) + beamcache(1,3)*beamcache(2,5) +  &
                beamcache(1,2)*beamcache(2,6) + beamcache(1,1)*beamcache(2,9)))* &
           Mi(2))/2._dp + (ason4pi**4* &
           (Ci**4*Gamma0**4*logR**4*beamcache(1,1)*beamcache(2,1) +  &
             2*Ci**3*Gamma0**3*logR**3* &
              (3*beta0*beamcache(1,1)*beamcache(2,1) +  &
                2*beamcache(1,3)*beamcache(2,1) +  &
                2*beamcache(1,1)*beamcache(2,3)) +  &
             Ci**2*Gamma0**2*logR**2* &
              (11*beta0**2*beamcache(1,1)*beamcache(2,1) +  &
                12*beta0*(beamcache(1,3)*beamcache(2,1) +  &
                   beamcache(1,1)*beamcache(2,3)) +  &
                12*(beamcache(1,6)*beamcache(2,1) +  &
                   beamcache(1,3)*beamcache(2,3) + beamcache(1,1)*beamcache(2,6) &
                   )) + 2*Ci*Gamma0*logR* &
              (3*beta0**3*beamcache(1,1)*beamcache(2,1) +  &
                4*beta0**2*(beamcache(1,3)*beamcache(2,1) +  &
                   beamcache(1,1)*beamcache(2,3)) +  &
                6*beta0*(beamcache(1,6)*beamcache(2,1) +  &
                   beamcache(1,3)*beamcache(2,3) + beamcache(1,1)*beamcache(2,6) &
                   ) + 12*(beamcache(1,10)*beamcache(2,1) +  &
                   beamcache(1,6)*beamcache(2,3) +  &
                   beamcache(1,3)*beamcache(2,6) +  &
                   beamcache(1,1)*beamcache(2,10))) +  &
             24*(beamcache(1,11)*beamcache(2,1) +  &
                beamcache(1,10)*beamcache(2,3) +  &
                beamcache(1,6)*beamcache(2,6) +  &
                beamcache(1,3)*beamcache(2,10) + beamcache(1,1)*beamcache(2,11)) &
             )*Mi(4))/24._dp

                ! additional piece unique to gluon-gluon
                if (.not. initQuark) then
                    resnew = resnew + (alphasMu/4._dp/pi)**2*Mi(0)*( &
                        getbeam2(ih1,f1,1,0,x1,mu)*getbeam2(ih2,f2,1,0,x2,mu) )
                endif
        endif

        if (order >= 7) then
            resnew2 = resnew

        ! g^5
        resnew = resnew + ason4pi**3*(Ci*Gamma2*logR*beamcache(1,1)*beamcache(2,1) +  &
           Ci*Gamma1*logR*beamcache(1,2)*beamcache(2,1) +  &
           Ci*Gamma0*logR*beamcache(1,4)*beamcache(2,1) +  &
           beamcache(1,8)*beamcache(2,1) +  &
           Ci*Gamma1*logR*beamcache(1,1)*beamcache(2,2) +  &
           Ci*Gamma0*logR*beamcache(1,2)*beamcache(2,2) +  &
           beamcache(1,5)*beamcache(2,2) + beamcache(1,4)*beamcache(2,3) -  &
           gamma1rX*logR*(2*beta0*beamcache(1,1)*beamcache(2,1) +  &
              Ci*Gamma0*logR*beamcache(1,1)*beamcache(2,1) +  &
              beamcache(1,3)*beamcache(2,1) + beamcache(1,1)*beamcache(2,3)) +  &
           Ci*Gamma0*logR*beamcache(1,1)*beamcache(2,4) +  &
           beamcache(1,3)*beamcache(2,4) + beamcache(1,2)*beamcache(2,5) +  &
           beamcache(1,1)*beamcache(2,8))*Mi(1) +  &
        (ason4pi**4*(Ci**3*Gamma0**3*logR**3* &
              (beamcache(1,2)*beamcache(2,1) + beamcache(1,1)*beamcache(2,2)) +  &
             3*Ci**2*Gamma0**2*logR**2* &
              (beta1*beamcache(1,1)*beamcache(2,1) +  &
                Ci*Gamma1*logR*beamcache(1,1)*beamcache(2,1) +  &
                beta0*beamcache(1,2)*beamcache(2,1) +  &
                beamcache(1,5)*beamcache(2,1) +  &
                beta0*beamcache(1,1)*beamcache(2,2) +  &
                beamcache(1,3)*beamcache(2,2) + beamcache(1,2)*beamcache(2,3) +  &
                beamcache(1,1)*beamcache(2,5)) +  &
             Ci*Gamma0*logR*(5*beta0*beta1*beamcache(1,1)*beamcache(2,1) +  &
                2*beta0**2*(beamcache(1,2)*beamcache(2,1) +  &
                   beamcache(1,1)*beamcache(2,2)) +  &
                3*beta1*(beamcache(1,3)*beamcache(2,1) +  &
                   beamcache(1,1)*beamcache(2,3)) +  &
                3*Ci*Gamma1*logR* &
                 (3*beta0*beamcache(1,1)*beamcache(2,1) +  &
                   2*beamcache(1,3)*beamcache(2,1) +  &
                   2*beamcache(1,1)*beamcache(2,3)) +  &
                3*beta0*(beamcache(1,5)*beamcache(2,1) +  &
                   beamcache(1,3)*beamcache(2,2) +  &
                   beamcache(1,2)*beamcache(2,3) + beamcache(1,1)*beamcache(2,5) &
                   ) + 6*(beamcache(1,9)*beamcache(2,1) +  &
                   beamcache(1,6)*beamcache(2,2) +  &
                   beamcache(1,5)*beamcache(2,3) +  &
                   beamcache(1,3)*beamcache(2,5) +  &
                   beamcache(1,2)*beamcache(2,6) + beamcache(1,1)*beamcache(2,9) &
                   )) + 6*(beamcache(1,12)*beamcache(2,1) +  &
                beamcache(1,10)*beamcache(2,2) +  &
                beamcache(1,9)*beamcache(2,3) + beamcache(1,6)*beamcache(2,5) +  &
                beamcache(1,5)*beamcache(2,6) +  &
                Ci*Gamma1*logR*(beta0**2*beamcache(1,1)*beamcache(2,1) +  &
                   beta0*beamcache(1,3)*beamcache(2,1) +  &
                   beamcache(1,6)*beamcache(2,1) +  &
                   beta0*beamcache(1,1)*beamcache(2,3) +  &
                   beamcache(1,3)*beamcache(2,3) + beamcache(1,1)*beamcache(2,6) &
                   ) + beamcache(1,3)*beamcache(2,9) +  &
                beamcache(1,2)*beamcache(2,10) + beamcache(1,1)*beamcache(2,12)) &
             )*Mi(3))/6._dp + (ason4pi**5* &
           (Ci**5*Gamma0**5*logR**5*beamcache(1,1)*beamcache(2,1) +  &
             5*Ci**4*Gamma0**4*logR**4* &
              (2*beta0*beamcache(1,1)*beamcache(2,1) +  &
                beamcache(1,3)*beamcache(2,1) + beamcache(1,1)*beamcache(2,3)) &
              + 5*Ci**3*Gamma0**3*logR**3* &
              (7*beta0**2*beamcache(1,1)*beamcache(2,1) +  &
                6*beta0*(beamcache(1,3)*beamcache(2,1) +  &
                   beamcache(1,1)*beamcache(2,3)) +  &
                4*(beamcache(1,6)*beamcache(2,1) +  &
                   beamcache(1,3)*beamcache(2,3) + beamcache(1,1)*beamcache(2,6) &
                   )) + 5*Ci**2*Gamma0**2*logR**2* &
              (10*beta0**3*beamcache(1,1)*beamcache(2,1) +  &
                11*beta0**2*(beamcache(1,3)*beamcache(2,1) +  &
                   beamcache(1,1)*beamcache(2,3)) +  &
                12*beta0*(beamcache(1,6)*beamcache(2,1) +  &
                   beamcache(1,3)*beamcache(2,3) + beamcache(1,1)*beamcache(2,6) &
                   ) + 12*(beamcache(1,10)*beamcache(2,1) +  &
                   beamcache(1,6)*beamcache(2,3) +  &
                   beamcache(1,3)*beamcache(2,6) +  &
                   beamcache(1,1)*beamcache(2,10))) +  &
             2*Ci*Gamma0*logR*(12*beta0**4*beamcache(1,1)*beamcache(2,1) +  &
                15*beta0**3*(beamcache(1,3)*beamcache(2,1) +  &
                   beamcache(1,1)*beamcache(2,3)) +  &
                20*beta0**2*(beamcache(1,6)*beamcache(2,1) +  &
                   beamcache(1,3)*beamcache(2,3) + beamcache(1,1)*beamcache(2,6) &
                   ) + 30*beta0*(beamcache(1,10)*beamcache(2,1) +  &
                   beamcache(1,6)*beamcache(2,3) +  &
                   beamcache(1,3)*beamcache(2,6) +  &
                   beamcache(1,1)*beamcache(2,10)) +  &
                60*(beamcache(1,11)*beamcache(2,1) +  &
                   beamcache(1,10)*beamcache(2,3) +  &
                   beamcache(1,6)*beamcache(2,6) +  &
                   beamcache(1,3)*beamcache(2,10) +  &
                   beamcache(1,1)*beamcache(2,11))) +  &
             120*(beamcache(1,13)*beamcache(2,1) +  &
                beamcache(1,11)*beamcache(2,3) +  &
                beamcache(1,10)*beamcache(2,6) +  &
                beamcache(1,6)*beamcache(2,10) +  &
                beamcache(1,3)*beamcache(2,11) + beamcache(1,1)*beamcache(2,13)) &
             )*Mi(5))/120._dp


        ! g^6
        resnew = resnew + ason4pi**3*(-(gamma2rX*logR*beamcache(1,1)*beamcache(2,1)) +  &
           beamcache(1,7)*beamcache(2,1) + beamcache(1,4)*beamcache(2,2) -  &
           gamma1rX*logR*(beamcache(1,2)*beamcache(2,1) +  &
              beamcache(1,1)*beamcache(2,2)) + beamcache(1,2)*beamcache(2,4) +  &
           beamcache(1,1)*beamcache(2,7))*Mi(0) +  &
        (ason4pi**4*(-6*beta0**2*gamma1rX*logR*beamcache(1,1)*beamcache(2,1) +  &
             3*beta0*Ci*Gamma2*logR*beamcache(1,1)*beamcache(2,1) +  &
             Ci**2*Gamma1**2*logR**2*beamcache(1,1)*beamcache(2,1) -  &
             4*beta0*gamma1rX*logR*beamcache(1,3)*beamcache(2,1) +  &
             2*Ci*Gamma2*logR*beamcache(1,3)*beamcache(2,1) -  &
             2*gamma1rX*logR*beamcache(1,6)*beamcache(2,1) +  &
             2*beamcache(1,14)*beamcache(2,1) +  &
             2*beamcache(1,9)*beamcache(2,2) -  &
             4*beta0*gamma1rX*logR*beamcache(1,1)*beamcache(2,3) +  &
             2*Ci*Gamma2*logR*beamcache(1,1)*beamcache(2,3) -  &
             2*gamma1rX*logR*beamcache(1,3)*beamcache(2,3) +  &
             2*beamcache(1,8)*beamcache(2,3) +  &
             2*beamcache(1,6)*beamcache(2,4) +  &
             Ci**2*Gamma0**2*logR**2* &
              (-(gamma1rX*logR*beamcache(1,1)*beamcache(2,1)) +  &
                beamcache(1,4)*beamcache(2,1) + beamcache(1,2)*beamcache(2,2) +  &
                beamcache(1,1)*beamcache(2,4)) +  &
             2*beamcache(1,5)*beamcache(2,5) +  &
             2*Ci*Gamma1*logR*(beta1*beamcache(1,1)*beamcache(2,1) +  &
                beta0*beamcache(1,2)*beamcache(2,1) +  &
                beamcache(1,5)*beamcache(2,1) +  &
                beta0*beamcache(1,1)*beamcache(2,2) +  &
                beamcache(1,3)*beamcache(2,2) +  &
                Ci*Gamma0*logR*(beamcache(1,2)*beamcache(2,1) +  &
                   beamcache(1,1)*beamcache(2,2)) +  &
                beamcache(1,2)*beamcache(2,3) + beamcache(1,1)*beamcache(2,5)) &
              - 2*gamma1rX*logR*beamcache(1,1)*beamcache(2,6) +  &
             2*beamcache(1,4)*beamcache(2,6) +  &
             2*beamcache(1,3)*beamcache(2,8) +  &
             Ci*Gamma0*logR*(beta2*beamcache(1,1)*beamcache(2,1) +  &
                2*Ci*Gamma2*logR*beamcache(1,1)*beamcache(2,1) +  &
                beta1*beamcache(1,2)*beamcache(2,1) +  &
                beta0*beamcache(1,4)*beamcache(2,1) +  &
                2*beamcache(1,8)*beamcache(2,1) +  &
                beta1*beamcache(1,1)*beamcache(2,2) +  &
                beta0*beamcache(1,2)*beamcache(2,2) +  &
                2*beamcache(1,5)*beamcache(2,2) +  &
                2*beamcache(1,4)*beamcache(2,3) -  &
                gamma1rX*logR*(5*beta0*beamcache(1,1)*beamcache(2,1) +  &
                   2*beamcache(1,3)*beamcache(2,1) +  &
                   2*beamcache(1,1)*beamcache(2,3)) +  &
                beta0*beamcache(1,1)*beamcache(2,4) +  &
                2*beamcache(1,3)*beamcache(2,4) +  &
                2*beamcache(1,2)*beamcache(2,5) +  &
                2*beamcache(1,1)*beamcache(2,8)) +  &
             2*beamcache(1,2)*beamcache(2,9) + 2*beamcache(1,1)*beamcache(2,14)) &
            *Mi(2))/2._dp + (ason4pi**5* &
           (Ci**4*Gamma0**4*logR**4* &
              (beamcache(1,2)*beamcache(2,1) + beamcache(1,1)*beamcache(2,2)) +  &
             2*Ci**3*Gamma0**3*logR**3* &
              (3*beta1*beamcache(1,1)*beamcache(2,1) +  &
                2*Ci*Gamma1*logR*beamcache(1,1)*beamcache(2,1) +  &
                3*beta0*beamcache(1,2)*beamcache(2,1) +  &
                2*beamcache(1,5)*beamcache(2,1) +  &
                3*beta0*beamcache(1,1)*beamcache(2,2) +  &
                2*beamcache(1,3)*beamcache(2,2) +  &
                2*beamcache(1,2)*beamcache(2,3) +  &
                2*beamcache(1,1)*beamcache(2,5)) +  &
             Ci**2*Gamma0**2*logR**2* &
              (26*beta0*beta1*beamcache(1,1)*beamcache(2,1) +  &
                11*beta0**2*(beamcache(1,2)*beamcache(2,1) +  &
                   beamcache(1,1)*beamcache(2,2)) +  &
                12*Ci*Gamma1*logR* &
                 (2*beta0*beamcache(1,1)*beamcache(2,1) +  &
                   beamcache(1,3)*beamcache(2,1) + beamcache(1,1)*beamcache(2,3) &
                   ) + 12*beta0*(beamcache(1,5)*beamcache(2,1) +  &
                   beamcache(1,3)*beamcache(2,2) +  &
                   beamcache(1,2)*beamcache(2,3) + beamcache(1,1)*beamcache(2,5) &
                   ) + 12*(beta1*beamcache(1,3)*beamcache(2,1) +  &
                   beamcache(1,9)*beamcache(2,1) +  &
                   beamcache(1,6)*beamcache(2,2) +  &
                   beta1*beamcache(1,1)*beamcache(2,3) +  &
                   beamcache(1,5)*beamcache(2,3) +  &
                   beamcache(1,3)*beamcache(2,5) +  &
                   beamcache(1,2)*beamcache(2,6) + beamcache(1,1)*beamcache(2,9) &
                   )) + 2*Ci*Gamma0*logR* &
              (13*beta0**2*beta1*beamcache(1,1)*beamcache(2,1) +  &
                3*beta0**3*(beamcache(1,2)*beamcache(2,1) +  &
                   beamcache(1,1)*beamcache(2,2)) +  &
                10*beta0*beta1*(beamcache(1,3)*beamcache(2,1) +  &
                   beamcache(1,1)*beamcache(2,3)) +  &
                4*beta0**2*(beamcache(1,5)*beamcache(2,1) +  &
                   beamcache(1,3)*beamcache(2,2) +  &
                   beamcache(1,2)*beamcache(2,3) + beamcache(1,1)*beamcache(2,5) &
                   ) + 6*beta1*(beamcache(1,6)*beamcache(2,1) +  &
                   beamcache(1,3)*beamcache(2,3) + beamcache(1,1)*beamcache(2,6) &
                   ) + 2*Ci*Gamma1*logR* &
                 (11*beta0**2*beamcache(1,1)*beamcache(2,1) +  &
                   9*beta0*(beamcache(1,3)*beamcache(2,1) +  &
                      beamcache(1,1)*beamcache(2,3)) +  &
                   6*(beamcache(1,6)*beamcache(2,1) +  &
                      beamcache(1,3)*beamcache(2,3) +  &
                      beamcache(1,1)*beamcache(2,6))) +  &
                6*beta0*(beamcache(1,9)*beamcache(2,1) +  &
                   beamcache(1,6)*beamcache(2,2) +  &
                   beamcache(1,5)*beamcache(2,3) +  &
                   beamcache(1,3)*beamcache(2,5) +  &
                   beamcache(1,2)*beamcache(2,6) + beamcache(1,1)*beamcache(2,9) &
                   ) + 12*(beamcache(1,12)*beamcache(2,1) +  &
                   beamcache(1,10)*beamcache(2,2) +  &
                   beamcache(1,9)*beamcache(2,3) +  &
                   beamcache(1,6)*beamcache(2,5) +  &
                   beamcache(1,5)*beamcache(2,6) +  &
                   beamcache(1,3)*beamcache(2,9) +  &
                   beamcache(1,2)*beamcache(2,10) +  &
                   beamcache(1,1)*beamcache(2,12))) +  &
             24*(beamcache(1,15)*beamcache(2,1) +  &
                beamcache(1,11)*beamcache(2,2) +  &
                beamcache(1,12)*beamcache(2,3) +  &
                beamcache(1,10)*beamcache(2,5) +  &
                beamcache(1,9)*beamcache(2,6) + beamcache(1,6)*beamcache(2,9) +  &
                beamcache(1,5)*beamcache(2,10) +  &
                Ci*Gamma1*logR*(beta0**3*beamcache(1,1)*beamcache(2,1) +  &
                   beamcache(1,10)*beamcache(2,1) +  &
                   beamcache(1,6)*beamcache(2,3) +  &
                   beta0**2*(beamcache(1,3)*beamcache(2,1) +  &
                      beamcache(1,1)*beamcache(2,3)) +  &
                   beamcache(1,3)*beamcache(2,6) +  &
                   beta0*(beamcache(1,6)*beamcache(2,1) +  &
                      beamcache(1,3)*beamcache(2,3) +  &
                      beamcache(1,1)*beamcache(2,6)) +  &
                   beamcache(1,1)*beamcache(2,10)) +  &
                beamcache(1,2)*beamcache(2,11) +  &
                beamcache(1,3)*beamcache(2,12) + beamcache(1,1)*beamcache(2,15)) &
             )*Mi(4))/24._dp + (ason4pi**6* &
           (Ci**6*Gamma0**6*logR**6*beamcache(1,1)*beamcache(2,1) +  &
             3*Ci**5*Gamma0**5*logR**5* &
              (5*beta0*beamcache(1,1)*beamcache(2,1) +  &
                2*beamcache(1,3)*beamcache(2,1) +  &
                2*beamcache(1,1)*beamcache(2,3)) +  &
             5*Ci**4*Gamma0**4*logR**4* &
              (17*beta0**2*beamcache(1,1)*beamcache(2,1) +  &
                12*beta0*(beamcache(1,3)*beamcache(2,1) +  &
                   beamcache(1,1)*beamcache(2,3)) +  &
                6*(beamcache(1,6)*beamcache(2,1) +  &
                   beamcache(1,3)*beamcache(2,3) + beamcache(1,1)*beamcache(2,6) &
                   )) + 15*Ci**3*Gamma0**3*logR**3* &
              (15*beta0**3*beamcache(1,1)*beamcache(2,1) +  &
                14*beta0**2*(beamcache(1,3)*beamcache(2,1) +  &
                   beamcache(1,1)*beamcache(2,3)) +  &
                12*beta0*(beamcache(1,6)*beamcache(2,1) +  &
                   beamcache(1,3)*beamcache(2,3) + beamcache(1,1)*beamcache(2,6) &
                   ) + 8*(beamcache(1,10)*beamcache(2,1) +  &
                   beamcache(1,6)*beamcache(2,3) +  &
                   beamcache(1,3)*beamcache(2,6) +  &
                   beamcache(1,1)*beamcache(2,10))) +  &
             2*Ci**2*Gamma0**2*logR**2* &
              (137*beta0**4*beamcache(1,1)*beamcache(2,1) +  &
                150*beta0**3*(beamcache(1,3)*beamcache(2,1) +  &
                   beamcache(1,1)*beamcache(2,3)) +  &
                165*beta0**2*(beamcache(1,6)*beamcache(2,1) +  &
                   beamcache(1,3)*beamcache(2,3) + beamcache(1,1)*beamcache(2,6) &
                   ) + 180*beta0* &
                 (beamcache(1,10)*beamcache(2,1) +  &
                   beamcache(1,6)*beamcache(2,3) +  &
                   beamcache(1,3)*beamcache(2,6) +  &
                   beamcache(1,1)*beamcache(2,10)) +  &
                180*(beamcache(1,11)*beamcache(2,1) +  &
                   beamcache(1,10)*beamcache(2,3) +  &
                   beamcache(1,6)*beamcache(2,6) +  &
                   beamcache(1,3)*beamcache(2,10) +  &
                   beamcache(1,1)*beamcache(2,11))) +  &
             12*Ci*Gamma0*logR*(10*beta0**5*beamcache(1,1)*beamcache(2,1) +  &
                12*beta0**4*(beamcache(1,3)*beamcache(2,1) +  &
                   beamcache(1,1)*beamcache(2,3)) +  &
                15*beta0**3*(beamcache(1,6)*beamcache(2,1) +  &
                   beamcache(1,3)*beamcache(2,3) + beamcache(1,1)*beamcache(2,6) &
                   ) + 20*beta0**2* &
                 (beamcache(1,10)*beamcache(2,1) +  &
                   beamcache(1,6)*beamcache(2,3) +  &
                   beamcache(1,3)*beamcache(2,6) +  &
                   beamcache(1,1)*beamcache(2,10)) +  &
                30*beta0*(beamcache(1,11)*beamcache(2,1) +  &
                   beamcache(1,10)*beamcache(2,3) +  &
                   beamcache(1,6)*beamcache(2,6) +  &
                   beamcache(1,3)*beamcache(2,10) +  &
                   beamcache(1,1)*beamcache(2,11)) +  &
                60*(beamcache(1,13)*beamcache(2,1) +  &
                   beamcache(1,11)*beamcache(2,3) +  &
                   beamcache(1,10)*beamcache(2,6) +  &
                   beamcache(1,6)*beamcache(2,10) +  &
                   beamcache(1,3)*beamcache(2,11) +  &
                   beamcache(1,1)*beamcache(2,13))) +  &
             720*(beamcache(1,16)*beamcache(2,1) +  &
                beamcache(1,13)*beamcache(2,3) +  &
                beamcache(1,11)*beamcache(2,6) +  &
                beamcache(1,10)*beamcache(2,10) +  &
                beamcache(1,6)*beamcache(2,11) +  &
                beamcache(1,3)*beamcache(2,13) + beamcache(1,1)*beamcache(2,16)) &
             )*Mi(6))/720._dp

        endif

        res = hardEvolutionThres(q2,muH**2,alphasMuH,alphasMu,order,initQuark, &
            nf, nfmax, mub, alphas4p, alphas4m, muc, alphas3p, alphas3m)*resnew

        endif

        if (present(resexp)) then

        LM = log(q2/qt**2)
        LMU = log(muH**2/qt**2)

        resexp = 0._dp
        resexp1 = 0._dp
        resexp2 = 0._dp
        resexp3 = 0._dp

        ! beamfunctions are evaluated at muF

        beamcache(:,:) = 0._dp

        beamcache(1,1) = getbeam(ih1,f1,0,0,x1,facscale,1)
        beamcache(2,1) = getbeam(ih2,f2,0,0,x2,facscale,2)
        beamcache(1,3) = getbeam(ih1,f1,1,1,x1,facscale,1)
        beamcache(2,3) = getbeam(ih2,f2,1,1,x2,facscale,2)

        if (order >= 5) then
            beamcache(1,2) = getbeam(ih1,f1,1,0,x1,facscale,1)
            beamcache(1,5) = getbeam(ih1,f1,2,1,x1,facscale,1)
            beamcache(1,6) = getbeam(ih1,f1,2,2,x1,facscale,1)
            beamcache(2,2) = getbeam(ih2,f2,1,0,x2,facscale,2)
            beamcache(2,5) = getbeam(ih2,f2,2,1,x2,facscale,2)
            beamcache(2,6) = getbeam(ih2,f2,2,2,x2,facscale,2)
        endif

        if (order >= 7) then
            beamcache(1,4) = getbeam(ih1,f1,2,0,x1,facscale,1)
            beamcache(1,8) = getbeam(ih1,f1,3,1,x1,facscale,1)
            beamcache(1,9) = getbeam(ih1,f1,3,2,x1,facscale,1)
            beamcache(1,10) = getbeam(ih1,f1,3,3,x1,facscale,1)
            beamcache(2,4) = getbeam(ih2,f2,2,0,x2,facscale,2)
            beamcache(2,8) = getbeam(ih2,f2,3,1,x2,facscale,2)
            beamcache(2,9) = getbeam(ih2,f2,3,2,x2,facscale,2)
            beamcache(2,10) = getbeam(ih2,f2,3,3,x2,facscale,2)


            !! TO SORT

            beamcache(1,7) = getbeam(ih1,f1,3,0,x1,facscale,1)
            beamcache(2,7) = getbeam(ih2,f2,3,0,x2,facscale,2)

        endif

        logR = 0._dp

        resexp1 = 2*gammai0*beamcache(1,1)*beamcache(2,1) + &
            Ci*Gamma0*LM*beamcache(1,1)*beamcache(2,1) - &
            beamcache(1,3)*beamcache(2,1) - beamcache(1,1)*beamcache(2,3)

        resexp = resexp1 * alphasMuMatch/4._dp/pi / qt**2

        if (order >= 5) then
            resexp2 = -(gamma1rX*beamcache(1,1)*beamcache(2,1)) + &
                2*gammai1*beamcache(1,1)*beamcache(2,1) + &
                2*gammai0*hard1*beamcache(1,1)*beamcache(2,1) + &
                Ci*Gamma1*LM*beamcache(1,1)*beamcache(2,1) + &
                Ci*Gamma0*hard1*LM*beamcache(1,1)*beamcache(2,1) + &
                2*beta0*gammai0*LMU*beamcache(1,1)*beamcache(2,1) - &
                4*gammai0**2*LMU*beamcache(1,1)*beamcache(2,1) + &
                beta0*Ci*Gamma0*LM*LMU*beamcache(1,1)*beamcache(2,1) - &
                4*Ci*Gamma0*gammai0*LM*LMU*beamcache(1,1)*beamcache(2,1) - &
                Ci**2*Gamma0**2*LM**2*LMU*beamcache(1,1)*beamcache(2,1) + &
                Ci*Gamma0*gammai0*LMU**2*beamcache(1,1)*beamcache(2,1) + &
                (Ci**2*Gamma0**2*LM*LMU**2*beamcache(1,1)*beamcache(2,1))/2._dp + &
                2*Ci**2*Gamma0**2*Zeta3*beamcache(1,1)*beamcache(2,1) + &
                2*gammai0*beamcache(1,2)*beamcache(2,1) + &
                Ci*Gamma0*LM*beamcache(1,2)*beamcache(2,1) - &
                hard1*beamcache(1,3)*beamcache(2,1) + &
                4*gammai0*LMU*beamcache(1,3)*beamcache(2,1) + &
                2*Ci*Gamma0*LM*LMU*beamcache(1,3)*beamcache(2,1) - &
                (Ci*Gamma0*LMU**2*beamcache(1,3)*beamcache(2,1))/2._dp - &
                beamcache(1,5)*beamcache(2,1) - 2*LMU*beamcache(1,6)*beamcache(2,1) + &
                2*gammai0*beamcache(1,1)*beamcache(2,2) + &
                Ci*Gamma0*LM*beamcache(1,1)*beamcache(2,2) - &
                beamcache(1,3)*beamcache(2,2) - hard1*beamcache(1,1)*beamcache(2,3) + &
                4*gammai0*LMU*beamcache(1,1)*beamcache(2,3) + &
                2*Ci*Gamma0*LM*LMU*beamcache(1,1)*beamcache(2,3) - &
                (Ci*Gamma0*LMU**2*beamcache(1,1)*beamcache(2,3))/2._dp - &
                beamcache(1,2)*beamcache(2,3) - 2*LMU*beamcache(1,3)*beamcache(2,3) - &
                beamcache(1,1)*beamcache(2,5) - 2*LMU*beamcache(1,1)*beamcache(2,6)
            if (resexp_linPC) then
                ! only return NNLO coefficient
                resexp = resexp2 * (alphasMuMatch/4._dp/pi)**2 / qt**2
            else
                resexp = resexp + resexp2 * (alphasMuMatch/4._dp/pi)**2 / qt**2
            endif
        endif

        if (order >= 7) then
            resexp3 = -(gamma2rX*beamcache(1,1)*beamcache(2,1)) + &
                2*gammai2*beamcache(1,1)*beamcache(2,1) - &
                gamma1rX*hard1*beamcache(1,1)*beamcache(2,1) + &
                2*gammai1*hard1*beamcache(1,1)*beamcache(2,1) + &
                2*gammai0*hard2*beamcache(1,1)*beamcache(2,1) - &
                2*beta0*gamma1rX*LM*beamcache(1,1)*beamcache(2,1) + &
                Ci*Gamma2*LM*beamcache(1,1)*beamcache(2,1) + &
                2*gamma1rX*gammai0*LM*beamcache(1,1)*beamcache(2,1) + &
                Ci*Gamma1*hard1*LM*beamcache(1,1)*beamcache(2,1) + &
                Ci*Gamma0*hard2*LM*beamcache(1,1)*beamcache(2,1) + &
                Ci*Gamma0*gamma1rX*LM**2*beamcache(1,1)*beamcache(2,1) - &
                2*beta0*gamma1rX*LMU*beamcache(1,1)*beamcache(2,1) + &
                2*beta1*gammai0*LMU*beamcache(1,1)*beamcache(2,1) + &
                2*gamma1rX*gammai0*LMU*beamcache(1,1)*beamcache(2,1) + &
                4*beta0*gammai1*LMU*beamcache(1,1)*beamcache(2,1) - &
                8*gammai0*gammai1*LMU*beamcache(1,1)*beamcache(2,1) + &
                2*beta0*gammai0*hard1*LMU*beamcache(1,1)*beamcache(2,1) - &
                4*gammai0**2*hard1*LMU*beamcache(1,1)*beamcache(2,1) + &
                beta1*Ci*Gamma0*LM*LMU*beamcache(1,1)*beamcache(2,1) + &
                2*beta0*Ci*Gamma1*LM*LMU*beamcache(1,1)*beamcache(2,1) + &
                Ci*Gamma0*gamma1rX*LM*LMU*beamcache(1,1)*beamcache(2,1) - &
                4*Ci*Gamma1*gammai0*LM*LMU*beamcache(1,1)*beamcache(2,1) - &
                4*Ci*Gamma0*gammai1*LM*LMU*beamcache(1,1)*beamcache(2,1) + &
                beta0*Ci*Gamma0*hard1*LM*LMU*beamcache(1,1)*beamcache(2,1) - &
                4*Ci*Gamma0*gammai0*hard1*LM*LMU*beamcache(1,1)*beamcache(2,1) - &
                2*Ci**2*Gamma0*Gamma1*LM**2*LMU*beamcache(1,1)*beamcache(2,1) - &
                Ci**2*Gamma0**2*hard1*LM**2*LMU*beamcache(1,1)*beamcache(2,1) - &
                (Ci*Gamma0*gamma1rX*LMU**2*beamcache(1,1)*beamcache(2,1))/2._dp + &
                2*beta0**2*gammai0*LMU**2*beamcache(1,1)*beamcache(2,1) + &
                Ci*Gamma1*gammai0*LMU**2*beamcache(1,1)*beamcache(2,1) - &
                6*beta0*gammai0**2*LMU**2*beamcache(1,1)*beamcache(2,1) + &
                4*gammai0**3*LMU**2*beamcache(1,1)*beamcache(2,1) + &
                Ci*Gamma0*gammai1*LMU**2*beamcache(1,1)*beamcache(2,1) + &
                Ci*Gamma0*gammai0*hard1*LMU**2*beamcache(1,1)*beamcache(2,1) + &
                beta0**2*Ci*Gamma0*LM*LMU**2*beamcache(1,1)*beamcache(2,1) + &
                Ci**2*Gamma0*Gamma1*LM*LMU**2*beamcache(1,1)*beamcache(2,1) - &
                6*beta0*Ci*Gamma0*gammai0*LM*LMU**2*beamcache(1,1)*beamcache(2,1) + &
                6*Ci*Gamma0*gammai0**2*LM*LMU**2*beamcache(1,1)*beamcache(2,1) + &
                (Ci**2*Gamma0**2*hard1*LM*LMU**2*beamcache(1,1)*beamcache(2,1))/2._dp - &
                (3*beta0*Ci**2*Gamma0**2*LM**2*LMU**2*beamcache(1,1)*beamcache(2,1))/&
                 2._dp + 3*Ci**2*Gamma0**2*gammai0*LM**2*LMU**2*beamcache(1,1)*&
                 beamcache(2,1) + (Ci**3*Gamma0**3*LM**3*LMU**2*beamcache(1,1)*&
                   beamcache(2,1))/2._dp + (4*beta0*Ci*Gamma0*gammai0*LMU**3*&
                   beamcache(1,1)*beamcache(2,1))/3._dp - &
                2*Ci*Gamma0*gammai0**2*LMU**3*beamcache(1,1)*beamcache(2,1) + &
                (2*beta0*Ci**2*Gamma0**2*LM*LMU**3*beamcache(1,1)*beamcache(2,1))/3._dp - &
                2*Ci**2*Gamma0**2*gammai0*LM*LMU**3*beamcache(1,1)*beamcache(2,1) - &
                (Ci**3*Gamma0**3*LM**2*LMU**3*beamcache(1,1)*beamcache(2,1))/2._dp + &
                (Ci**2*Gamma0**2*gammai0*LMU**4*beamcache(1,1)*beamcache(2,1))/4._dp + &
                (Ci**3*Gamma0**3*LM*LMU**4*beamcache(1,1)*beamcache(2,1))/8._dp - &
                4*beta0**2*Ci*Gamma0*Zeta3*beamcache(1,1)*beamcache(2,1) + &
                4*Ci**2*Gamma0*Gamma1*Zeta3*beamcache(1,1)*beamcache(2,1) + &
                (56*beta0*Ci*Gamma0*gammai0*Zeta3*beamcache(1,1)*beamcache(2,1))/3._dp - &
                16*Ci*Gamma0*gammai0**2*Zeta3*beamcache(1,1)*beamcache(2,1) + &
                2*Ci**2*Gamma0**2*hard1*Zeta3*beamcache(1,1)*beamcache(2,1) + &
                (28*beta0*Ci**2*Gamma0**2*LM*Zeta3*beamcache(1,1)*beamcache(2,1))/3._dp - &
                16*Ci**2*Gamma0**2*gammai0*LM*Zeta3*beamcache(1,1)*beamcache(2,1) - &
                4*Ci**3*Gamma0**3*LM**2*Zeta3*beamcache(1,1)*beamcache(2,1) + &
                4*beta0*Ci**2*Gamma0**2*LMU*Zeta3*beamcache(1,1)*beamcache(2,1) - &
                4*Ci**2*Gamma0**2*gammai0*LMU*Zeta3*beamcache(1,1)*beamcache(2,1) - &
                2*Ci**3*Gamma0**3*LM*LMU*Zeta3*beamcache(1,1)*beamcache(2,1) + &
                Ci**3*Gamma0**3*LMU**2*Zeta3*beamcache(1,1)*beamcache(2,1) - &
                6*Ci**3*Gamma0**3*Zeta5*beamcache(1,1)*beamcache(2,1) - &
                gamma1rX*beamcache(1,2)*beamcache(2,1) + &
                2*gammai1*beamcache(1,2)*beamcache(2,1) + &
                2*gammai0*hard1*beamcache(1,2)*beamcache(2,1) + &
                Ci*Gamma1*LM*beamcache(1,2)*beamcache(2,1) + &
                Ci*Gamma0*hard1*LM*beamcache(1,2)*beamcache(2,1) + &
                2*beta0*gammai0*LMU*beamcache(1,2)*beamcache(2,1) - &
                4*gammai0**2*LMU*beamcache(1,2)*beamcache(2,1) + &
                beta0*Ci*Gamma0*LM*LMU*beamcache(1,2)*beamcache(2,1) - &
                4*Ci*Gamma0*gammai0*LM*LMU*beamcache(1,2)*beamcache(2,1) - &
                Ci**2*Gamma0**2*LM**2*LMU*beamcache(1,2)*beamcache(2,1) + &
                Ci*Gamma0*gammai0*LMU**2*beamcache(1,2)*beamcache(2,1) + &
                (Ci**2*Gamma0**2*LM*LMU**2*beamcache(1,2)*beamcache(2,1))/2._dp + &
                2*Ci**2*Gamma0**2*Zeta3*beamcache(1,2)*beamcache(2,1) - &
                hard2*beamcache(1,3)*beamcache(2,1) - &
                gamma1rX*LM*beamcache(1,3)*beamcache(2,1) - &
                gamma1rX*LMU*beamcache(1,3)*beamcache(2,1) + &
                4*gammai1*LMU*beamcache(1,3)*beamcache(2,1) + &
                4*gammai0*hard1*LMU*beamcache(1,3)*beamcache(2,1) + &
                2*Ci*Gamma1*LM*LMU*beamcache(1,3)*beamcache(2,1) + &
                2*Ci*Gamma0*hard1*LM*LMU*beamcache(1,3)*beamcache(2,1) - &
                (Ci*Gamma1*LMU**2*beamcache(1,3)*beamcache(2,1))/2._dp + &
                3*beta0*gammai0*LMU**2*beamcache(1,3)*beamcache(2,1) - &
                6*gammai0**2*LMU**2*beamcache(1,3)*beamcache(2,1) - &
                (Ci*Gamma0*hard1*LMU**2*beamcache(1,3)*beamcache(2,1))/2._dp + &
                (3*beta0*Ci*Gamma0*LM*LMU**2*beamcache(1,3)*beamcache(2,1))/2._dp - &
                6*Ci*Gamma0*gammai0*LM*LMU**2*beamcache(1,3)*beamcache(2,1) - &
                (3*Ci**2*Gamma0**2*LM**2*LMU**2*beamcache(1,3)*beamcache(2,1))/2._dp - &
                (beta0*Ci*Gamma0*LMU**3*beamcache(1,3)*beamcache(2,1))/6._dp + &
                2*Ci*Gamma0*gammai0*LMU**3*beamcache(1,3)*beamcache(2,1) + &
                Ci**2*Gamma0**2*LM*LMU**3*beamcache(1,3)*beamcache(2,1) - &
                (Ci**2*Gamma0**2*LMU**4*beamcache(1,3)*beamcache(2,1))/8._dp - &
                (16*beta0*Ci*Gamma0*Zeta3*beamcache(1,3)*beamcache(2,1))/3._dp + &
                16*Ci*Gamma0*gammai0*Zeta3*beamcache(1,3)*beamcache(2,1) + &
                8*Ci**2*Gamma0**2*LM*Zeta3*beamcache(1,3)*beamcache(2,1) + &
                2*Ci**2*Gamma0**2*LMU*Zeta3*beamcache(1,3)*beamcache(2,1) + &
                2*gammai0*beamcache(1,4)*beamcache(2,1) + &
                Ci*Gamma0*LM*beamcache(1,4)*beamcache(2,1) - &
                hard1*beamcache(1,5)*beamcache(2,1) + &
                4*gammai0*LMU*beamcache(1,5)*beamcache(2,1) + &
                2*Ci*Gamma0*LM*LMU*beamcache(1,5)*beamcache(2,1) - &
                (Ci*Gamma0*LMU**2*beamcache(1,5)*beamcache(2,1))/2._dp - &
                2*hard1*LMU*beamcache(1,6)*beamcache(2,1) + &
                6*gammai0*LMU**2*beamcache(1,6)*beamcache(2,1) + &
                3*Ci*Gamma0*LM*LMU**2*beamcache(1,6)*beamcache(2,1) - &
                Ci*Gamma0*LMU**3*beamcache(1,6)*beamcache(2,1) - &
                8*Ci*Gamma0*Zeta3*beamcache(1,6)*beamcache(2,1) - &
                beamcache(1,8)*beamcache(2,1) - 2*LMU*beamcache(1,9)*beamcache(2,1) - &
                3*LMU**2*beamcache(1,10)*beamcache(2,1) - &
                gamma1rX*beamcache(1,1)*beamcache(2,2) + &
                2*gammai1*beamcache(1,1)*beamcache(2,2) + &
                2*gammai0*hard1*beamcache(1,1)*beamcache(2,2) + &
                Ci*Gamma1*LM*beamcache(1,1)*beamcache(2,2) + &
                Ci*Gamma0*hard1*LM*beamcache(1,1)*beamcache(2,2) + &
                2*beta0*gammai0*LMU*beamcache(1,1)*beamcache(2,2) - &
                4*gammai0**2*LMU*beamcache(1,1)*beamcache(2,2) + &
                beta0*Ci*Gamma0*LM*LMU*beamcache(1,1)*beamcache(2,2) - &
                4*Ci*Gamma0*gammai0*LM*LMU*beamcache(1,1)*beamcache(2,2) - &
                Ci**2*Gamma0**2*LM**2*LMU*beamcache(1,1)*beamcache(2,2) + &
                Ci*Gamma0*gammai0*LMU**2*beamcache(1,1)*beamcache(2,2) + &
                (Ci**2*Gamma0**2*LM*LMU**2*beamcache(1,1)*beamcache(2,2))/2._dp + &
                2*Ci**2*Gamma0**2*Zeta3*beamcache(1,1)*beamcache(2,2) + &
                2*gammai0*beamcache(1,2)*beamcache(2,2) + &
                Ci*Gamma0*LM*beamcache(1,2)*beamcache(2,2) - &
                hard1*beamcache(1,3)*beamcache(2,2) + &
                4*gammai0*LMU*beamcache(1,3)*beamcache(2,2) + &
                2*Ci*Gamma0*LM*LMU*beamcache(1,3)*beamcache(2,2) - &
                (Ci*Gamma0*LMU**2*beamcache(1,3)*beamcache(2,2))/2._dp - &
                beamcache(1,5)*beamcache(2,2) - 2*LMU*beamcache(1,6)*beamcache(2,2) - &
                hard2*beamcache(1,1)*beamcache(2,3) - &
                gamma1rX*LM*beamcache(1,1)*beamcache(2,3) - &
                gamma1rX*LMU*beamcache(1,1)*beamcache(2,3) + &
                4*gammai1*LMU*beamcache(1,1)*beamcache(2,3) + &
                4*gammai0*hard1*LMU*beamcache(1,1)*beamcache(2,3) + &
                2*Ci*Gamma1*LM*LMU*beamcache(1,1)*beamcache(2,3) + &
                2*Ci*Gamma0*hard1*LM*LMU*beamcache(1,1)*beamcache(2,3) - &
                (Ci*Gamma1*LMU**2*beamcache(1,1)*beamcache(2,3))/2._dp + &
                3*beta0*gammai0*LMU**2*beamcache(1,1)*beamcache(2,3) - &
                6*gammai0**2*LMU**2*beamcache(1,1)*beamcache(2,3) - &
                (Ci*Gamma0*hard1*LMU**2*beamcache(1,1)*beamcache(2,3))/2._dp + &
                (3*beta0*Ci*Gamma0*LM*LMU**2*beamcache(1,1)*beamcache(2,3))/2._dp - &
                6*Ci*Gamma0*gammai0*LM*LMU**2*beamcache(1,1)*beamcache(2,3) - &
                (3*Ci**2*Gamma0**2*LM**2*LMU**2*beamcache(1,1)*beamcache(2,3))/2._dp - &
                (beta0*Ci*Gamma0*LMU**3*beamcache(1,1)*beamcache(2,3))/6._dp + &
                2*Ci*Gamma0*gammai0*LMU**3*beamcache(1,1)*beamcache(2,3) + &
                Ci**2*Gamma0**2*LM*LMU**3*beamcache(1,1)*beamcache(2,3) - &
                (Ci**2*Gamma0**2*LMU**4*beamcache(1,1)*beamcache(2,3))/8._dp - &
                (16*beta0*Ci*Gamma0*Zeta3*beamcache(1,1)*beamcache(2,3))/3._dp + &
                16*Ci*Gamma0*gammai0*Zeta3*beamcache(1,1)*beamcache(2,3) + &
                8*Ci**2*Gamma0**2*LM*Zeta3*beamcache(1,1)*beamcache(2,3) + &
                2*Ci**2*Gamma0**2*LMU*Zeta3*beamcache(1,1)*beamcache(2,3) - &
                hard1*beamcache(1,2)*beamcache(2,3) + &
                4*gammai0*LMU*beamcache(1,2)*beamcache(2,3) + &
                2*Ci*Gamma0*LM*LMU*beamcache(1,2)*beamcache(2,3) - &
                (Ci*Gamma0*LMU**2*beamcache(1,2)*beamcache(2,3))/2._dp - &
                2*hard1*LMU*beamcache(1,3)*beamcache(2,3) + &
                6*gammai0*LMU**2*beamcache(1,3)*beamcache(2,3) + &
                3*Ci*Gamma0*LM*LMU**2*beamcache(1,3)*beamcache(2,3) - &
                Ci*Gamma0*LMU**3*beamcache(1,3)*beamcache(2,3) - &
                8*Ci*Gamma0*Zeta3*beamcache(1,3)*beamcache(2,3) - &
                beamcache(1,4)*beamcache(2,3) - 2*LMU*beamcache(1,5)*beamcache(2,3) - &
                3*LMU**2*beamcache(1,6)*beamcache(2,3) + &
                2*gammai0*beamcache(1,1)*beamcache(2,4) + &
                Ci*Gamma0*LM*beamcache(1,1)*beamcache(2,4) - &
                beamcache(1,3)*beamcache(2,4) - hard1*beamcache(1,1)*beamcache(2,5) + &
                4*gammai0*LMU*beamcache(1,1)*beamcache(2,5) + &
                2*Ci*Gamma0*LM*LMU*beamcache(1,1)*beamcache(2,5) - &
                (Ci*Gamma0*LMU**2*beamcache(1,1)*beamcache(2,5))/2._dp - &
                beamcache(1,2)*beamcache(2,5) - 2*LMU*beamcache(1,3)*beamcache(2,5) - &
                2*hard1*LMU*beamcache(1,1)*beamcache(2,6) + &
                6*gammai0*LMU**2*beamcache(1,1)*beamcache(2,6) + &
                3*Ci*Gamma0*LM*LMU**2*beamcache(1,1)*beamcache(2,6) - &
                Ci*Gamma0*LMU**3*beamcache(1,1)*beamcache(2,6) - &
                8*Ci*Gamma0*Zeta3*beamcache(1,1)*beamcache(2,6) - &
                2*LMU*beamcache(1,2)*beamcache(2,6) - &
                3*LMU**2*beamcache(1,3)*beamcache(2,6) - &
                beamcache(1,1)*beamcache(2,8) - 2*LMU*beamcache(1,1)*beamcache(2,9) - &
                3*LMU**2*beamcache(1,1)*beamcache(2,10)

            if (resexp_linPC) then
                ! only return N3LO coefficient
                resexp = resexp3 * (alphasMuMatch/4._dp/pi)**3 / qt**2
            else
                resexp = resexp3 * (alphasMuMatch/4._dp/pi)**3 / qt**2
            endif
        endif


        endif

    end subroutine

    function gen_res(r,qt,p,wt, y_in)
      use qtResummation_params
      use types
      implicit none
      include 'src/Inc/constants.f'
      include 'src/Inc/mxpart.f'
      include 'src/Inc/mxdim.f'
      include 'src/Inc/breit.f'
      include 'src/Inc/limits.f'
      include 'src/Inc/energy.f'
      include 'src/Inc/x1x2.f'

      logical :: gen_res
      real(dp), intent(in) :: r(mxdim)
      real(dp), intent(in) :: qt
      real(dp), intent(out) :: p(mxpart,4)
      real(dp), intent(out) :: wt
      real(dp), optional, intent(in) :: y_in

      real(dp), parameter :: wt0 = 1._dp/2._dp/pi
      real(dp) :: Qsqmin, Qsqmax, wtbw, wt12, wt34
      real(dp) :: tau, Q(4), Qsq
      real(dp) :: y, ymin, ymax, xa, xb

      real(dp) :: massvec

      p(:,:) = 0._dp

      Qsqmin = max(wsqmin,1._dp)
      Qsqmax = min(wsqmax, sqrts**2*0.999999_dp)

      if (enable_dsigma_dQ) then
          !!! BENCHMARK !!! This fixes Qsq = mZ^2 and calculates
          ! cross-section as dsigma/dQ
          Qsq = mass3**2
          wtbw = 2*mass3
      else
          call breitw(r(1), Qsqmin, Qsqmax, mass3, width3, Qsq, wtbw)
      endif


      ! one approximation
      tau = (Qsq + qt**2)/sqrts**2

      ! "better" approximation
      !tau = ((sqrts**2 + Qsq) - sqrt((Qsq-sqrts**2)**2 - 4._dp * qt**2 * sqrts**2)) &
                !/ (2._dp * sqrts * sqrt(Qsq + qt**2))
      !tau = tau**2

      ymin = log(sqrt(tau))
      ymax = -ymin
      if (present(y_in)) then 
          y = y_in
          wt12 = 1._dp
      elseif (enable_fixed_y) then
          y = fixed_y
          wt12 = 1._dp
      else
          y = ymin + r(2)*(ymax-ymin)
          wt12 = ymax-ymin
      endif
      xx(1) = min(sqrt(tau)*exp(y), 1._dp)
      xx(2) = min(sqrt(tau)*exp(-y), 1._dp)

      p(1,4) = 0.5_dp*sqrts*xx(1)
      p(1,3) = p(1,4)
      p(2,4) = 0.5_dp*sqrts*xx(2)
      p(2,3) = -p(2,4)

      wt12 = 2._dp*pi*wt12/sqrts**2

      Q = p(1,:) + p(2,:)

      !!! BENCHMARK !!!
      if (enable_fixed_y .and. enable_dsigma_dQ .and. (qtminRes == qtmaxRes)) then
          call phi3m0(r(1), r(2), Q, p(3,:), p(4,:), wt34, *99)
      else
          call phi3m0(r(3), r(4), Q, p(3,:), p(4,:), wt34, *99)
      endif

      p(1,:) = -p(1,:)
      p(2,:) = -p(2,:)


      wt = wt0*wt12*wtbw*wt34

      ! for debugging
      !wt = ymax-ymin

      gen_res = .true.

      return

 99   continue
      wt = 0
      gen_res = .false.

    end function

end module
