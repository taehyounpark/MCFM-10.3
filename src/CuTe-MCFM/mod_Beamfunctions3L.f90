!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module Beamfunctions3L
    use LHAPDF
    use types
    use constants
    implicit none

    public :: getbeam
    public :: getbeam2

    ! for debugging
    public :: ibar_select
    public :: set_params, set_nf

!   public :: check_Ibar_select

    logical, public, save :: usegrid, makegrid

    private

    ! these are used by the integrands
    integer, save :: flavor
!$omp threadprivate(flavor)
    real(dp), save :: xi, mu
!$omp threadprivate(xi,mu)
    integer, save :: flavor_i, powAs, powLperp, ih
!$omp threadprivate(flavor_i, powAs, powLperp, ih)
    integer, save :: nf
!$omp threadprivate(nf)

    integer, save :: usebeam
!$omp threadprivate(usebeam)

    real(dp), parameter :: integrationPrecision = 1e-4_dp


    integer, parameter :: n1=-1, n2=1, nw = 6

    complex(dp), save :: Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), &
                    Hc4(n1:n2,n1:n2,n1:n2,n1:n2),Hc5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2), &
                    Hc6(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2,n1:n2)
!$omp threadprivate(Hc1,Hc2,Hc3,Hc4,Hc5,Hc6)
    real(dp), save :: Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2), &
                    Hr4(n1:n2,n1:n2,n1:n2,n1:n2),Hr5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2), &
                    Hr6(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2,n1:n2)
!$omp threadprivate(Hr1,Hr2,Hr3,Hr4,Hr5,Hr6)
    real(dp), save :: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), &
                    Hi4(n1:n2,n1:n2,n1:n2,n1:n2),Hi5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2), &
                    Hi6(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2,n1:n2)
!$omp threadprivate(Hi1,Hi2,Hi3,Hi4,Hi5,Hi6)

    contains

    subroutine set_nf(nf_in)
        implicit none
        integer, intent(in) :: nf_in

        nf = nf_in
    end subroutine

    subroutine set_params(ih_in, flavor_i_in, powAs_in, powLperp_in, nf_in)
        implicit none
        integer, intent(in) :: ih_in
        integer, intent(in) :: flavor_i_in
        integer, intent(in) :: powAs_in, powLperp_in
        integer, intent(in) :: nf_in

        ih = ih_in
        flavor_i = flavor_i_in
        powAs = powAs_in
        powLperp = powLperp_in
        nf = nf_in
    end subroutine

    function getbeam2(ih_in, flavor_i_in, powAs_in, powLperp_in, xi_in, mu_in)
        use iso_fortran_env
        use Quadpack
        implicit none

        integer, intent(in) :: ih_in
        integer, intent(in) :: flavor_i_in
        integer, intent(in) :: powAs_in, powLperp_in
        real(dp), intent(in) :: xi_in, mu_in
        real(dp) :: getbeam2

        real(dp) :: absErr
        integer(int64) :: neval, ierr
        real(dp) :: beam_contrib
        integer :: flavor_j
        real(dp) :: pdfone

        if (flavor_i_in /= 0) then
            getbeam2 = 0._dp
        endif

        if (usegrid .and. (.not. makegrid)) then
            getbeam2 = fdist_one_beam2(ih_in, powAs_in, powLperp_in, xi_in, mu_in, flavor_i_in)
            return
        endif

        ih = ih_in
        flavor_i = flavor_i_in
        powAs = powAs_in
        powLperp = powLperp_in
        xi = xi_in
        mu = mu_in

        nf = getnumflavors(mu_in)

        getbeam2 = 0._dp

        ! perform getbeam_integrand integral
        beam_contrib = 0._dp
        call qags(getbeam2_integrand, xi, 1._dp, &
            integrationPrecision, integrationPrecision, &
            beam_contrib, absErr, neval, ierr)
        getbeam2 = beam_contrib

    end function

    function getbeam2_integrand(z)
        implicit none

        real(dp), intent(in) :: z
        real(dp) :: getbeam2_integrand

        integer :: flavor_j
        real(dp) :: pdfone, pdfz

        getbeam2_integrand = 0._dp

        do flavor_j=-nf,nf
            pdfz = fdist_one(ih, xi/z, mu, flavor_j)

            if (flavor_i == 0 .and. flavor_j == 0) then
                getbeam2_integrand = getbeam2_integrand + &
                    (-4*CA*(1-z)/z) / z * pdfz
            elseif (flavor_i == 0 .and. flavor_j /= 0) then
                getbeam2_integrand = getbeam2_integrand + &
                    (-4*CF*(1-z)/z) / z * pdfz
            endif
        enddo
    end function

    ! this gets the beam function component for a specific as and Lperp component
    function getbeam(ih_in, flavor_i_in, powAs_in, powLperp_in, xi_in, mu_in, ibeam_in)
        use iso_fortran_env
        use Quadpack
        implicit none

        integer, intent(in) :: ih_in
        integer, intent(in) :: flavor_i_in
        integer, intent(in) :: powAs_in, powLperp_in
        real(dp), intent(in) :: xi_in, mu_in
        integer, optional, intent(in) :: ibeam_in
        real(dp) :: getbeam

        real(dp) :: absErr
        integer(int64) :: neval, ierr
        real(dp) :: beam_contrib
        integer :: flavor_j
        real(dp) :: pdfone
        integer :: ibeam

        if (present(ibeam_in)) then
            ibeam = ibeam_in
        else
            ibeam = 1
        endif
        ! save globally
        usebeam = ibeam

        if (usegrid .and. (.not. makegrid)) then
            getbeam = fdist_one_beam(ih_in, powAs_in, powLperp_in, xi_in, mu_in, flavor_i_in, ibeam)
            return
        endif

        ih = ih_in
        flavor_i = flavor_i_in
        powAs = powAs_in
        powLperp = powLperp_in
        xi = xi_in
        mu = mu_in
        
        nf = getnumflavors(mu_in)

        getbeam = 0._dp


        ! perform getbeam_integrand integral
        beam_contrib = 0._dp
        call qags(getbeam_integrand, xi, 1._dp, &
            integrationPrecision, integrationPrecision, &
            beam_contrib, absErr, neval, ierr)
        getbeam = getbeam + beam_contrib

        do flavor_j=-nf,nf
             pdfone = fdist_one(ih, xi, mu, flavor_j, ibeam)

            ! delta(1-z) contribution
            getbeam = getbeam + &
                Ibar_select(flavor_j, 2, 1._dp) * pdfone

            ! boundary term for (1/(1-z))_+
            if (powAs > 0) then
                if (xi < 1._dp) then
                    getbeam = getbeam + &
                        Ibar_select(flavor_j, 3, 1._dp) * pdfone * log(1._dp - xi)
                endif
            endif

            ! boundary term for (log(1-z)/(1-z))_+
            if (powAs > 1) then
                if (xi < 1._dp) then
                    getbeam = getbeam + &
                        Ibar_select(flavor_j, 4, 1._dp) * pdfone * log(1._dp - xi)**2 / 2._dp
                endif
            endif

            ! boundary term for ((log(1-z))^2/(1-z))_+
            if (powAs > 2) then
                if (xi < 1._dp) then
                    getbeam = getbeam + &
                        Ibar_select(flavor_j, 5, 1._dp) * pdfone * log(1._dp - xi)**3 / 3._dp
                endif
            endif

            ! boundary term for ((log(1-z))^3/(1-z))_+
            if (powAs > 3) then
                if (xi < 1._dp) then
                    getbeam = getbeam + &
                        Ibar_select(flavor_j, 6, 1._dp) * pdfone * log(1._dp - xi)**4 / 4._dp
                endif
            endif

            ! boundary term for ((log(1-z))^4/(1-z))_+
            if (powAs > 4) then
                if (xi < 1._dp) then
                    getbeam = getbeam + &
                        Ibar_select(flavor_j, 7, 1._dp) * pdfone * log(1._dp - xi)**5 / 5._dp
                endif
            endif

            ! boundary term for ((log(1-z))^5/(1-z))_+
            if (powAs > 5) then
                if (xi < 1._dp) then
                    getbeam = getbeam + &
                        Ibar_select(flavor_j, 8, 1._dp) * pdfone * log(1._dp - xi)**6 / 6._dp
                endif
            endif


        enddo

    end function

    function getbeam_integrand(z)
        implicit none

        real(dp), intent(in) :: z
        real(dp) :: getbeam_integrand

        integer :: flavor_j
        real(dp) :: pdfone, pdfz

        getbeam_integrand = 0._dp

        ! we're only using exact HPL's for higher powers
        ! for lower powers the interpolated expressions from 2012.03256 are used
        if (powAs >= 4) then
            call inithpl(z)
        endif

        do flavor_j=-nf,nf
            pdfz = fdist_one(ih, xi/z, mu, flavor_j, usebeam)
            ! O(1) contribution
            getbeam_integrand = getbeam_integrand + &
                Ibar_select(flavor_j, 1, z)/z * pdfz

            ! delta(1-z), see above

            ! (1/(1-z))_+ contribution, boundary term see above
            if (z < 1._dp) then
                pdfone = fdist_one(ih, xi, mu, flavor_j, usebeam)

                getbeam_integrand = getbeam_integrand + &
                    1._dp/(1._dp - z) * &
                    (Ibar_select(flavor_j, 3, z)/z * pdfz &
                  - Ibar_select(flavor_j, 3, 1._dp) * pdfone )

                ! (log(1-z)/(1-z))_+ contribution, boundary term see above
                if (powAs > 1) then
                    getbeam_integrand = getbeam_integrand + &
                        log(1._dp-z)/(1._dp-z) * &
                        (Ibar_select(flavor_j, 4, z)/z * pdfz &
                      - Ibar_select(flavor_j, 4, 1._dp) * pdfone )
                endif

                ! ((log(1-z))^2/(1-z))_+ contribution, boundary term see above
                if(powAs > 2) then
                    getbeam_integrand = getbeam_integrand + &
                        log(1._dp-z)**2/(1._dp-z) * &
                        (Ibar_select(flavor_j, 5, z)/z * pdfz &
                      - Ibar_select(flavor_j, 5, 1._dp) * pdfone )
                endif

                ! ((log(1-z))^3/(1-z))_+ contribution, boundary term see above
                if(powAs > 3) then
                    getbeam_integrand = getbeam_integrand + &
                        log(1._dp-z)**3/(1._dp-z) * &
                        (Ibar_select(flavor_j, 6, z)/z * pdfz &
                      - Ibar_select(flavor_j, 6, 1._dp) * pdfone )

                endif

                ! ((log(1-z))^4/(1-z))_+ contribution, boundary term see above
                if(powAs > 4) then
                    getbeam_integrand = getbeam_integrand + &
                        log(1._dp-z)**4/(1._dp-z) * &
                        (Ibar_select(flavor_j, 7, z)/z * pdfz &
                      - Ibar_select(flavor_j, 7, 1._dp) * pdfone )

                endif

                ! ((log(1-z))^5/(1-z))_+ contribution, boundary term see above
                if(powAs > 5) then
                    getbeam_integrand = getbeam_integrand + &
                        log(1._dp-z)**5/(1._dp-z) * &
                        (Ibar_select(flavor_j, 8, z)/z * pdfz &
                      - Ibar_select(flavor_j, 8, 1._dp) * pdfone )

                endif

            endif

        enddo
    end function

!   subroutine check_Ibar_select()
!       implicit none
!       include 'ibar_checks2.f90'
!   end subroutine

    
    subroutine inithpl(x)
        use types
        use mod_hplog6
        implicit none
        real(dp), intent(in) :: x
        
        call hplog6(x,nw,Hc1,Hc2,Hc3,Hc4,Hc5,Hc6, &
                         Hr1,Hr2,Hr3,Hr4,Hr5,Hr6, &
                         Hi1,Hi2,Hi3,Hi4,Hi5,Hi6,n1,n2)

    end subroutine

    ! select appropriate Ibar_ functions based on i and j
    ! this has been exhaustively tested against MM code
    function Ibar_select(flavor_j, type_in, z)
        use constants
        implicit none

        integer, intent(in) :: flavor_j
        integer, intent(in) :: type_in
        real(dp), intent(in) :: z
        real(dp) :: Ibar_select

        real(dp) :: Li2, Li3, ReLi3

        Ibar_select = 0._dp

        if (flavor_i /= 0) then
            if (flavor_j == flavor_i) then
                ! Ibar_qq
#include "src/CuTe-MCFM/include_qq_3l.f90"
            elseif (flavor_j == 0) then
                ! Ibar_qg
#include "src/CuTe-MCFM/include_qg_3l.f90"
            elseif (flavor_j == -flavor_i) then
                ! Ibar_qbq
#include "src/CuTe-MCFM/include_qbq_3l.f90"
            elseif (flavor_i * flavor_j > 0) then
                ! Ibar_qpq
#include "src/CuTe-MCFM/include_qpq_3l.f90"
            elseif (flavor_i * flavor_j < 0) then
                ! Ibar_qpbq
#include "src/CuTe-MCFM/include_qpbq_3l.f90"
            else
                error stop "invalid flavor selection in Ibar_select"
            endif
        elseif (flavor_i == 0) then
            if (flavor_j == 0) then
                ! Ibar_gg
#include "src/CuTe-MCFM/include_gg_3l.f90"
            elseif (flavor_j /= 0) then
                ! Ibar_gq
#include "src/CuTe-MCFM/include_gq_3l.f90"
            else
                error stop "invalid flavor selection in Ibar_select"
            endif
        else
            error stop "invalid flavor selection in Ibar_select"
        endif

    end function

end module
