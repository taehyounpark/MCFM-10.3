!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
! The idea is to use the total cross information, its uncertainty
! and possibly taucut fitting and compare it with stored benchmarks.

! All benchmarks are supposed to run with the default 0.2 per mille precisiongoal

! Benchmarking also allows us to determine better default integration settings
! for some processes.
    

module MCFMBenchmark
    use types
    use MCFMStorage
    use MCFMPrint, only : chisqMax
    use m_config
    use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
        stdout=>output_unit, stderr=>error_unit

    public :: setupBenchmark, comparisonCode
    private

    real(dp), save :: storedCross, storedError
    type(CFG_T), save :: cfgBench

    contains

    ! routine to set up cuts for a finite cross section
    subroutine setupBenchmark()
        use m_config
        use parseinput
        implicit none
        include 'nproc.f'
        include 'kpart.f'

        character(len=cfg_string_len) :: benchmark

        if (cfg_var_configadded(cfg, "extra%benchmark")) then
            call cfg_get(cfg, "extra%benchmark", benchmark)
        else
            error stop "setupBenchmark called but extra%benchmark not set"
        endif

        storedCross = 0._dp
        storedError = 0._dp

        ! To add a new benchmark, add a case and use cfg_set appropriately.
        ! Keep storedCross and storedError as 0 for the first run.
        select case(benchmark)
            case ("1_WPlus_nlo")
                storedCross = 4221443.97329230_dp 
                storedError = 2656.87500656489_dp
            case ("1_W_nnlo")
                call cfg_set(cfg, "general%part", "nnlo")

                storedCross = 4.20965e6_dp
                storedError = 0.55583E+04_dp
            case ("6_WMinus_nnlo")
                call cfg_set(cfg, "general%nproc", 6)

                storedCross = 3.27474e6_dp
                storedError = 0.36388E+04_dp
            case ("31_Z_nnlo")
                call cfg_set(cfg, "general%part", "nnlo")
                call cfg_set(cfg, "general%nproc", 31)
                call cfg_set(cfg, "masscuts%m34min", 40._dp)

                storedCross = 875.799e3_dp
                storedError = 0.85120e3_dp
            case ("61_WW_nlo")
                call cfg_set(cfg, "general%nproc", 61)
                call cfg_set(cfg, "scales%dynamicscale", "m(3456)")

                storedCross = 435.095081136187_dp
                storedError = 0.869225231444390_dp
            case ("71_WZ_nlo")
                call cfg_set(cfg, "general%nproc", 71)
                call cfg_set(cfg, "scales%dynamicscale", "m(3456)")
                call cfg_set(cfg, "masscuts%m56min", 40._dp)

                storedCross = 26.9180118435005_dp
                storedError = 5.353466059414824E-002_dp
            case ("81_ZZ_nlo")
                call cfg_set(cfg, "general%nproc", 81)
                call cfg_set(cfg, "scales%dynamicscale", "m(3456)")
                call cfg_set(cfg, "masscuts%m34min", 40._dp)
                call cfg_set(cfg, "masscuts%m56min", 40._dp)

                storedCross = 13.4128496832762_dp
                storedError = 1.771228715518188E-002_dp
            case ("91_WPlusH_nnlo")
                call cfg_set(cfg, "general%part", "nnlo")
                call cfg_set(cfg, "general%nproc", 91)
                call cfg_set(cfg, "scales%dynamicscale", "m(3456)")

                storedCross = 2.26184_dp
                storedError = 0.41537E-02_dp
            case ("96_WMinusH_nnlo")
                call cfg_set(cfg, "general%part", "nnlo")
                call cfg_set(cfg, "general%nproc", 96)
                call cfg_set(cfg, "scales%dynamicscale", "m(3456)")

                storedCross = 1.52580_dp
                storedError = 0.27710E-02_dp
            case ("110_ZH_nnlo")
                call cfg_set(cfg, "general%part", "nnlo")
                call cfg_set(cfg, "general%nproc", 110)
                call cfg_set(cfg, "scales%dynamicscale", "m(3456)")
                call cfg_set(cfg, "masscuts%m34min", 40._dp)

                storedCross = 0.841771_dp
                storedError = 0.10933E-02_dp
            case ("112_H_nnlo")
                call cfg_set(cfg, "general%part", "nnlo")

                storedCross = 1.87237e3_dp
                storedError = 0.21994E1_dp
            case ("285_GamGam_nnlo")
                call cfg_set(cfg, "general%part", "nnlo")

                storedCross = 43.5362e3_dp
                storedError = 0.79486e2_dp
            case ("300_ZGam_nnlo")
                call cfg_set(cfg, "general%part", "nnlo")
                call cfg_set(cfg, "general%nproc", 300)
                call cfg_set(cfg, "scales%dynamicscale", "m(345)")
                call cfg_set(cfg, "masscuts%m34min", 40._dp)
                call cfg_set(cfg, "photon%Rgalmin", 0.3_dp)

                storedCross = 525.522_dp
                storedError = 1.0391_dp
            case default
                error stop "ERROR: No setup for this benchmark"
        end select

    end subroutine

    ! failed check: comparisonCode = 1
    ! no benchmark available: comparisonCode = 2
    ! successful comparison, or not final result: comparisonCode = 0
    ! undefined: comparisonCode = 100

    function comparisonCode(stage)
        use m_config
        use parseinput
        implicit none
        include 'nproc.f'
        include 'kpart.f'
        include 'mpicommon.f'

        integer, intent(in) :: stage
        integer :: comparisonCode
        integer :: outputUnit

        real(dp) :: cross, error
        ! The check for the difference will be performed with respect
        ! to the integration error times this value of errorMargin
        real(dp), parameter :: errorMargin = 2._dp

        character(len=cfg_string_len) :: benchmark

        if (cfg_var_configadded(cfg, "extra%benchmark")) then
            call cfg_get(cfg, "extra%benchmark", benchmark)
        else
            error stop "comparisonCode called but extra%benchmark not set"
        endif

        comparisonCode = 100
        cross = finalSum%histCentral%histos(1)%xx(1)
        error = finalSum%histCentral%histos(1)%xxsq(1)

        if (storedCross == 0._dp .and. storedError == 0._dp) then
            write (stderr,*) ""
            write (stderr,*) "No benchmark comparison available for nproc = ", nproc
            write (stderr,*) "Calculated cross: ", cross
            write (stderr,*) "Calculated error: ", error
            if (chisqMax() > 1._dp) then
                write(stderr,*) "WARNING: maximum chi^2/it = ", chisqMax()
            endif
            write (stderr,*) ""
            comparisonCode = 2
            return
        endif

        if (stage == 1) then
            write (stdout,*) "INFO: Intermediate benchmark information"
            outputUnit = stdout
            comparisonCode = 0
        else
            if ( abs(storedCross - cross) > errorMargin*storedError ) then
                write (stderr,*) "ERROR: Found unusually large difference"
                outputUnit = stderr
                comparisonCode = 1
            else
                write (stdout,*) "INFO: Successful benchmark comparison"
                outputUnit = stdout
                comparisonCode = 0
            endif
        endif

        write (outputUnit,*) "Cross calculated, stored:", cross, storedCross
        write (outputUnit,*) "Difference relative to results: ", &
            abs(cross-storedCross)/abs(cross+storedCross)/2._dp
        write (outputUnit,*) "Difference relative to errors: ", &
            abs(cross-storedCross)/(error+storedError)
        write (outputUnit,*) "Error calculated, stored:", error, storedError
        ! for now just use it as a warning
        if (chisqMax() > 1._dp) then
            write(stderr,*) "WARNING: maximum chi^2/it = ", chisqMax()
        endif
        write (outputUnit,*) ""

    end function

end module
