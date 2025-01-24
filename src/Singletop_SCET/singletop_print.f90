!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module SingletopPrint
    use MCFMStorage, only: PartStorage
    use types
    implicit none

    type(PartStorage), public, save :: finalSumNaive
    type(PartStorage), public, save :: finalSumFO

    public :: finalizeStorageFixedOrder

    public :: printcross_singletop

    public :: writeAllHistogramsSingletop

    private

    ! these correction factors in percent of twidth0 are insensitive to the
    ! top-quark mass and can be reused around 172.5 +/- 2
    real(dp), save :: twidth1
!$omp threadprivate(twidth1)
    real(dp), save :: twidth2
!$omp threadprivate(twidth2)

    real(dp), save :: twidth1_2mt
!$omp threadprivate(twidth1_2mt)
    real(dp), save :: twidth2_2mt
!$omp threadprivate(twidth2_2mt)

    real(dp), save :: twidth1_mtdiv2
!$omp threadprivate(twidth1_mtdiv2)
    real(dp), save :: twidth2_mtdiv2
!$omp threadprivate(twidth2_mtdiv2)

    ! like scalevarmult arraw in mod_Scalevar.f90
    ! multiplication factors: [2d0, 0.5d0, 2d0, 0.5d0, 1d0, 1d0, 2d0, 0.5d0]
    real(dp), save :: twidth1arr(1:8)
!$omp threadprivate(twidth1arr)

    real(dp), save :: twidth2arr(1:8)
!$omp threadprivate(twidth2arr)


    contains

    subroutine varianceToStdDev(stor)
        use types
        use SCET
        use Scalevar
        use PDFerrors
        use Integration
        use MCFMStorage
        implicit none

        type(PartStorage), intent(inout) :: stor

        integer :: nplotmax
        common/nplotmax/nplotmax

        integer :: k,m

        ! variance to standard deviation
        do k=1,nplotmax
            associate ( fin => stor%histCentral%histos(k) )
                if (fin%initialized()) then
                    fin%xxsq = sqrt(fin%xxsq)
                endif
            end associate

            if (doScalevar) then
                do m=1,maxscalevar
                    associate ( fin => stor%histScalevar(m)%histos(k) )
                        if (fin%initialized()) then
                            fin%xxsq = sqrt(fin%xxsq)
                        endif
                    end associate
                enddo
            endif

            if (maxPDFsets > 0) then
                do m=1,maxPDFsets
                    associate ( fin => stor%histPDFerrors(m)%histos(k) )
                        if (fin%initialized()) then
                            fin%xxsq = sqrt(fin%xxsq)
                        endif
                    end associate
                enddo
            endif

            do m=1,size(tcutarray)
                associate ( fin => stor%histTaucut(m)%histos(k) )
                    if (fin%initialized()) then
                        fin%xxsq = sqrt(fin%xxsq)
                    endif
                end associate
            enddo
        enddo

    end subroutine

    subroutine addData(histTarget, histSource, mult, tw1, tw2, tw1sq, naive)
        use types
        use SCET
        use Scalevar
        use PDFerrors
        use Integration
        use MCFMStorage
        use singletop2_scale_m, only: use_DDIS
        implicit none
        include 'kpart.f'
        include 'masses.f'

        type(PartStorage), intent(inout) :: histTarget
        type(PartStorage), intent(in) :: histSource
        real(dp), intent(in), optional :: mult
        logical, intent(in), optional :: tw1, tw2, tw1sq, naive

        integer :: nplotmax
        common/nplotmax/nplotmax

        integer :: k,l,m
        real(dp) :: mult_use

        integer :: truescalevar

        if (present(mult)) then
            mult_use = mult
        else
            mult_use = 1._dp
        endif

        !twidth0 = lotopdecaywidth(mt,0d0,wmass,0d0)

        ! if you modify this, also modify it below
        ! for scale variation, and afterwards
        if (present(tw1)) then
            if (tw1) then
                mult_use = -twidth1
            endif
        elseif (present(tw2)) then
            if (tw2) then
                mult_use = -twidth2
            endif
        elseif (present(tw1sq)) then
            if (tw1sq) then
                mult_use = (twidth1)**2
            endif
        elseif (present(naive)) then
            if (naive) then
                if (origKpart == ktota .or. origKpart == ksnlo) then
                    ! to correct to nlo decay width
                    mult_use = 1d0/(1d0+twidth1)
                elseif (origKpart == knnlo) then
                    mult_use = 1d0/(1d0+twidth1+twidth2)
                endif
            endif
        endif

        do k=1,nplotmax
            associate ( histo => histSource%histCentral%histos(k), &
                        fin => histTarget%histCentral%histos(k) )
                if ((fin%initialized() .eqv. .false.) .and. &
                    (histo%initialized() .eqv. .true.)) then
                    ! sum not initialized, but contribution
                    fin = histo
                    call fin%reset()
                endif
                if (histo%initialized()) then
                    do l=0,fin%nbins+1
                        if (histo%xx(l) /= 0._dp) then
                            fin%xx(l) = fin%xx(l) + &
                                histo%xx(l) / histo%xxsq(l) * mult_use
                            fin%xxsq(l) = fin%xxsq(l) + &
                                1._dp/histo%xxsq(l) * abs(mult_use)
                        endif
                    enddo
                endif
            end associate

            if (doScalevar) then
                if (use_DDIS) then
                    truescalevar = (maxscalevar-1)/2
                else
                    truescalevar = maxscalevar
                endif

                do m=1,maxscalevar
                    if (m <= truescalevar) then
                        ! set twidth1 based on twidth1arr(m)
                        if (present(tw1)) then
                            if (tw1) then
                                mult_use = -twidth1arr(m)
                            endif
                        elseif (present(tw2)) then
                            if (tw2) then
                                mult_use = -twidth2arr(m)
                            endif
                        elseif (present(tw1sq)) then
                            if (tw1sq) then
                                mult_use = (twidth1arr(m))**2
                            endif
                        elseif (present(naive)) then
                            if (naive) then
                                if (origKpart == ktota .or. origKpart == ksnlo) then
                                    ! to correct to nlo decay width
                                    mult_use = 1d0/(1d0+twidth1arr(m))
                                elseif (origKpart == knnlo) then
                                    mult_use = 1d0/(1d0+twidth1arr(m)+twidth2arr(m))
                                endif
                            endif
                        endif

                    elseif (m == truescalevar + 1) then
                        ! fixed mt
                        if (present(tw1)) then
                            if (tw1) then
                                mult_use = -twidth1
                            endif
                        elseif (present(tw2)) then
                            if (tw2) then
                                mult_use = -twidth2
                            endif
                        elseif (present(tw1sq)) then
                            if (tw1sq) then
                                mult_use = (twidth1)**2
                            endif
                        elseif (present(naive)) then
                            if (naive) then
                                if (origKpart == ktota .or. origKpart == ksnlo) then
                                    ! to correct to nlo decay width
                                    mult_use = 1d0/(1d0+twidth1)
                                elseif (origKpart == knnlo) then
                                    mult_use = 1d0/(1d0+twidth1+twidth2)
                                endif
                            endif
                        endif
                    elseif (m > truescalevar + 1) then
                        ! set twidth based on twidth1arr(m - truescalevar - 1)

                        if (present(tw1)) then
                            if (tw1) then
                                mult_use = -twidth1arr(m-truescalevar-1)
                            endif
                        elseif (present(tw2)) then
                            if (tw2) then
                                mult_use = -twidth2arr(m-truescalevar-1)
                            endif
                        elseif (present(tw1sq)) then
                            if (tw1sq) then
                                mult_use = (twidth1arr(m-truescalevar-1))**2
                            endif
                        elseif (present(naive)) then
                            if (naive) then
                                if (origKpart == ktota .or. origKpart == ksnlo) then
                                    ! to correct to nlo decay width
                                    mult_use = 1d0/(1d0+twidth1arr(m-truescalevar-1))
                                elseif (origKpart == knnlo) then
                                    mult_use = 1d0/(1d0+twidth1arr(m-truescalevar-1) &
                                                       +twidth2arr(m-truescalevar-1))
                                endif
                            endif
                        endif
                    endif
                    associate (histo => histSource%histScalevar(m)%histos(k), &
                            fin => histTarget%histScalevar(m)%histos(k) )
                        if ((fin%initialized() .eqv. .false.) .and. (histo%initialized() .eqv. .true.)) then
                            ! sum not initialized, but contribution
                            fin = histo
                        endif
                        if (histo%initialized()) then
                            do l=0,fin%nbins+1
                                if (histo%xx(l) /= 0._dp) then
                                    fin%xx(l) = fin%xx(l) + &
                                        histo%xx(l) / histo%xxsq(l) * mult_use
                                    fin%xxsq(l) = fin%xxsq(l) + &
                                        1._dp/histo%xxsq(l) * abs(mult_use)
                                endif
                            enddo
                        endif
                    end associate
                enddo

                ! reset mult_use factors
            endif

            if (maxPDFsets > 0) then
                do m=1,maxPDFsets
                    associate (histo => histSource%histPDFerrors(m)%histos(k), &
                            fin => histTarget%histPDFerrors(m)%histos(k) )
                        if ((fin%initialized() .eqv. .false.) .and. (histo%initialized() .eqv. .true.)) then
                            ! sum not initialized, but contribution
                            fin = histo
                        endif
                        if (histo%initialized()) then
                            do l=0,fin%nbins+1
                                if (histo%xx(l) /= 0._dp) then
                                    fin%xx(l) = fin%xx(l) + &
                                        histo%xx(l) / histo%xxsq(l) * mult_use
                                    fin%xxsq(l) = fin%xxsq(l) + &
                                        1._dp/histo%xxsq(l) * abs(mult_use)
                                endif
                            enddo
                        endif
                    end associate
                enddo
            endif

            do m=1,size(tcutarray)
                associate (histo => histSource%histTaucut(m)%histos(k), &
                        fin => histTarget%histTaucut(m)%histos(k) )
                   if ((fin%initialized() .eqv. .false.) .and. (histo%initialized() .eqv. .true.)) then
                       ! sum not initialized, but contribution
                       fin = histo
                   endif
                    if (histo%initialized()) then
                        do l=0,fin%nbins+1
                            if (histo%xx(l) /= 0._dp) then
                                fin%xx(l) = fin%xx(l) + &
                                    histo%xx(l) / histo%xxsq(l) * mult_use
                                fin%xxsq(l) = fin%xxsq(l) + &
                                    1._dp/histo%xxsq(l) * abs(mult_use)
                            endif
                        enddo
                    endif
                end associate
            enddo
        enddo

    end subroutine

    ! like our normal finalizeStorage, but combines pieces in production
    ! and decay with the width in a fixed-order way
    subroutine finalizeStorageFixedOrder()
        use types
        use SCET
        use Scalevar
        use PDFerrors
        use Integration
        use MCFMStorage, only: iterationStorage
        use singletop2_nnlo, only: singletop2_nnlo_enable_heavy_decay, &
            singletop2_nnlo_enable_lxh, singletop2_nnlo_fully_inclusive, &
            singletop2_nnlo_enable_lxd, singletop2_nnlo_enable_hxd
        use topwidth
        implicit none
        include 'mpicommon.f'
        include 'kpart.f'
        include 'masses.f'

        ! process iterationStorage for presentation
        integer :: nplotmax
        common/nplotmax/nplotmax
        integer :: i,j
        integer :: forder

        integer :: iused, jused

        ! find a storage block that has been used to initialize finalSum
        ! and also count number of contributions
        do i=1,maxParts
            do j=1,maxIps
                if (iterationStorage(i,j)%used .eqv. .true.) then
                    iused = i
                    jused = j
                endif
            enddo
        enddo

        finalSumNaive = iterationStorage(iused,jused)
        call finalSumNaive%reset

        ! strict fixed-order expansion
        finalSumFO = finalSumNaive
        call finalSumFO%reset


        ! set up top decay width shared variables

        ! renormalized with mu=mt
        twidth1 = nloratiotopdecay(mt, 0d0, wmass, 0d0, mt)
        twidth1_2mt = nloratiotopdecay(mt, 0d0, wmass, 0d0, 2._dp*mt)
        twidth1_mtdiv2 = nloratiotopdecay(mt, 0d0, wmass, 0d0, mt/2._dp)

        twidth2 = nnlotopdecay(mt, wmass, mt)
        twidth2_2mt = nnlotopdecay(mt, wmass, 2._dp*mt)
        twidth2_mtdiv2 = nnlotopdecay(mt, wmass, mt/2._dp)

        twidth1arr(1:8) = [twidth1_2mt, twidth1_mtdiv2, &
            twidth1_2mt, twidth1_mtdiv2, &
            twidth1, twidth1, &
            twidth1_2mt, twidth1_mtdiv2]

        twidth2arr(1:8) = [twidth2_2mt, twidth2_mtdiv2, &
            twidth2_2mt, twidth2_mtdiv2, &
            twidth2, twidth2, &
            twidth2_2mt, twidth2_mtdiv2]


        ! for strict fixed-order expansion
        ! only for ktota and knnlo the expansion below makes sense
        if (origKpart == ktota .or. origKpart == ksnlo) then
            forder = 1
        elseif (origKpart == knnlo) then
            forder = 2
        else
            forder = 0
        endif

        if (iterationStorage(lord,1)%used .eqv. .true.) then
            call addData(finalSumNaive, iterationStorage(lord,1), naive=.true.)
            call addData(finalSumFO, iterationStorage(lord,1))
        endif

        do j=1,7
            if (iterationStorage(snloAbove,j)%used .eqv. .true.) then
                if (forder > 0) then
                    call addData(finalSumFO, iterationStorage(snloAbove,j))
                endif
                if (forder > 1 .and. (.not. singletop2_nnlo_fully_inclusive)) then
                    call addData(finalSumFO, iterationStorage(snloAbove,j), tw1=.true.)
                endif
            endif
        enddo

        do j=1,5
            if (iterationStorage(snloBelow,j)%used .eqv. .true.) then
                if (forder > 0) then
                    call addData(finalSumFO, iterationStorage(snloBelow,j))
                endif
                if (forder > 1 .and. (.not. singletop2_nnlo_fully_inclusive)) then
                    call addData(finalSumFO, iterationStorage(snloBelow,j), tw1=.true.)
                endif
            endif
        enddo

        ! light line nlocoeff
        if (iterationStorage(nloReal,1)%used .eqv. .true.) then
            call addData(finalSumNaive,iterationStorage(nloReal,1), naive=.true.)
            if (forder > 0) then
                call addData(finalSumFO, iterationStorage(nloReal,1))
            endif
            if (forder > 1 .and. (.not. singletop2_nnlo_fully_inclusive)) then
                call addData(finalSumFO, iterationStorage(nloReal,1), tw1=.true.)
            endif
        endif
        if (iterationStorage(nloVirt,1)%used .eqv. .true.) then
            call addData(finalSumNaive,iterationStorage(nloVirt,1), naive=.true.)
            if (forder > 0) then
                call addData(finalSumFO, iterationStorage(nloVirt,1))
            endif
            
            if (forder > 1 .and. (.not. singletop2_nnlo_fully_inclusive)) then
                call addData(finalSumFO, iterationStorage(nloVirt,1), tw1=.true.)
            endif
        endif

        ! heavy line nlocoeff
        if (iterationStorage(nloReal,2)%used .eqv. .true.) then
            call addData(finalSumNaive, iterationStorage(nloReal,2), naive=.true.)
            if (forder > 0) then
                call addData(finalSumFO, iterationStorage(nloReal,2))
            endif

            if (forder > 1 .and. (.not. singletop2_nnlo_fully_inclusive)) then
                call addData(finalSumFO, iterationStorage(nloReal,2), tw1=.true.)
            endif
        endif
        if (iterationStorage(nloVirt,2)%used .eqv. .true.) then
            call addData(finalSumNaive, iterationStorage(nloVirt,2), naive=.true.)
            if (forder > 0) then
                call addData(finalSumFO, iterationStorage(nloVirt,2))
            endif

            if (forder > 1 .and. (.not. singletop2_nnlo_fully_inclusive)) then
                call addData(finalSumFO, iterationStorage(nloVirt,2), tw1=.true.)
            endif
        endif

        ! decay nlocoeff
        if (iterationStorage(nloReal,3)%used .eqv. .true.) then
            call addData(finalSumNaive, iterationStorage(nloReal,3), naive=.true.)
            if (forder > 0) then
                call addData(finalSumFO, iterationStorage(nloReal,3))
            endif

            if (forder > 1) then
                call addData(finalSumFO, iterationStorage(nloReal,3), tw1=.true.)
            endif
        endif
        if (iterationStorage(nloVirt,3)%used .eqv. .true.) then
            call addData(finalSumNaive, iterationStorage(nloVirt,3), naive=.true.)
            if (forder > 0) then
                call addData(finalSumFO, iterationStorage(nloVirt,3))
            endif

            if (forder > 1) then
                call addData(finalSumFO, iterationStorage(nloVirt,3), tw1=.true.)
            endif
        endif

        if ((.not. singletop2_nnlo_fully_inclusive) .and. &
                    singletop2_nnlo_enable_heavy_decay ) then
            ! additional pieces for decay; comment out if only interested in
            ! production
            if (forder > 0) then
                if (iterationStorage(lord,1)%used .eqv. .true.) then
                    call addData(finalSumFO, iterationStorage(lord,1), tw1=.true.)
                endif
            endif
        endif

        !! NNLO pieces

        !!! nnlo production pieces dsigma2 * dgamma0

        do j=1,4
            if (iterationStorage(nnloBelow,j)%used .eqv. .true.) then
                call addData(finalSumNaive, iterationStorage(nnloBelow,j), naive=.true.)
                if (forder > 1) then
                    call addData(finalSumFO, iterationStorage(nnloBelow,j))
                endif
            endif
        enddo

        do j=1,4
            if (iterationStorage(nnloVirtAbove,j)%used .eqv. .true.) then
                call addData(finalSumNaive, iterationStorage(nnloVirtAbove,j), naive=.true.)
                if (forder > 1) then
                    call addData(finalSumFO, iterationStorage(nnloVirtAbove,j))
                endif
            endif
        enddo

        ! real above for heavy production
        do j=1,6
            if (iterationStorage(nnloRealAbove,j)%used .eqv. .true.) then
                call addData(finalSumNaive, iterationStorage(nnloRealAbove,j), naive=.true.) 
                if (forder > 1) then
                    call addData(finalSumFO, iterationStorage(nnloRealAbove,j))
                endif
            endif
        enddo

        ! real above for light production
#define SINGLETOP_NNLO_REALCHANNELS 12
#if SINGLETOP_NNLO_REALCHANNELS == 66
        do j=8,73
#else
        do j=8,19
#endif
            if (iterationStorage(nnloRealAbove,j)%used .eqv. .true.) then
                call addData(finalSumNaive, iterationStorage(nnloRealAbove,j), naive=.true.)
                if (forder > 1) then
                    call addData(finalSumFO, iterationStorage(nnloRealAbove,j))
                endif
            endif
        enddo

        ! proper heavy x light interference
        if (forder > 1 .and. singletop2_nnlo_enable_lxh) then
            if (iterationStorage(nloVirt,4)%used .eqv. .true.) then
                call addData(finalSumNaive, iterationStorage(nloVirt,4), naive=.true.)
                call addData(finalSumFO, iterationStorage(nloVirt,4))
            endif
            if (iterationStorage(nloVirt,5)%used .eqv. .true.) then
                call addData(finalSumNaive, iterationStorage(nloVirt,5), naive=.true.)
                call addData(finalSumFO, iterationStorage(nloVirt,5))
            endif

            if (iterationStorage(nloReal,4)%used .eqv. .true.) then
                call addData(finalSumNaive, iterationStorage(nloReal,4), naive=.true.)
                call addData(finalSumFO, iterationStorage(nloReal,4))
            endif
            if (iterationStorage(nloReal,5)%used .eqv. .true.) then
                call addData(finalSumNaive, iterationStorage(nloReal,5), naive=.true.)
                call addData(finalSumFO, iterationStorage(nloReal,5))
            endif
        endif


        if (.not. singletop2_nnlo_fully_inclusive) then
            if (forder > 1 .and. singletop2_nnlo_enable_lxd) then
                if (iterationStorage(nloVirt,6)%used .eqv. .true.) then
                    call addData(finalSumNaive, iterationStorage(nloVirt,6), naive=.true.)
                    call addData(finalSumFO, iterationStorage(nloVirt,6))
                endif
                if (iterationStorage(nloVirt,7)%used .eqv. .true.) then
                    call addData(finalSumNaive, iterationStorage(nloVirt,7), naive=.true.)
                    call addData(finalSumFO, iterationStorage(nloVirt,7))
                endif

                if (iterationStorage(nloReal,6)%used .eqv. .true.) then
                    call addData(finalSumNaive, iterationStorage(nloReal,6), naive=.true.)
                    call addData(finalSumFO, iterationStorage(nloReal,6))
                endif
                if (iterationStorage(nloReal,7)%used .eqv. .true.) then
                    call addData(finalSumNaive, iterationStorage(nloReal,7), naive=.true.)
                    call addData(finalSumFO, iterationStorage(nloReal,7))
                endif
            endif

            if (forder > 1 .and. singletop2_nnlo_enable_hxd) then
                if (iterationStorage(nloVirt,8)%used .eqv. .true.) then
                    call addData(finalSumNaive, iterationStorage(nloVirt,8), naive=.true.)
                    call addData(finalSumFO, iterationStorage(nloVirt,8))
                endif
                if (iterationStorage(nloVirt,9)%used .eqv. .true.) then
                    call addData(finalSumNaive, iterationStorage(nloVirt,9), naive=.true.)
                    call addData(finalSumFO, iterationStorage(nloVirt,9))
                endif

                if (iterationStorage(nloReal,8)%used .eqv. .true.) then
                    call addData(finalSumNaive, iterationStorage(nloReal,8), naive=.true.)
                    call addData(finalSumFO, iterationStorage(nloReal,8))
                endif
                if (iterationStorage(nloReal,9)%used .eqv. .true.) then
                    call addData(finalSumNaive, iterationStorage(nloReal,9), naive=.true.)
                    call addData(finalSumFO, iterationStorage(nloReal,9))
                endif

            endif
        endif


        ! nnlo decay, dsigma0 * dgamma2
        if (iterationStorage(nnloBelow,5)%used .eqv. .true.) then
            call addData(finalSumNaive, iterationStorage(nnloBelow,5), naive=.true.)
            if (forder > 1) then
                call addData(finalSumFO, iterationStorage(nnloBelow,5))
            endif
        endif
        if (iterationStorage(nnloVirtAbove,5)%used .eqv. .true.) then
            call addData(finalSumNaive, iterationStorage(nnloVirtAbove,5), naive=.true.)
            if (forder > 1) then
                call addData(finalSumFO, iterationStorage(nnloVirtAbove,5))
            endif
        endif
        if (iterationStorage(nnloRealAbove,7)%used .eqv. .true.) then
            call addData(finalSumNaive, iterationStorage(nnloRealAbove,7), naive=.true.)
            if (forder > 1) then
                call addData(finalSumFO, iterationStorage(nnloRealAbove,7))
            endif
        endif

        if ((.not. singletop2_nnlo_fully_inclusive) .and. &
                singletop2_nnlo_enable_heavy_decay ) then
            ! additional pieces for decay; comment out if only interested in
            ! production
            if (forder > 1) then
                if (iterationStorage(lord,1)%used .eqv. .true.) then
                    call addData(finalSumFO, iterationStorage(lord,1), tw2=.true.)
                    call addData(finalSumFO, iterationStorage(lord,1), tw1sq=.true.)
                endif
            endif
        endif


        call varianceToStdDev(finalSumNaive)
        call varianceToStdDev(finalSumFO)

    end subroutine

    subroutine printcross_singletop()
        use PDFerrors
        use MCFMPrint, only: formatCross, printScaleUncertainties
        use topwidth
        use LHAPDF
        implicit none
        include 'kpart.f'
        include 'masses.f'
        real(dp) :: centralCross, centralError, cross, error
        integer, allocatable :: centralBins(:)
        integer :: j, accum
        real(dp) :: twidth0

        allocate(centralBins(numPDFsets))

            ! the first set is our MC central value, saved separately
        accum = 1
        do j=2,numPDFsets
            if (doPDFerrors) then
                accum = accum + lhapdf_number(trim(PDFnames(j-1)))
            else
                accum = accum + 1
            endif
            
            centralBins(j) = accum
        enddo

        centralBins(2:) = centralBins(2:) - 1

        if (numPDFsets > 1) then
            write (*,'(A)') "=== Printing central cross section values for all PDF sets ==="
        endif

        write (*,'(A,A,A,I3,A)') "=== Result for PDF set ", trim(PDFnames(1)), " member ", &
            PDFmembers(1), " ==="

        centralCross = finalSumNaive%histCentral%histos(1)%xx(1)
        centralError = finalSumNaive%histCentral%histos(1)%xxsq(1)
        write (6,'(A,A)') "Naive result (sum) is ", formatCross(centralCross, centralError)

        twidth0 = lotopdecaywidth(mt,0d0,wmass,0d0)

        if (origKpart == ktota .or. origKpart == ksnlo) then
            write (*,'(A,F6.4)') "Using NLO decay width ", twidth0*(1d0+twidth1)
        elseif (origKpart == knnlo) then
            write (*,'(A,F6.4)') "Using NNLO decay width ", twidth0*(1d0+twidth1+twidth2)
        else
            write (*,'(A,F6.4)') "Using LO decay width ", twidth0
        endif

        if (origKpart == ktota .or. origKpart == ksnlo .or. origKpart == knnlo) then
            centralCross = finalSumFO%histCentral%histos(1)%xx(1)
            centralError = finalSumFO%histCentral%histos(1)%xxsq(1)
            write (6,'(A,A)') "Strict fixed-order expanded result is ", &
                formatCross(centralCross, centralError)
        endif

        write (*,*) "=== Scale uncertainties for fixed-order expansion ==="
        call printScaleUncertainties(finalSumFO)

        write (*,*) "=== Scale uncertainties for naive expansion ==="
        call printScaleUncertainties(finalSumNaive)

        do j=2,numPDFsets
            cross = finalSumNaive%histPDFerrors(centralBins(j))%histos(1)%xx(1)
            error = finalSumNaive%histPDFerrors(centralBins(j))%histos(1)%xxsq(1)

            write (*,*) ""
            write (6,'(A,A,A,A)') "Naive difference from ", trim(PDFnames(j)), &
                " to first PDF set is ", formatCross(cross,error)
        enddo

        if (origKpart == klord .or. origKpart == ktota .or. origKpart == ksnlo .or. origKpart == knnlo) then
            do j=2,numPDFsets
                cross = finalSumFO%histPDFerrors(centralBins(j))%histos(1)%xx(1)
                error = finalSumFO%histPDFerrors(centralBins(j))%histos(1)%xxsq(1)

                write (*,*) ""
                write (6,'(A,A,A,A)') "Fixed-order difference from ", trim(PDFnames(j)), &
                    " to first PDF set is ", formatCross(cross,error)
            enddo
        endif

    end subroutine

    subroutine writeAllHistogramsSingletop()
        use parseinput
        use Superhisto
        use PDFerrors
        use MCFMPrint
        implicit none

        logical :: writetxt, writetop
        call cfg_get(cfg, "histogram%writetxt", writetxt)
        call cfg_get(cfg, "histogram%writetop", writetop)

        if (writetxt) then
            call writeHistograms(shwrite,"fo.txt",finalSumFO)

            call writeHistograms(shwrite,"naive.txt",finalSumNaive)

            if (doPDFerrors) then
                call writeHistogramsPDFerrors(shwritepdf,"fo.txt",finalSumFO)
            endif
        endif
        if (writetop) then
            call writeHistograms(shwritetop,"fo.txt",finalSumFO)
        endif

    end subroutine



end module
