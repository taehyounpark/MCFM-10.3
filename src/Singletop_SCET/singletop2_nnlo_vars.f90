!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module singletop2_nnlo_vars
    implicit none
    logical, public, save :: singletop2_nnlo_enable_light = .true.
    logical, public, save :: singletop2_nnlo_enable_heavy_prod = .true.
    logical, public, save :: singletop2_nnlo_enable_heavy_decay = .true.

    ! flags for interferences
    logical, public, save :: singletop2_nnlo_enable_lxh = .true.
    logical, public, save :: singletop2_nnlo_enable_lxd = .true.
    logical, public, save :: singletop2_nnlo_enable_hxd = .true.

    logical, public, save :: singletop2_nnlo_fully_inclusive = .false.

    integer, public, parameter :: max_corr_on_beam = 2
    integer, public, parameter :: max_bcontrib = 5

    ! replaces ipsgen 1,2 or 3
    integer, public :: currentContrib = 0

    ! bcontribs are
    ! 1: main b
    ! 2: additional b
    ! 3: additional b~
    ! 4: additional pair b b~

    integer, public, save :: corr_on_beam = 1
!$omp threadprivate(corr_on_beam)

    integer, public, save :: b_contrib = 1
!$omp threadprivate(b_contrib)

    ! which tau phase space region is selected 
    integer, public, save :: taups = 0

    integer, public, save :: maxbeams, beams_enabled(2)

    integer, public, parameter :: quarkChannel = 1, gluonChannel = 2
    integer, public, save :: partons_enabled

    logical, public, save :: maskb1(-5:5), maskb2(-5:5)

    logical, public, save :: usemask = .false.

    public :: setup_singletop

    private

    contains

    ! this routine maps the triplet (nproc, currentPart, currentIps) to
    ! * currentContrib (light, heavy, decay)
    ! * partons enabled (mask, partons_enabled)
    ! * corr_on_beam (beams_enabled)

    ! some information is redundant, like mask and partons_enabled
    ! but we want different granularity in different pieces

    ! we do not want to depend at all anymore on "ipsgen" in any of the
    ! singletop codes, but just rely on the derived variables above
    subroutine setup_singletop(ips)
!        use MCFMStorage, only: currentPart, currentIps
        use Integration
        implicit none
        include 'nproc.f'
        include 'kpart.f'
        include 'mpicommon.f'

        integer, intent(in) :: ips

        maskb1(:) = .false.
        maskb2(:) = .false.

        if (ips == 0) then
            ! this happens when chooser is called in in parseinput
!           if (rank == 0) then
!               write (*,*) "WARNING: setup_singletop called with ips = 0"
!           endif
            return
        endif

        if (nproc == 1650) then
            if (kpart == klord) then
                if (any(ips == [1,2])) then
                     currentContrib = 1
                     maxbeams = 1
                     beams_enabled(1) = 1
                     usemask = .true.
                     maskb1(:) = .true.
                     maskb2(5) = .true.
                     partons_enabled = gluonChannel + quarkChannel
                     taups = ips
                     return
                 elseif (any(ips == [3,4])) then
                     currentContrib = 1
                     maxbeams = 1
                     beams_enabled(1) = 2
                     usemask = .true.
                     maskb1(5) = .true.
                     maskb2(:) = .true.
                     partons_enabled = gluonChannel + quarkChannel
                     taups = ips - 2
                     return
                elseif (any(ips == [5])) then
                     currentContrib = 2
                     maxbeams = 1
                     beams_enabled(1) = 1
                     usemask = .true.
                     maskb1([0,5]) = .true.
                     maskb2(:) = .true.
                     partons_enabled = gluonChannel + quarkChannel
                     taups = 0
                     return
                elseif (any(ips == [6])) then
                     currentContrib = 2
                     maxbeams = 1
                     beams_enabled(1) = 2
                     usemask = .false.
                     maskb1(:) = .true.
                     maskb2([0,5]) = .true.
                     partons_enabled = gluonChannel + quarkChannel
                     taups = 0
                     return
                 elseif (ips == 7) then
                     currentContrib = 3
                     maxbeams = 2
                     beams_enabled(1:2) = [1,2]
                     usemask = .true.
                     maskb1([-3,-1,2,4,5]) = .true.
                     maskb2([-3,-1,2,4,5]) = .true.
                     partons_enabled = gluonChannel + quarkChannel
                     taups = 0
                     return
                 endif

            elseif (origKpart == knnlo .and. kpart == kvirt) then
                if (ips == 1) then
                    ! light line RV, *b
                    currentContrib = 1
                    maxbeams = 1
                    beams_enabled(1) = 1
                    usemask = .false.
                    maskb1(:) = .true.
                    maskb2(5) = .true.
                    partons_enabled = quarkChannel + gluonChannel
                    taups = 0
                    return
                elseif (ips == 2) then
                    currentContrib = 1
                    maxbeams = 1
                    beams_enabled(1) = 2
                    usemask = .false.
                    maskb1(5) = .true.
                    maskb2(:) = .true.
                    partons_enabled = quarkChannel + gluonChannel
                    taups = 0
                    return
                elseif (ips == 3) then
                    currentContrib = 2
                    maxbeams = 1
                    beams_enabled(1) = 1
                    usemask = .true.
                    maskb1(:) = .true.
                    maskb2([-3,-1,2,4]) = .true.
                    partons_enabled = quarkChannel + gluonChannel
                    taups = 0
                    return
                elseif (ips == 4) then
                    currentContrib = 2
                    maxbeams = 1
                    beams_enabled(1) = 2
                    usemask = .true.
                    maskb1([-3,-1,2,4]) = .true.
                    maskb2(:) = .true.
                    partons_enabled = quarkChannel + gluonChannel
                    taups = 0
                    return
                elseif (ips == 5) then
                    currentContrib = 3
                    maxbeams = 2
                    beams_enabled(1:2) = [1,2]
                    usemask = .true.
                    maskb1([-3,-1,2,4,5]) = .true.
                    maskb2([-3,-1,2,4,5]) = .true.
                    partons_enabled = quarkChannel + gluonChannel
                    taups = 0
                    return
                endif
            elseif (kpart == kvirt) then
                if (ips == 1) then
                    ! 165 nlo, virt, all together
                    currentContrib = 1
                    maxbeams = 2
                    beams_enabled(:) = [1,2]
                    usemask = .false.
                    maskb1(:) = .true.
                    maskb2(:) = .true.
                    partons_enabled = quarkChannel + gluonChannel
                    taups = 0
                    return
                elseif (ips == 2) then
                    ! 165 nlo, virt, all together
                    currentContrib = 2
                    maxbeams = 2
                    beams_enabled(:) = [1,2]
                    usemask = .false.
                    maskb1(:) = .true.
                    maskb2(:) = .true.
                    partons_enabled = quarkChannel + gluonChannel
                    taups = 0
                    return
                elseif (ips == 3) then
                    currentContrib = 3
                    maxbeams = 2
                    beams_enabled(:) = [1,2]
                    usemask = .false.
                    maskb1(:) = .true.
                    maskb2(:) = .true.
                    partons_enabled = quarkChannel + gluonChannel
                    taups = 0
                    return
                endif
            elseif (origKpart == knnlo .and. kpart == kreal) then
#define SINGLETOP_NNLO_REALCHANNELS 12
#define OFFSET 7
#if SINGLETOP_NNLO_REALCHANNELS == 66
                if (any(ips == [(i, i=1+OFFSET,33+OFFSET, 1)])) then
                    ! light line, beam 1
                    currentContrib = 1
                    maxbeams = 1
                    beams_enabled(1) = 1

                    usemask = .true.
                    maskb2(5) = .true.

                    if (ips < 12+OFFSET) then
                        maskb1(ips-6-OFFSET) = .true.
                        if (ips-6 == OFFSET) then
                            partons_enabled = gluonChannel
                        else
                            partons_enabled = quarkChannel
                        endif
                        taups = 1
                    elseif (ips < 23+OFFSET) then
                        maskb1(ips-6-11-OFFSET) = .true.
                        if (ips-6-11 == OFFSET) then
                            partons_enabled = gluonChannel
                        else
                            partons_enabled = quarkChannel
                        endif
                        taups = 2
                    else
                        maskb1(ips-6-22-OFFSET) = .true.
                        if (ips-6-22 == OFFSET) then
                            partons_enabled = gluonChannel
                        else
                            partons_enabled = quarkChannel
                        endif
                        taups = 3
                    endif

                    return
                elseif (any(ips == [(i, i=34+OFFSET,66+OFFSET, 1)])) then
                    ! light line beam 2
                    currentContrib = 1
                    maxbeams = 1
                    beams_enabled(1) = 2

                    usemask = .true.
                    maskb1(5) = .true.

                    partons_enabled = quarkChannel

                    if (ips < 45+OFFSET) then
                        maskb2(ips-6-33-OFFSET) = .true.
                        if (ips-6-33 == OFFSET) then
                            partons_enabled = gluonChannel
                        else
                            partons_enabled = quarkChannel
                        endif
                        taups = 1
                    elseif (ips < 56+OFFSET) then
                        maskb2(ips-6-33-11-OFFSET) = .true.
                        if (ips-6-33-11 == OFFSET) then
                            partons_enabled = gluonChannel
                        else
                            partons_enabled = quarkChannel
                        endif
                        taups = 2
                    else
                        maskb2(ips-6-33-22-OFFSET) = .true.
                        if (ips-6-33-22 == OFFSET) then
                            partons_enabled = gluonChannel
                        else
                            partons_enabled = quarkChannel
                        endif
                        taups = 3
                    endif
                    return
#else
                if (any(ips == [1+OFFSET,2+OFFSET,3+OFFSET])) then
                    ! light line, qb
                    currentContrib = 1
                    maxbeams = 1
                    beams_enabled(1) = 1

                    usemask = .true.
                    maskb1(-5:-1) = .true.
                    maskb1(1:5) = .true.

                    maskb2(5) = .true.

                    partons_enabled = quarkChannel

                    taups = ips - OFFSET
                    return

                elseif (any(ips == [4+OFFSET,5+OFFSET,6+OFFSET])) then
                    ! light line gb
                    currentContrib = 1
                    maxbeams = 1
                    beams_enabled(1) = 1

                    usemask = .true.
                    maskb1(0) = .true.
                    maskb2(5) = .true.

                    partons_enabled = gluonChannel 

                    taups = ips - (3+OFFSET)
                    return

                elseif (any(ips == [7+OFFSET,8+OFFSET,9+OFFSET])) then
                    ! light line bq
                    currentContrib = 1
                    maxbeams = 1
                    beams_enabled(1) = 2

                    usemask = .true.
                    maskb1(5) = .true.

                    maskb2(-5:-1) = .true.
                    maskb2(1:5) = .true.

                    partons_enabled = quarkChannel

                    taups = ips - (6+OFFSET)
                    return

                elseif (any(ips == [10+OFFSET,11+OFFSET,12+OFFSET])) then
                    ! light line bg
                    currentContrib = 1
                    maxbeams = 1
                    beams_enabled(1) = 2

                    usemask = .true.
                    maskb1(5) = .true.
                    maskb2(0) = .true.

                    partons_enabled = gluonChannel 

                    taups = ips - (9+OFFSET)
                    return

#endif
                elseif (ips == 1) then
                    currentContrib = 2
                    maxbeams = 1
                    beams_enabled(1) = 2
                    usemask = .true.
                    maskb1([-1,-3,2,4]) = .true.
                    maskb2(-5:-1) = .true.
                    maskb2(1:4) = .true.
                    partons_enabled = quarkChannel
                    taups = 0
                    return
                elseif (ips == 2) then
                    currentContrib = 2
                    maxbeams = 1
                    beams_enabled(1) = 2
                    usemask = .true.
                    maskb1([-1,-3,2,4]) = .true.
                    maskb2(5) = .true.
                    partons_enabled = quarkChannel
                    taups = 0
                    return
                elseif (ips == 3) then
                    currentContrib = 2
                    maxbeams = 1
                    beams_enabled(1) = 1
                    usemask = .true.
                    maskb1(-5:-1) = .true.
                    maskb1(1:4) = .true.
                    maskb2([-1,-3,2,4]) = .true.
                    partons_enabled = quarkChannel
                    taups = 0
                    return
                elseif (ips == 4) then
                    currentContrib = 2
                    maxbeams = 1
                    beams_enabled(1) = 1
                    usemask = .true.
                    maskb1(5) = .true.
                    maskb2([-1,-3,2,4]) = .true.
                    partons_enabled = quarkChannel
                    taups = 0
                    return
                elseif (ips == 5) then
                    currentContrib = 2
                    maxbeams = 1
                    beams_enabled(1) = 2
                    usemask = .true.
                    maskb1([-1,-3,2,4]) = .true.
                    maskb2(0) = .true.
                    partons_enabled = gluonChannel
                    taups = 0
                    return
                elseif (ips == 6) then
                    currentContrib = 2
                    maxbeams = 1
                    beams_enabled(1) = 1
                    usemask = .true.
                    maskb1(0) = .true.
                    maskb2([-1,-3,2,4]) = .true.
                    partons_enabled = gluonChannel
                    taups = 0
                    return
                elseif (ips == 7) then
                    currentContrib = 3
                    maxbeams = 2
                    beams_enabled(1:2) = [1,2]
                    usemask = .true.
                    maskb1([-3,-1,2,4,5]) = .true.
                    maskb2([-3,-1,2,4,5]) = .true.
                    partons_enabled = quarkChannel + gluonChannel
                    taups = 0
                    return
                endif
            elseif (kpart == kreal) then
                if (ips == 1) then
                    ! 165 nlo, real, quarks, just beam 1
                    currentContrib = 1
                    maxbeams = 1
                    beams_enabled(1) = 1
                    usemask = .true.
                    maskb1(-5:-1) = .true.
                    maskb1(1:5) = .true.
                    maskb2(5) = .true.
                    partons_enabled = quarkChannel
                    taups = 0
                    return
                elseif (ips == 2) then
                    ! 165 nlo, real, quarks, just beam 2
                    currentContrib = 1
                    maxbeams = 1
                    beams_enabled(1) = 2
                    usemask = .true.
                    maskb1(5) = .true.
                    maskb2(-5:-1) = .true.
                    maskb2(1:5) = .true.
                    partons_enabled = quarkChannel
                    taups = 0
                    return
                elseif (ips == 3) then
                    ! 165 nlo, real, gluon, just beam 1
                    currentContrib = 1
                    maxbeams = 1
                    beams_enabled(1) = 1
                    usemask = .true.
                    maskb1(0) = .true.
                    maskb2(5) = .true.
                    partons_enabled = gluonChannel
                    taups = 0
                    return
                elseif (ips == 4) then
                    ! 165 nlo, real, gluon, just beam 2
                    currentContrib = 1
                    maxbeams = 1
                    beams_enabled(1) = 2
                    usemask = .true.
                    maskb1(5) = .true.
                    maskb2(0) = .true.
                    partons_enabled = gluonChannel
                    taups = 0
                    return
                elseif (ips == 5) then
                    currentContrib = 2
                    maxbeams = 1
                    beams_enabled(1) = 2
                    usemask = .true.
                    maskb1([-1,-3,2,4]) = .true.
                    maskb2(-5:-1) = .true.
                    maskb2(1:5) = .true.
                    partons_enabled = quarkChannel
                    taups = 0
                    return
                elseif (ips == 6) then
                    currentContrib = 2
                    maxbeams = 1
                    beams_enabled(1) = 1
                    usemask = .true.
                    maskb1(-5:-1) = .true.
                    maskb1(1:5) = .true.
                    maskb2([-1,-3,2,4]) = .true.
                    partons_enabled = quarkChannel
                    taups = 0
                    return
                elseif (ips == 7) then
                    currentContrib = 2
                    maxbeams = 1
                    beams_enabled(1) = 2
                    usemask = .true.
                    maskb1([-1,-3,2,4]) = .true.
                    maskb2(0) = .true.
                    partons_enabled = gluonChannel
                    taups = 0
                    return
                elseif (ips == 8) then
                    currentContrib = 2
                    maxbeams = 1
                    beams_enabled(1) = 1
                    usemask = .true.
                    maskb1(0) = .true.
                    maskb2([-1,-3,2,4]) = .true.
                    partons_enabled = gluonChannel
                    taups = 0
                    return
                elseif (ips == 9) then
                    currentContrib = 3
                    maxbeams = 2
                    beams_enabled(:) = [1,2]
                    usemask = .false.
                    maskb1(:) = .true.
                    maskb2(:) = .true.
                    partons_enabled = quarkChannel + gluonChannel
                    taups = 0
                    return
                endif
            endif


        elseif (nproc == 1610) then
            ! it even makes sense to split up the lo into at least beams
            if (kpart == klord) then
                currentContrib = 1
                maxbeams = 2
                beams_enabled(:) = [1,2]
                usemask = .false.
                maskb1(:) = .true.
                maskb2(:) = .true.
                partons_enabled = gluonChannel + quarkChannel
                taups = 0
                return
            elseif (kpart == kreal) then
                if (ips == 1) then
                     currentContrib = 1
                     maxbeams = 2
                     beams_enabled(:) = [1,2]
                     usemask = .false.
                     maskb1(:) = .true.
                     maskb2(:) = .true.
                     partons_enabled = gluonChannel + quarkChannel
                     taups = 0
                     return
                elseif (ips == 2) then
                    currentContrib = 2
                    maxbeams = 2
                    beams_enabled(:) = [1,2]
                    usemask = .false.
                    maskb1(:) = .true.
                    maskb2(:) = .true.
                    partons_enabled = gluonChannel + quarkChannel
                    taups = 0
                    return
                elseif (ips == 3) then
                    currentContrib = 3
                    maxbeams = 2
                    beams_enabled(:) = [1,2]
                    usemask = .true.
                    maskb1([-3,-1,2,4,5]) = .true.
                    maskb2([-3,-1,2,4,5]) = .true.
                    partons_enabled = gluonChannel + quarkChannel
                    taups = 0
                    return
                elseif (ips == 4) then
                    ! lxh interference, real on light, real on heavy
                    currentContrib = 4
                    maxbeams = 2
                    beams_enabled(:) = [1,2]
                    usemask = .false.
                    maskb1(:) = .true.
                    maskb2(:) = .true.
                    partons_enabled = gluonChannel + quarkChannel
                    taups = 0
                    return
                elseif (ips == 5) then
                    ! lxh interference, real on light, virt on heavy
                    currentContrib = 4
                    maxbeams = 2
                    beams_enabled(:) = [1,2]
                    usemask = .false.
                    maskb1(:) = .true.
                    maskb2(:) = .true.
                    partons_enabled = gluonChannel + quarkChannel
                    taups = 0
                    return
                elseif (any(ips == [6,7])) then
                    currentContrib = 5
                    maxbeams = 2
                    beams_enabled(:) = [1,2]
                    usemask = .false.
                    maskb1(:) = .true.
                    maskb2(:) = .true.
                    partons_enabled = gluonChannel + quarkChannel
                    taups = 0
                    return
                elseif (any(ips == [8,9])) then
                    currentContrib = 6
                    maxbeams = 2
                    beams_enabled(:) = [1,2]
                    usemask = .false.
                    maskb1(:) = .true.
                    maskb2(:) = .true.
                    partons_enabled = gluonChannel + quarkChannel
                    taups = 0
                    return
                endif
            elseif (kpart == kvirt) then
                if (ips == 1) then
                     currentContrib = 1
                     maxbeams = 2
                     beams_enabled(:) = [1,2]
                     usemask = .false.
                     maskb1(:) = .true.
                     maskb2(:) = .true.
                     partons_enabled = gluonChannel + quarkChannel
                     taups = 0
                     coeffonly = .true.
                     return
                elseif (ips == 2) then
                    currentContrib = 2
                    maxbeams = 2
                    beams_enabled(:) = [1,2]
                    usemask = .false.
                    maskb1(:) = .true.
                    maskb2(:) = .true.
                    partons_enabled = gluonChannel + quarkChannel
                    taups = 0
                    coeffonly = .true.
                    return
                elseif (ips == 3) then
                    currentContrib = 3
                    maxbeams = 2
                    beams_enabled(:) = [1,2]
                    usemask = .true.
                    maskb1([-3,-1,2,4,5]) = .true.
                    maskb2([-3,-1,2,4,5]) = .true.
                    partons_enabled = gluonChannel + quarkChannel
                    taups = 0
                    coeffonly = .true.
                    return
                elseif (ips == 4) then
                    currentContrib = 4
                    maxbeams = 2
                    beams_enabled(:) = [1,2]
                    usemask = .false.
                    maskb1(:) = .true.
                    maskb2(:) = .true.
                    partons_enabled = gluonChannel + quarkChannel
                    taups = 0
                    coeffonly = .true.
                    return
                elseif (ips == 5) then
                    currentContrib = 4
                    maxbeams = 2
                    beams_enabled(:) = [1,2]
                    usemask = .false.
                    maskb1(:) = .true.
                    maskb2(:) = .true.
                    partons_enabled = gluonChannel + quarkChannel
                    taups = 0
                    coeffonly = .true.
                    return
                elseif (any(ips == [6,7])) then
                    currentContrib = 5
                    maxbeams = 2
                    beams_enabled(:) = [1,2]
                    usemask = .false.
                    maskb1(:) = .true.
                    maskb2(:) = .true.
                    partons_enabled = gluonChannel + quarkChannel
                    taups = 0
                    coeffonly = .true.
                    return
                elseif (any(ips == [8,9])) then
                    currentContrib = 6
                    maxbeams = 2
                    beams_enabled(:) = [1,2]
                    usemask = .false.
                    maskb1(:) = .true.
                    maskb2(:) = .true.
                    partons_enabled = gluonChannel + quarkChannel
                    taups = 0
                    coeffonly = .true.
                    return
                endif
            elseif (origKpart == ksnlo) then
                if (ips == 1) then
                     currentContrib = 1
                     maxbeams = 1
                     beams_enabled(1) = 1
                     usemask = .true.
                     maskb1(:) = .true.
                     maskb2(5) = .true.
                     partons_enabled = gluonChannel + quarkChannel
                     taups = 0

                     coeffonly = .true.
                     return
                 elseif (ips == 2) then
                     currentContrib = 1
                     maxbeams = 1
                     beams_enabled(1) = 2
                     usemask = .true.
                     maskb1(5) = .true.
                     maskb2(:) = .true.
                     partons_enabled = gluonChannel + quarkChannel
                     taups = 0

                     coeffonly = .true.
                     return
                elseif (ips == 3) then
                    currentContrib = 2
                    maxbeams = 1
                    beams_enabled(1) = 2
                    usemask = .false.
                    maskb1(:) = .true.
                    maskb2(:) = .true.
                    partons_enabled = gluonChannel + quarkChannel

                    taups = 0
                    coeffonly = .true.
                    return
                elseif (ips == 4) then
                    currentContrib = 2
                    maxbeams = 1
                    beams_enabled(1) = 1
                    usemask = .false.
                    maskb1(:) = .true.
                    maskb2(:) = .true.
                    partons_enabled = gluonChannel + quarkChannel

                    taups = 0
                    coeffonly = .true.
                    return
                elseif (ips == 5) then
                    currentContrib = 3
                    maxbeams = 2
                    beams_enabled(1:2) = [1,2]
                    usemask = .true.
                    maskb1([-3,-1,2,4,5]) = .true.
                    maskb2([-3,-1,2,4,5]) = .true.
                    partons_enabled = gluonChannel + quarkChannel

                    taups = 0
                    coeffonly = .true.
                    return
                endif
            elseif (origKpart == knnlo) then
                if (ips == 1) then
                     currentContrib = 1
                     maxbeams = 1
                     beams_enabled(1) = 1
                     usemask = .true.
                     maskb1(:) = .true.
                     maskb2(5) = .true.
                     partons_enabled = gluonChannel + quarkChannel
                     taups = 0

                     coeffonly = .true.
                     return
                 elseif (ips == 2) then
                     currentContrib = 1
                     maxbeams = 1
                     beams_enabled(1) = 2
                     usemask = .true.
                     maskb1(5) = .true.
                     maskb2(:) = .true.
                     partons_enabled = gluonChannel + quarkChannel
                     taups = 0

                     coeffonly = .true.
                     return
                elseif (ips == 3) then
                     currentContrib = 2
                     maxbeams = 1
                     beams_enabled(1) = 1
                     usemask = .true.
                     maskb1(:) = .true.
                     maskb2([-3,-1,2,4]) = .true.
                     partons_enabled = gluonChannel + quarkChannel
                     taups = 0

                     coeffonly = .true.
                     return
                elseif (ips == 4) then
                     currentContrib = 2
                     maxbeams = 1
                     beams_enabled(1) = 2
                     usemask = .true.
                     maskb1([-3,-1,2,4]) = .true.
                     maskb2(:) = .true.
                     partons_enabled = gluonChannel + quarkChannel
                     taups = 0

                     coeffonly = .true.
                     return
                elseif (ips == 5) then
                     currentContrib = 3
                     maxbeams = 2
                     beams_enabled(1:2) = [1,2]
                     usemask = .true.
                     maskb1([-3,-1,2,4,5]) = .true.
                     maskb2([-3,-1,2,4,5]) = .true.
                     partons_enabled = gluonChannel + quarkChannel
                     taups = 0

                     coeffonly = .true.
                     return

                endif
            endif
        endif

        write (*,*) "IPS = ", ips
        write (*,*) "kpart = ", kpart
        write (*,*) "origKpart = ", origKpart
        write (*,*) "nproc = ", nproc
        error stop "please setup part in setup_singletop"

    end subroutine

end module
