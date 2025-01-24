!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module Integration
    use types
    implicit none

    enum, bind(c)
          enumerator :: nloReal=1, nloVirt=2, &
               nnloWilsonReal = 3, nnloWilsonVirt = 4, &
               nnloBelow=5, nnloVirtAbove=6, nnloRealAbove=7, &
               snloBelow = 8, snloAbove = 9, lord=10, nloFrag=11, &
               nloRealExtra = 12, nloResummed=13, nloResAbove = 14, &
               nnloResVirtAbove=15, nnloResRealAbove=16, nloResexp = 17, &
               qtnloBelow = 18, qtnnloBelow = 19, qtn3loBelow = 20, &
               nloResVetoReal = 21, nloResVetoVirt = 22, &
               nnloResVetoBelow = 23, nnloResVetoRealAbove = 24, nnloResVetoVirtAbove = 25
    endenum

    ! adjust maxParts to be the maximum number used above
    integer, parameter :: maxParts = 25

#define SINGLETOP_NNLO_REALCHANNELS 12
#if SINGLETOP_NNLO_REALCHANNELS == 66
    integer, parameter :: maxIPS = 73
#else
    integer, parameter :: maxIPS = 19
#endif

    integer, parameter :: ndim_incr(maxParts) = [3,1,3,1,2,1,3,2,0,0,1,3,2,0,1,3,2,0,0,0,3,1,2,3,1]

    character(*), parameter :: partStrings(maxParts) = [ &
        "NLO real           ", &
        "NLO virt           ", &
        "NNLO Wilson real   ", &
        "NNLO Wilson virt   ", &
        "NNLO below cut     ", &
        "NNLO virt above cut", &
        "NNLO real above cut", &
        "SNLO below cut     ", &
        "SNLO above cut     ", &
        "LO                 ", &
        "NLO frag           ", &
        "NLO real extra     ", &
        "Resummed           ", &
        "NLO Resummed above ", &
        "NNLO Res. virtabove", &
        "NNLO Res. realabove", &
        "Res. f.o. expansion", &
        "qT NLO below       ", &
        "qT NNLO below      ", &
        "qT N3LO below      ", &
        "NLO Res. veto real ", &
        "NLO Res. veto virt ", &
        "NNLO Res.veto below", &
        "NNLO Res.veto r.abv", &
        "NNLO Res.veto v.abv"]

    real(dp), save :: warmupChisqGoal

end module
