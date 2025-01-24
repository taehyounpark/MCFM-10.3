!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine setupscet(nprocabove)
      implicit none
c---- Routine to setup variables for performing SCET calculation;
c---- inspects value of nproc (passed via common block) and
c---- returns nprocabove, the corresponding process with one
c---- additional jet;  stops with error message if appropriate.
      include 'types.f'
      include 'nproc.f'
      include 'taucut.f'
      integer:: nprocabove

      if     ((nproc == 1) .or. (nproc == 6)
     &   .or. (nproc == 31) .or. (nproc == 32) .or. (nproc == 33)) then
        nprocabove=nproc+10
        ntau=0
      elseif (nproc == 61) then
        nprocabove=461
        ntau=0
      elseif (nproc == 71) then
        nprocabove=471
        ntau=0
      elseif (nproc == 76) then
        nprocabove=476
        ntau=0
      elseif (nproc == 81) then
        nprocabove=481
        ntau=0
      elseif (nproc == 82) then
        nprocabove=482
        ntau=0
      ! W+H(tau tau) and W-H(tau tau)
      elseif ((nproc == 91) .or. (nproc == 96)) then
        nprocabove=nproc+519
        ntau=0
      elseif (nproc == 900) then
        nprocabove=609
        ntau=0
      elseif ((nproc == 92) .or. (nproc == 97)) then
        nprocabove=nproc+519
        ntau=0
      elseif ((nproc == 93) .or. (nproc == 98)) then
        nprocabove=nproc+519
        ntau=0
      elseif ((nproc == 94) .or. (nproc == 99)) then
        nprocabove=nproc+519
        ntau=0
      elseif (nproc == 101) then
        nprocabove=621
        ntau=0
      elseif (nproc == 104) then
        nprocabove=622
        ntau=0
      elseif (nproc == 106) then
        nprocabove=623
        ntau=0
      elseif (nproc == 110) then
        nprocabove=620
        ntau=0
      elseif (nproc == 111) then
        nprocabove=203
        ntau=0
      elseif (nproc == 112) then
        nprocabove=204
        ntau=0
      elseif (nproc == 119) then
        nprocabove=210
        ntau=0
      elseif (nproc == 120) then
        nprocabove=205
        ntau=0
      ! q qbar -> gam gam
      elseif ((nproc == 285)) then
        nprocabove=nproc+1
        ntau=0
      ! gg -> gam gam
      elseif (nproc == 2851) then
        nprocabove = 2861
        ntau = 0
      elseif (nproc ==290 ) then
        nprocabove=292
        ntau=0
      elseif (nproc ==295 ) then
        nprocabove=297
        ntau=0
      ! Zgam ( e+ e-)
      elseif (nproc == 300) then
        nprocabove=302
        ntau=0
      ! Zgam (nu nu)
      elseif (nproc == 305) then
        nprocabove=307
        ntau=0
      elseif (nproc == 3000) then
        nprocabove=3002
        ntau=0
      elseif (nproc == 1610) then ! single-top nnlo
        nprocabove = 1650
        ntau=0
      elseif (nproc == 41) then ! Z+jet
        nprocabove = 44
        ntau = 1
! These processes currently untested
!      elseif (nproc == 11) then ! W+ jet
!        nprocabove = 22
!        ntau = 1
!      elseif (nproc == 16) then ! W- jet
!        nprocabove = 27
!        ntau = 1
      elseif (nproc == 204) then
        nprocabove = 272
        ntau = 1
      elseif (nproc == 210) then
        nprocabove = 270
        ntau = 1
      else
        write(6,*) 'This process cannot be computed at NNLO'
        write(6,*) 'or at NLO using SCET or QTCUT.'
        stop
      endif

      return
      end

