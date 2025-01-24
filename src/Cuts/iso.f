!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

c--- This is a driver function for photon isolation, the relevant photon isolation
c--- critera are selected and applied

      function iso(p,phot_id,isub,nd)
       implicit none
      include 'types.f'
      logical:: iso

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'frag.f'
      include 'useet.f'
      include 'phot_dip.f'
      include 'z_dip.f'
      include 'lastphot.f'
      include 'mpicommon.f'
      include 'first.f'
      real(dp):: p(mxpart,4)
      integer:: isub,nd,imode
      integer:: phot_id ! refers to which photon we are isolating
      logical:: photo_iso_phys,photo_iso_z!,photo_iso,iso_old

c--- imode: set the behaviour of the isolation requirement (default: imode=0)
c---   imode=0 : scaling cut for epsilon_h < 1, fixed cut for epsilon_h > 1
c---   imode=1 : scaling cut for all epsilon_h
c---   imode=2 : fixed cut for all epsilon_h
      imode=0

      iso=.true.

c Check if no isolation required
      if ((abs(cone_ang) < 1.e-4_dp).or.(abs(epsilon_h) < 1.e-4_dp)) then
        if (first) then
!$omp master
        if (rank == 0) then
        write(6,*)'****************************************************'
        write(6,*)'*                                                  *'
        write(6,*)'*         No photon isolation cuts applied         *'
        write(6,*)'*                                                  *'
        write(6,*)'****************************************************'
        endif
!$omp end master
        first=.false.
      endif
        return
      endif

c--- for imode=0, decide which cut to use depending on epsilon_h
c--   epsilon_h < 1 : epsilon_h corresponds to a pt fraction i.e. pt(had) < epsilon_h pt_gamma
c--   epsilon_h > 1 : treat it as an E_t max in cone   i.e. pt(had) < epsilon_h
       if (imode == 0) then
         if (epsilon_h < 0.9999_dp) then
         imode=1
         else
         imode=2
         endif
       endif

c===== NEW ISOLATION : WORK OUT WHICH ROUTINE TO CALL

c========== are we doing fragmentation and is this the
c========== photon produced by fragmentation?
       if((fragint_mode) .and. (phot_id == lastphot)) then
          iso=photo_iso_z(p,phot_id,z_frag,imode,isub)
c========== are we doing a photon dipole and is this the
c========== photon produced by fragmentation?
       elseif(phot_dip(nd) .and. (phot_id == lastphot)) then
          iso=photo_iso_z(p,phot_id,z_dip(nd),imode,isub)
c========== we are in default mode, BORN, virt, real, non-photon dipole
       else
          iso=photo_iso_phys(p,phot_id,imode,isub)
       endif

c===== OLD ISOLATION
c       iso_old=photo_iso(p,isub,phot_dip(nd),phot_id,nd,imode)
c       if (iso .neqv. iso_old) then
c         write(6,*)'WARNING:',iso,'(new) vs',iso_old,' (old) for nd=',nd
c       else
c         write(6,*) 'OKAY'
c       endif

c--- write out isolation parameters
       if    (first .and. (imode == 1)) then
!$omp master
        if (rank == 0) then
        write(6,*)'************** Photons Isolated     ****************'
        write(6,*)'*                                                  *'
        write(6,99)'*    E_t(had) in cone',cone_ang,' < ',epsilon_h,
     &   ' E_t(phot)     *'
        write(6,*)'*                                                  *'
        write(6,*)'****************************************************'
        endif
!$omp end master
        first=.false.
      elseif (first .and. (imode == 2)) then
!$omp master
        if (rank == 0) then
        write(6,*)'************** Photons Isolated     ****************'
        write(6,*)'*                                                  *'
        write(6,96)'* E_t (had) in cone',cone_ang,' < ',epsilon_h,
     &   'GeV    *'
        write(6,*)'*                                                  *'
        write(6,*)'****************************************************'
        endif
!$omp end master
        first=.false.
      endif

c      endif

      return

 99   format(1x,a21,f6.2,a3,f6.2,a16)
 96   format(1x,a19,f6.2,a4,f6.2,a17)
      end
