!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      !----------------------------------------------------------------------
      !      This replaces the original includedipole
      !      It calls the original and then the user one
      logical function includedipole(nd,ptrans)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      real(dp) ptrans(mxpart,4)
      integer nd
      logical :: mcfm_includedipole

      ! it looks like all (nd  /=  0) automatically have their momenta
      ! stored in the ptilde common; do this also for nd=0
      if (nd  ==  0) call storeptilde(nd,ptrans)

      ! first call the original MCFM includedipole
      includedipole = mcfm_includedipole(nd,ptrans)

      end

      logical function mcfm_includedipole(nd,ptrans)
     &                 result(mcfmincdipole)
       use types
       use m_gencuts
       use qtResummation_params, only: qtcutoff
       use nnlo_z1jet, only: passed_taucut_1jet
       use SCET, only: scetreweight, useQT
       use ptveto, only: jetptveto, usept
c--- This function returns TRUE if the specified point ptrans,
c--- corresponding to dipole nd (nd=0 => real radiation),
c--- should be included
      implicit none
      include 'mxpart.f'
      include 'constants.f'
      include 'clustering.f'
      include 'npart.f'
      include 'ptilde.f'
      include 'jetlabel.f'
      include 'kprocess.f'
      include 'frag.f'
      include 'phot_dip.f'
      include 'nqcdjets.f'
      include 'nproc.f'
      include 'notag.f'
      include 'taucut.f'
      include 'hdecaymode.f'
      include 'ewcorr.f'
      include 'energy.f'
      include 'runstring.f'
      include 'ipsgen.f'
      include 'kpart.f'
      include 'Rcut.f'
      include 'makecuts.f'
      include 'limits.f'
!      include 'hybridiso.f'

      real(dp) ptrans(mxpart,4),pjet(mxpart,4),pttwo,ptthree,qt,pt
      real(dp) :: ptfour,ptsix
      integer j,nd,isub,nfail
      logical failedgencuts,photoncuts,filterWbbmas,
     &     photonfailed,filterW_bjet,is_photon,photoncuts_ew,
     &     passedcuts_wgamma_ew,passedcuts_wgammajet_ew,passedcuts_w_ew
      integer count_photo,nphotons
      logical:: iso, passed_iso, passed_taucut

      real(dp) :: pt34

c--- default: include this contribution
      mcfmincdipole=.true.

c--- isub=1 for dipole subtractions, isub=0 for real radiation
      if (nd  >  0) then
        isub=1
      else
        isub=0
      endif

c Special cuts for processes with QED corrections
      if ((kcase == kWln_ew) .or. (kcase == kWln_aq)) then
c Perform generic EW Wgamma cuts only
        if (passedcuts_w_ew(isub,ptrans) .eqv. .false.) then
          mcfmincdipole=.false.
          return
        else
          pjet(:,:)=ptrans(:,:)
          jets=0
          ptildejet(1:npart+2,1:4,nd)=pjet(1:npart+2,1:4)
        endif
        return
      endif
      if (kcase == kWga_ew) then
c Perform generic EW Wgamma cuts only
        if (passedcuts_wgamma_ew(isub,ptrans) .eqv. .false.) then
          mcfmincdipole=.false.
          return
        else
          pjet(:,:)=ptrans(:,:)
          jets=0
          ptildejet(1:npart+2,1:4,nd)=pjet(1:npart+2,1:4)
        endif
        return
      endif
      if (kcase == kWgajew) then
c Perform generic EW Wgamma cuts (including a jet) only
        if (passedcuts_wgammajet_ew(isub,ptrans) .eqv. .false.) then
          mcfmincdipole=.false.
          return
        else
          pjet(:,:)=ptrans(:,:)
          jets=1
          ptildejet(1:npart+2,1:4,nd)=pjet(1:npart+2,1:4)
        endif
        return
      endif
c end special cuts

      nphotons=count_photo()
      if (nphotons  >  0) then
c--- Photons: Frixione isolation cuts if no fragmentation included
         if (frag .eqv. .false.) then
            nfail=0
            do j=3,mxpart
               if(is_photon(j)) then
!                  if (hybridiso) then
!                    call do_hybridiso(ptrans,passed_iso,j,isub)
!                  else
                    call frix(ptrans,passed_iso,j,isub)
!                  endif
                  if(passed_iso .eqv. .false.) then
                    if (kcase == kWgajew) then
                      nfail=nfail+1
                    else
                      mcfmincdipole=.false.
                      return
                    endif
                  endif
               endif
            enddo
            if (kcase == kWgajew) then
              if (((nd == 0)  .and. (nfail == 2))
     &        .or. ((nd > 0) .and.( nfail == 1))) then
                mcfmincdipole=.false.
                return
              endif
            endif
            call genclustphotons(ptrans,Rcut,pjet,isub)
         else
c--- Photons: not Frixione, need fragmentation and isolation
c---- do not want to cluster partons inside of jet cone, allow
c---- isolation to describe these regions, therefore use
c---- genclustphotons here
            call genclustphotons(ptrans,Rcut,pjet,isub)
c            call genclust2(ptrans,Rcut,pjet,isub)
c---  Isolate photon
            do j=3,mxpart
            if (is_photon(j)) then
               if (iso(ptrans,j,isub,nd) .eqv. .false.)then
                  mcfmincdipole=.false.
                  return
               endif
            endif
            enddo
         endif
c--- check the photon cuts
         if ((kcase == kWga_ew) .or. (kcase == kWgajew)) then
           photonfailed=photoncuts_ew(isub,pjet)
         else
           photonfailed=photoncuts(pjet)
         endif
         if (photonfailed) then
            mcfmincdipole=.false.
            return
         endif
      else
c--- No photons: the usual case
         call genclust2(ptrans,Rcut,pjet,isub,nd)
      endif


c--- perform mass cuts
      call masscuts(pjet,*999)

c--- fill ptilde array as persistent storage for the jet momenta
      ! GPS written in compact f90 form - in case we wish to copy it elsewhere
      ! (e.g. earlier, so that jets are defined even when includedipole is false)
      ptildejet(1:npart+2,1:4,nd)=pjet(1:npart+2,1:4)

c--- for the Wbb process, we divide up contributions in a specific way;
c--- therefore filter events using special code and skip normal jet testing
c--- NOTE: only for process numbers > 400 (20 and 25 should be handled normally)
      if ((kcase==kWbbmas) .and. (nproc  >  400)) then
        mcfmincdipole=filterWbbmas()
      if (mcfmincdipole .eqv. .false.) return
        if (makecuts) then
          failedgencuts=gencuts(pjet,jets)
          if (failedgencuts) mcfmincdipole=.false.
        endif
      goto 99
      endif


c--- for the Wb+X process, we divide up contributions in a specific way;
c--- therefore filter events using special code and skip normal jet testing
      if (kcase==kW_bjet) then
        mcfmincdipole=filterW_bjet()
      if (mcfmincdipole .eqv. .false.) return
        if (makecuts) then
          failedgencuts=gencuts(pjet,jets)
          if (failedgencuts) mcfmincdipole=.false.
        endif
      goto 99
      endif

      ! Z with qT subtractions.
      ! could probably be removed now and merged with general case below
      if (usescet .and. useqt .and. any(nprocbelow == [31,32])) then
          call makeqtcut(ptrans,pjet,isub,passed_taucut,nd)
          mcfmincdipole = passed_taucut
          if (mcfmincdipole .eqv. .false.) return
      ! Z+jet
      elseif (usescet .and.  any(nprocbelow == [11,16,41,204,210])) then
          ! qT cut for for qT resummation or N3LO
          ! WARNING: for this qT cutoff we must either use tau minimization or 
          ! for the hardest jet axis with ptjetmin=0.0, etajetmax=99.9 set in
          ! the input file to ensure that maketaucut has access to a hardest jet
          pt34 = pttwo(3,4,ptrans)
          !if (pt34 < 1._dp .or. (pt34 < pt34min .or. pt34 > pt34max)) then
          !if (jets < 1) then ! when interested in Z+jet process

          ! with a pt cut we discard all events outside the range right away
          if (pt34min > 0._dp) then
            if (pt34 < pt34min .or. pt34 > pt34max) then
                passed_taucut = .false.
                mcfmincdipole = .false.
                return
            endif
          else
            ! otherwise we're in V+jet mode with a 1-jet requirement
            if (jets < 1) then
              passed_taucut = .false.
              mcfmincdipole = .false.
              return
            endif
          endif

          if (abovecut) then

              ! minimization
              !passed_taucut = passed_taucut_1jet(ptrans,scetreweight,taucut)
              !mcfmincdipole = passed_taucut

              ! hardest jet axis, seems to have better power corrections
              call maketaucut(ptrans,pjet,jets,isub,passed_taucut,nd)
              mcfmincdipole = passed_taucut

              if (mcfmincdipole .eqv. .false.) return
          else
              passed_taucut = .true.
              mcfmincdipole = .true.
          endif
      !elseif (nproc == 44) then
      !    !!! debug 
      !    !!! Z+2jet part with taucut
      !    if (jets < 1) then ! when interested in Z+jet process
      !        passed_taucut = .false.
      !        mcfmincdipole = .false.
      !        return
      !    else
      !        passed_taucut = passed_taucut_1jet(ptrans,scetreweight,taucut)
      !        mcfmincdipole = passed_taucut
      !        if (mcfmincdipole .eqv. .false.) return
      !    endif
      elseif (usescet .and. (nproc /= 1610) .and. (nproc /= 1650)) then
        if (useQT .or. useqt_nnlo) then
          call makeQTcut(ptrans,pjet,isub,passed_taucut,nd)
        elseif (usept) then
          call makeptcut(ptrans,pjet,isub,passed_taucut,nd)
        else
          call maketaucut(ptrans,pjet,jets,isub,passed_taucut,nd)
        endif

        mcfmincdipole=passed_taucut
c           write(6,*) 'includedipole: nd,mcfmincdipole',nd,mcfmincdipole
c           if (passed_taucut .eqv. .false.) write(6,*) 'tau failed: ',nd
        if (mcfmincdipole .eqv. .false.) return
        if ((inclusive .eqv. .false.) .and. (jets  /=  ntau)) then
          mcfmincdipole = .false.
          return
        endif
      elseif (origkpart == kresummed) then
          ! cut on boson transverse momentum
          ! but fully inclusive on any real radiation
          if (abovecut) then
              if (any(nprocbelow == [1,6,31,32,33,111,112,119,285,2851])) then
                  qt=pttwo(3,4,ptrans)
              elseif (any(nprocbelow == [120,290,295,300,305])) then
                  qt=ptthree(3,4,5,ptrans)
              elseif (any(nprocbelow == [61,71,711,72,76,761,77,81,82,900,91,92,93,96,97,98,101,104,110])) then
                  qt=ptfour(3,4,5,6,ptrans)
              elseif (any(nprocbelow == [106])) then
                  qt=ptsix(3,4,5,6,7,8,ptrans)
              else
                  error stop "implement small qT cut in includedipole.f"
              endif

              if ((qt < qtcutoff)) then
                  mcfmincdipole = .false.
                  return
              endif
          endif
      else
        ! for V+jet calculation with minimum V transverse momentum (for matching with N4LL)
        if (pt34min > 0._dp) then
            pt34 = pttwo(3,4,ptrans)
            if ((pt34 < qtcutoff) .or. (pt34 < pt34min) .or. (pt34 > pt34max)) then
                mcfmincdipole = .false.
                return
            endif
        else
            !--- for a normal calculation,      
            !--- if the number of jets is not correct, then do not include dipole
            if ((clustering .and. (jets  /=  nqcdjets-notag)
     &             .and. (inclusive .eqv. .false.)) .or.
     &          (clustering .and. (jets  <  nqcdjets-notag)
     &             .and. (inclusive .eqv. .true.))) then
                mcfmincdipole=.false.
                return
            endif
        endif
      endif

      if (jetptveto < 1.e4_dp) then
! cut on maximum jet pt (jet pt veto)
       if (any(nprocbelow == [1,6,31,32,33,111,112,119,285,2851])) then
           qt=max(pt(5,pjet),pt(6,pjet),0._dp)
       elseif (any(nprocbelow == [120,290,295,300,305])) then
           qt=max(pt(6,pjet),pt(7,pjet),0._dp)
       elseif (any(nprocbelow == [61,71,711,72,76,761,77,81,82,900,91,92,93,96,97,98,101,104,110,8211])) then
           qt=max(pt(7,pjet),pt(8,pjet),0._dp)
       elseif (any(nprocbelow == [106])) then
           qt=max(pt(9,pjet),pt(10,pjet),0._dp)
       else
           error stop "implement pt(veto) cut in includedipole.f"
       endif

       if (qt > jetptveto) then
           mcfmincdipole = .false.
           return
       endif
      endif


c--- check the lepton cuts, if necessary
      if (makecuts) then
        failedgencuts=gencuts(pjet,jets)
        if (failedgencuts) then
          mcfmincdipole=.false.
          return
        endif
      endif

   99 continue

      return

  999 continue
      mcfmincdipole=.false.
      return

      end

