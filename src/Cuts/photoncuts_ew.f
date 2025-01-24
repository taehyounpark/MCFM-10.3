!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function photoncuts_ew(isub,pjet)
      implicit none
      include 'types.f'
      logical:: photoncuts_ew
c***********************************************************************
c   Author: J.M. Campbell, 24th January 2011                           *
c       and C. Williams                                                *
c   This routine imposes the photon cuts that are specified in the     *
c   input file. The cuts are applied here (rather than in gencuts.f)   *
c   since they may be necessary for a finite cross section             *
c   and thus should be applied regardless of the value of "makecuts"   *
c                                                                      *
c   Return TRUE if this point FAILS the cuts                           *
c   Modified March 11 to apply mass cuts to m34 for gamgam (CW)        *
c                                                                      *
c***********************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'limits.f'
      include 'kprocess.f'
      include 'leptcuts.f'
      include 'z_dip.f'
      include 'jetlabel.f'
      include 'first.f'
      include 'mpicommon.f'
      include 'runstring.f'
      include 'npart.f'
      logical:: is_lepton,is_photon,is_hadronic
      integer:: i,j,isub,nunobs
      real(dp):: aygam
      integer:: countlept,leptindex(mxpart),countgamm,gammindex(mxpart),
     & countjet,jetindex(mxpart)
      real(dp):: pjet(mxpart,4),pt,etarap,R,pt1,pt2,eta1,eta2
      save countlept,countgamm,countjet,
     & leptindex,gammindex,jetindex
      logical, save :: CMScrack
!$omp threadprivate(countlept,countgamm,countjet)
!$omp threadprivate(leptindex,gammindex,jetindex)
!$omp threadprivate(CMScrack)

      photoncuts_ew=.false.

      if (first) then
      first=.false.
      CMScrack=(index(runstring,'CMScrack') > 0)
c--- write-out the cuts we are using
!$omp master
      if (rank  == 0) then
      write(6,*)
      write(6,*)  '****************** Photon cuts *********************'
      write(6,*)  '*                                                  *'
      if (gammptmax > 0.99e6_dp) then
        write(6,99) '*   pt(photon 1)         >   ',gammptmin,
     &                '                *'
      else
        write(6,97) gammptmin,'    pt(photon)     ',gammptmax,'GeV'
      endif
      write(6,99) '*   pt(photon 2)         >   ',gammpt2,
     &                '                *'
      if((kcase==ktrigam) .or. (kcase==kfourga))then
         write(6,99) '*   pt(photon 3)         >   ',gammpt3,
     &        '                *'
      endif
      if(kcase==kfourga) then
         write(6,99) '*   pt(photon 4)         >   ',gammpt3,
     &        '                *'
      endif
      write(6,97) gammrapmin,'   |eta(photon)|   ',gammrapmax,'   '
      write(6,99) '*   R(photon,lepton)     >   ',Rgalmin,
     &                '                *'
      write(6,99) '*   R(photon,photon)     >   ',Rgagamin,
     &                '                *'
      write(6,99) '*   R(photon,jet)        >   ',Rgajetmin,
     &                '                *'
      if(CMScrack) then
       write(6,*) '*   excluding rapidity window 1.44 < eta < 1.57    *'
      endif
      write(6,*)  '*                                                  *'
      write(6,*)  '****************************************************'
      if((kcase==kgamgam) .or. (kcase==kgg2gam) .or. (kcase==kgmgmjt)) then
       write(6,*)
       write(6,*) '************* M(gam,gam) mass cuts *****************'
       write(6,*) '*                                                  *'
       write(6,98) sqrt(wsqmin),'m34',sqrt(wsqmax)
       write(6,*) '****************************************************'
      endif
      endif
!$omp end master
c--- initialize counters and arrays that will be used to perform cuts
      countlept=0
      countgamm=0
      countjet=0
      do j=3,mxpart
        if (is_lepton(j)) then
          countlept=countlept+1
          leptindex(countlept)=j
        endif
        if (is_photon(j)) then
          countgamm=countgamm+1
          gammindex(countgamm)=j
        endif
        if (is_hadronic(j)) then
          countjet=countjet+1
          jetindex(countjet)=j
        endif
      enddo
      endif

      countgamm=npart-2-isub
      if (kcase == kWgajew) then
        countgamm=2-isub
      endif

c==== option to exclude CMS crack corresponding to 1.44 < eta < 1.57
      if (CMScrack) then
         nunobs=0
         do i=1,countgamm
            aygam=abs(etarap(gammindex(i),pjet))
            if((aygam > 1.44_dp) .and. (aygam < 1.57_dp)) then
               nunobs=nunobs+1
            endif
         enddo
         if (countgamm - nunobs < 1) then
           photoncuts_ew=.true.
           return
         endif
      endif

c     Basic pt and rapidity cuts for photon and lepton-photon separation
      if (countgamm == 1) then
          eta1=abs(etarap(gammindex(1),pjet))
          pt1=pt(gammindex(1),pjet)
          if ( (pt1 < gammptmin) .or. (pt1 > gammptmax) .or.
     &    (eta1 > gammrapmax) .or. (eta1 < gammrapmin)) then
            photoncuts_ew=.true.
            return
          endif
          if ((eta1 > gammvetomin) .and.
     &        (eta1 < gammvetomax)) then
            photoncuts_ew=.true.
            return
          endif
          if (R(pjet,gammindex(1),leptindex(1)) < Rgalmin) then
            photoncuts_ew=.true.
            return
          endif
      endif
      if (countgamm == 2) then
        pt1=pt(gammindex(1),pjet)
        pt2=pt(gammindex(2),pjet)
        eta1=abs(etarap(gammindex(1),pjet))
        eta2=abs(etarap(gammindex(2),pjet))
        nunobs=0
        if ( ( pt1 < gammptmin) .or. (pt1 > gammptmax) .or.
     &       (eta1 > gammrapmax) .or. (eta1 < gammrapmin) .or.
     &       ((eta1 > gammvetomin) .and. (eta1 < gammvetomax)) .or.
     &       (R(pjet,gammindex(1),leptindex(1)) < Rgalmin) ) then
          nunobs=nunobs+1
        endif
        if ( ( pt2 < gammptmin) .or. (pt2 > gammptmax) .or.
     &       (eta2 > gammrapmax) .or. (eta2 < gammrapmin) .or.
     &       ((eta2 > gammvetomin) .and. (eta2 < gammvetomax)) .or.
     &       (R(pjet,gammindex(2),leptindex(1)) < Rgalmin) ) then
          nunobs=nunobs+1
        endif
        if (countgamm - nunobs < 1) then
          photoncuts_ew=.true.
          return
        endif
      endif

cc--- photon-photon separation
c      if (countgamm >= 2) then
c        do j=1,countgamm
c        do k=j+1,countgamm
c          if (R(pjet,gammindex(j),gammindex(k)) < Rgagamin) then
c            photoncuts_ew=.true.
c            return
c          endif
c        enddo
c        enddo
c      endif

cc--- jet-photon separation (if there are 1 or more jets and photons)
c      if ((jets >= 1) .and. (countgamm >= 1)) then
c        do j=1,countgamm
c        do k=1,jets
c          if (R(pjet,gammindex(j),jetindex(k)) < Rgajetmin) then
c            photoncuts_ew=.true.
c            return
c          endif
c        enddo
c        enddo
c      endif


      return




c--- Lines below here are commented out for now;
c---  may be restored at a later date.


c--- jet-photon separation (if there are 1 or more jets and photons)
c      if ((njets > 0) .and. (countgamm > 0)) then
c        do j=1,countgamm
c        do k=1,njets
c          if (R(pjet,gammindex(j),jetindex(k)) < Rjlmin) then
c            gencuts=.true.
c            return
c          endif
c        enddo
c        enddo
c      endif

c--- DEBUG: removed all isolation
cc--- photon/hadron isolation
c      if ((njets > 0) .and. (countgamm > 0)) then
c        do j=1,countgamm
cc--- Frixione cut, hep-ph/9801442
c          if (njets > 2) then
c            write(6,*) 'Photon-hadron isolation not coded for njets > 2'
c            stop
c          endif
c          do k=1,njets
c            delta(k)=R(pjet,gammindex(j),jetindex(k))
c            ptjet(k)=pt(jetindex(k),pjet)
c          enddo
c          pntr=1
c          if (delta(2) < delta(1)) pntr=2
c          if (njets == 1) delta(2)=gammcone
c          discr=(1._dp-cos(delta(pntr)))/(1._dp-cos(gammcone))
c          if (ptjet(pntr) > discr*pt(gammindex(j),pjet)) then
c            gencuts=.true.
c            return
c          endif
c          if (njets >= 2) then
c            discr=(1._dp-cos(delta(3-pntr)))/(1._dp-cos(gammcone))
c            if (ptjet(1)+ptjet(2) > discr*pt(gammindex(j),pjet))then
c              gencuts=.true.
c              return
c            endif
c          endif

c--- this block was already removed
c--- optional jet-veto
c          if (   (pt(4+countgamm+1,pjet) > 50._dp)
c     &     .and. (abs(etarap(4+countgamm+1,pjet)) < 2.5_dp)) then
c            gencuts=.true.
c          endif
c--- de-Florian,Signer cut
c          do nu=1,2
c            sumjetpt(nu)=0._dp
c          enddo
c          do k=1,njets
c            if (R(pjet,gammindex(j),4+countgamm+k) < gammcone) then
c              do nu=1,2
c              sumjetpt(nu)=sumjetpt(nu)+pjet(4+countgamm+k,nu)
c              enddo
c            endif
c          enddo
c          if ( sqrt(sumjetpt(1)**2+sumjetpt(2)**2)
c     &    > gammcut*pt(gammindex(j),pjet) ) then
c            gencuts=.true.
c          endif
c--- this block was already removed

c        enddo
c      endif
c--- DEBUG: removed all isolation

 99   format(1x,a29,f6.2,a17)
 98   format(' *          ',f8.2,'  <   ',a3,'  < ',f8.2,'           *')
 97   format(' *  ',f9.3,' < ',a18,' < ',f9.3,a4,'  *')
      end




