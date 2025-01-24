!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function passedcuts_wgammajet_ew(isub,p)
      implicit none
      include 'types.f'
      logical:: passedcuts_wgammajet_ew
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
      include 'mxpart.f'
      include 'leptcuts.f'
      include 'jetcuts.f'
      include 'mpicommon.f'
      include 'npart.f'
      include 'nwz.f'
!      include 'hybridiso.f'
      logical:: passed_iso,unobs(2)
      integer:: i,j,inu,isub,nunobs,ilep,countgamm
      real(dp):: p(mxpart,4),pt,yrap,Rpure,ptpure,yrappure,mtrans,
     & plep(4),pgamm1(4),pgamm2(4),pmet(4),Rl1,Rl2
      real(dp):: ptgamm1,ptlep,mtranselgam
      common/ew_observables/ptgamm1,ptlep,mtranselgam
!$omp threadprivate(/ew_observables/)

      passedcuts_wgammajet_ew=.false.

      ptgamm1=0._dp
      ptlep=0._dp
      mtranselgam=0._dp

c Work out how many photons are in the event
      countgamm=npart-3-isub

c Identify lepton and neutrino
      if (nwz == +1) then
        inu=3
        ilep=4
      else
        inu=4
        ilep=3
      endif
      pmet(:)=p(inu,:)

c First ensure all photons are isolated
c       if (hybridiso) then
c         call do_hybridiso(p,passed_iso,5,isub)
c       else
c         call frix(p,passed_iso,5,isub)
c       endif
c       if (passed_iso .eqv. .false.) return
c       if (countgamm == 2) then
c         if (hybridiso) then
c           call do_hybridiso(p,passed_iso,7,isub)
c         else
c           call frix(p,passed_iso,7,isub)
c         endif
c         if (passed_iso .eqv. .false.) return
c       endif

c Recombine lepton and photons
      if (countgamm == 2) then
        Rl1=Rpure(p(ilep,:),p(5,:))
        Rl2=Rpure(p(ilep,:),p(7,:))
        if (Rl1 < Rl2) then
          pgamm1(:)=p(5,:)
          pgamm2(:)=p(7,:)
        else
          pgamm1(:)=p(7,:)
          pgamm2(:)=p(5,:)
          Rl1=Rl2
        endif
        if (Rl1 < 0.1_dp) then
          plep(:)=p(ilep,:)+pgamm1(:)
          pgamm1(:)=pgamm2(:)
          countgamm=countgamm-1
        else
          plep(:)=p(ilep,:)
        endif
        Rl2=Rpure(plep,pgamm2)
        if (Rl2 < 0.1_dp) then
          plep(:)=plep(:)+pgamm2(:)
          countgamm=countgamm-1
        endif
        if (countgamm == 0) return ! all photons recombined into lepton
      else
        plep(:)=p(ilep,:)
        pgamm1(:)=p(5,:)
        Rl1=Rpure(plep,pgamm1)
        if (Rl1 < 0.1_dp) return ! photon recombined into lepton
      endif

c Combine photons
      if (countgamm == 2) then
        if (Rpure(pgamm1,pgamm2) < 0.1_dp) then
          pgamm1(:)=pgamm1(:)+pgamm2(:)
          countgamm=countgamm-1
        endif
      endif

c      write(6,*) 'checking pt, rap, etc.'
c Need to implement isolation here ....
      if (ptpure(plep) < leptptmin) return  ! lepton pt
      if (abs(yrappure(plep)) > leptrapmax) return  ! lepton rapidity
      if (ptpure(pmet) < misspt) return  ! missing pt

c Check how many photons are actually observed
c (note that for Pasold et al, gammptmin = 15)
      nunobs=0
      unobs(:)=.false.
      if ((ptpure(pgamm1) < gammptmin) .or. (abs(yrappure(pgamm1)) > gammrapmax)
     & .or. (Rpure(plep,pgamm1) < Rgalmin)) then
        nunobs=nunobs+1
        unobs(1)=.true.
      endif
      if (countgamm == 2) then
        if ((ptpure(pgamm2) < gammptmin) .or. (abs(yrappure(pgamm2)) > gammrapmax)
     &   .or. (Rpure(plep,pgamm2) < Rgalmin)) then
          nunobs=nunobs+1
          unobs(2)=.true.
        endif
      endif
c      write(6,*) 'checking nunobs',countgamm,nunobs
      if (countgamm - nunobs < 1) return

c      write(6,*) 'checking MT'
c Transverse mass cut (M_T(l,nu))
      mtrans=(plep(1)*pmet(1)+plep(2)*pmet(2))/sqrt((plep(1)**2+plep(2)**2)
     &       *(pmet(1)**2+pmet(2)**2))
      mtrans=2._dp*sqrt(plep(1)**2+plep(2)**2)*sqrt(pmet(1)**2+pmet(2)**2)
     &       *(1._dp-mtrans)
      mtrans=sqrt(max(mtrans,zip))
      if (mtrans < mtrans34cut)  return

c! Now check isolation condition for photons
      do i=1,countgamm
        if (unobs(i) .eqv. .false.) then
          if (i == 1) then
            if (abs(pgamm1(4)/p(5,4)-1._dp) < 1.e-8_dp) then
!              if (hybridiso) then
!                call do_hybridiso(p,passed_iso,5,isub)
!              else
                call frix(p,passed_iso,5,isub)
!              endif
            else
!              if (hybridiso) then
!                call do_hybridiso(p,passed_iso,7,isub)
!              else
                call frix(p,passed_iso,7,isub)
!              endif
            endif
c            if (Rpure(p(6,:),pgamm1) < Rgajetmin) then
c              passed_iso=.false.
c            else
c              passed_iso=.true.
c            endif
          else
            if (abs(pgamm2(4)/p(5,4)-1._dp) < 1.e-8_dp) then
!              if (hybridiso) then
!                call do_hybridiso(p,passed_iso,5,isub)
!              else
                call frix(p,passed_iso,5,isub)
!              endif
            else
!              if (hybridiso) then
!                call do_hybridiso(p,passed_iso,7,isub)
!              else
                call frix(p,passed_iso,7,isub)
!              endif
            endif
c            if (Rpure(p(6,:),pgamm2) < Rgajetmin) then
c              passed_iso=.false.
c            else
c              passed_iso=.true.
c            endif
          endif
          if (passed_iso .eqv. .false.) nunobs=nunobs+1
        endif
      enddo
      if (countgamm - nunobs < 1) return

c      call frix(p,passed_iso,5,isub)
c      if (countgamm == 2) then
c        call frix(p,passed_iso2,7,isub)
c        if ((passed_iso .eqv. .false.) .and. (passed_iso2 .eqv. .false.)) return
c      else
c        if (passed_iso .eqv. .false.) return
c      endif

c Check jet cuts
      if ((pt(6,p) < ptjetmin) .or. (abs(yrap(6,p)) > etajetmax)) return

      passedcuts_wgammajet_ew = .true.

c fill quantities to plot
      ptlep=ptpure(plep)
      if (countgamm == 1) then
        ptgamm1=ptpure(pgamm1)
      else
        if (nunobs == 0) then
          if (ptpure(pgamm1) > ptpure(pgamm2)) then
            ptgamm1=ptpure(pgamm1)
          else
            ptgamm1=ptpure(pgamm2)
            pgamm1(:)=pgamm2(:)
          endif
        else
          if ((ptpure(pgamm1) < gammptmin) .or. (abs(yrappure(pgamm1)) > gammrapmax)
     &     .or. (Rpure(plep,pgamm1) < Rgalmin)) then
            ptgamm1=ptpure(pgamm2)
            pgamm1(:)=pgamm2(:)
          else
            ptgamm1=ptpure(pgamm1)
          endif
        endif
      endif
      mtranselgam=(plep(4)+pgamm1(4))**2-(plep(1)+pgamm1(1))**2
     &           -(plep(2)+pgamm1(2))**2-(plep(3)+pgamm1(3))**2
      mtranselgam=mtranselgam+(plep(1)+pgamm1(1))**2+(plep(2)+pgamm1(2))**2
      mtranselgam=sqrt(max(mtranselgam,zip))+sqrt(pmet(1)**2+pmet(2)**2)
      mtranselgam=mtranselgam**2
      do j=1,2
         mtranselgam=mtranselgam-(plep(j)+pgamm1(j)+pmet(j))**2
      enddo
      mtranselgam=sqrt(max(mtranselgam,zip))

      return
      end




