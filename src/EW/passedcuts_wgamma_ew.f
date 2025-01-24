!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function passedcuts_wgamma_ew(isub,p)
      implicit none
      include 'types.f'
      logical:: passedcuts_wgamma_ew
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
      include 'mpicommon.f'
      include 'npart.f'
      include 'nwz.f'
      integer:: j,inu,isub,nunobs,ilep,countgamm
      real(dp):: p(mxpart,4),Rpure,ptpure,yrappure,mtrans,
     & plep(4),pgamm1(4),pgamm2(4),pmet(4),Rl1,Rl2
      real(dp):: ptgamm1,ptlep,mtranselgam,ygamm1,ylepgamm1,Rlepgamm1
      common/ew_observables/ptgamm1,ptlep,mtranselgam
      common/ew_observables_extra/ygamm1,ylepgamm1,Rlepgamm1
!$omp threadprivate(/ew_observables/)
!$omp threadprivate(/ew_observables_extra/)

      passedcuts_wgamma_ew=.false.

      ptgamm1=0._dp
      ptlep=0._dp
      mtranselgam=0._dp

c Work out how many photons are in the event
      countgamm=npart-2-isub

c Identify lepton and neutrino
      if (nwz == +1) then
        inu=3
        ilep=4
      else
        inu=4
        ilep=3
      endif
      pmet(:)=p(inu,:)

c Recombine lepton and photons
      if (countgamm == 2) then
        Rl1=Rpure(p(ilep,:),p(5,:))
        Rl2=Rpure(p(ilep,:),p(6,:))
        if (Rl1 < Rl2) then
          pgamm1(:)=p(5,:)
          pgamm2(:)=p(6,:)
        else
          pgamm1(:)=p(6,:)
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
c (note that for wgamma_ew et al, gammptmin = 15)
      nunobs=0
      if ((ptpure(pgamm1) < gammptmin) .or. (abs(yrappure(pgamm1)) > gammrapmax)
     & .or. (Rpure(plep,pgamm1) < Rgalmin)) nunobs=nunobs+1
      if (countgamm == 2) then
        if ((ptpure(pgamm2) < gammptmin) .or. (abs(yrappure(pgamm2)) > gammrapmax)
     &   .or. (Rpure(plep,pgamm2) < Rgalmin)) nunobs=nunobs+1
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

      passedcuts_wgamma_ew = .true.

c fill quantities to plot
      ptlep=ptpure(plep)
      if (countgamm == 1) then
        ptgamm1=ptpure(pgamm1)
        ygamm1=yrappure(pgamm1)
        Rlepgamm1=Rpure(plep,pgamm1)
      else
        if (nunobs == 0) then
          if (ptpure(pgamm1) > ptpure(pgamm2)) then
            ptgamm1=ptpure(pgamm1)
           ygamm1=yrappure(pgamm1)
           Rlepgamm1=Rpure(plep,pgamm1)
          else
            ptgamm1=ptpure(pgamm2)
            pgamm1(:)=pgamm2(:)
            ygamm1=yrappure(pgamm1)
            Rlepgamm1=Rpure(plep,pgamm1)
          endif
        else
          if ((ptpure(pgamm1) < gammptmin) .or. (abs(yrappure(pgamm1)) > gammrapmax)
     &     .or. (Rpure(plep,pgamm1) < Rgalmin)) then
            ptgamm1=ptpure(pgamm2)
            pgamm1(:)=pgamm2(:)
            ygamm1=yrappure(pgamm1)
            Rlepgamm1=Rpure(plep,pgamm1)
          else
            ptgamm1=ptpure(pgamm1)
            ygamm1=yrappure(pgamm1)
            Rlepgamm1=Rpure(plep,pgamm1)
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
      ylepgamm1=nwz*(ygamm1-yrappure(plep))

      return
      end




