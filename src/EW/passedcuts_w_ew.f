!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function passedcuts_w_ew(isub,p)
      implicit none
      include 'types.f'
      logical:: passedcuts_w_ew
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
c      include 'nf.f'
      include 'mxpart.f'
c      include 'cplx.h'
c      include 'limits.f'
      include 'kprocess.f'
      include 'leptcuts.f'
c      include 'z_dip.f'
c      include 'jetlabel.f'
      include 'mpicommon.f'
c      include 'runstring.f'
      include 'npart.f'
      include 'nwz.f'
      integer:: j,inu,isub,ilep,countgamm
      real(dp):: p(mxpart,4),Rpure,ptpure,yrappure,mtrans,
     & plep(4),pgamm1(4),pmet(4),Rl1
      real(dp):: mtranselgam
      real(dp):: mln,puremass
      common/w_ew_observables/mln
!$omp threadprivate(/w_ew_observables/)

      passedcuts_w_ew=.false.

      mtranselgam=0._dp
      mln=1._dp

c Work out how many photons are in the event
      countgamm=npart-2-isub
      if (kcase == kWln_aq) countgamm=0

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
      plep(:)=p(ilep,:)
      pgamm1(:)=p(5,:)
c Optional recombination (only necessary for distributions)
      if (countgamm == 1) then
        Rl1=Rpure(plep,pgamm1)
        if (Rl1 < 0.1_dp) then
          plep(:)=plep(:)+pgamm1(:)
          pgamm1(:)=0._dp
          countgamm=countgamm-1
        elseif (Rl1 < 0.4_dp) then
          if (pgamm1(4) > 0.1_dp*plep(4)) return ! cut from 1606.02330
        endif
      endif

c      write(6,*) 'checking pt, rap, etc.'
c Need to implement isolation here ....
      if (ptpure(plep) < leptptmin) return  ! lepton pt
      if (abs(yrappure(plep)) > leptrapmax) return  ! lepton rapidity
      if (ptpure(pmet) < misspt) return  ! missing pt

c      write(6,*) 'checking MT'
c Transverse mass cut (M_T(l,nu))
      mtrans=(plep(1)*pmet(1)+plep(2)*pmet(2))/sqrt((plep(1)**2+plep(2)**2)
     &       *(pmet(1)**2+pmet(2)**2))
      mtrans=2._dp*sqrt(plep(1)**2+plep(2)**2)*sqrt(pmet(1)**2+pmet(2)**2)
     &       *(1._dp-mtrans)
      mtrans=sqrt(max(mtrans,zip))
      if (mtrans < mtrans34cut)  return

      passedcuts_w_ew = .true.

c fill quantities to plot
      mln=puremass(plep(:)+pmet(:))
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




