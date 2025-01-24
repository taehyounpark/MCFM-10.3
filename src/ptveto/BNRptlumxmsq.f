!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine BNRptlumxmsq(p,xx,z1,z2,ptveto,Q,R,order,central,
     & hard_routine,xmsq)
!     Calculate Luminosity in appropriate order
!     hard_routine is the hard routine for the appropriate order
!     p is the four-momentum
!     xx(2) are the Born x values
!     z1,z2 are integration variables
!     ptveto is the value of the ptcut
!     Q^2=x1*x2*S
!     R is jet cut variable
      use LHAPDF
      use SCET
      use MCFMStorage, only: currentPart
      use Integration
c Wrapper routine for new implementation of BNR resummation
      implicit none
!      include 'types.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'first.f'
      include 'scale.f'
      include 'facscale.f'
      include 'nflav.f'
      include 'kprocess.f'
      include 'beamtype.f'
      include 'taucut.f'
      include 'ipsgen.f'

      integer:: order,m
      real(dp),intent(in)::ptveto,R,Q
      real(dp),intent(out):: xmsq
      real(dp):: p(mxpart,4),xx(2),
     & tBa0(-5:5),tBb0(-5:5),
     & tBa1(-5:5),tBb1(-5:5),
     & tBa2(-5:5),tBb2(-5:5),
     & z1,z2,lnmuonpt,
     & msq0(-nf:nf,-nf:nf),msq1(-nf:nf,-nf:nf),msq2(-nf:nf,-nf:nf),
     & BNRptgetxmsq,Fanom(3),h(3),ptvetovar

      logical, intent(in) :: central
      real(dp) :: origtaucut
      external hard_routine

      if (first) then
        first=.false.
        call Gammafill(nflav)
      endif

      lnmuonpt=log(scale/ptveto)

      if ((currentpart == nnloResVetoBelow) .and. (ipsgen == 1)) then
! for the pt-veto case, this is the matching correction and msq2 is not
! included in the fixed-order expansion of the resummed result
        call hard_routine(p,1,msq0,msq1,msq2)
        msq2(:,:)=0._dp
      else
! regular case
        call hard_routine(p,order,msq0,msq1,msq2)
      endif

c to include gg->h
      if (kcase==kggfus0) then
        call Ffill(0,lnmuonpt,R,Fanom)
        call hfill(0,lnmuonpt,h)
      else
        call Ffill(1,lnmuonpt,R,Fanom)
        call hfill(1,lnmuonpt,h)
      endif

      lnmuonpt=log(facscale/ptveto)

      if (order >= 0) then
        call fdist(ih1,xx(1),facscale,tBa0,1)
        call fdist(ih2,xx(2),facscale,tBb0,2)
      endif
      if (order >= 1) then
        call BNRptbeam1(ih1,z1,xx(1),lnmuonpt,tBa1,1)
        call BNRptbeam1(ih2,z2,xx(2),lnmuonpt,tBb1,2)
      endif
      if (order >= 2) then
        call BNRptbeam2(ih1,z1,xx(1),lnmuonpt,R,tBa2,1)
        call BNRptbeam2(ih2,z2,xx(2),lnmuonpt,R,tBb2,2)
      endif

      xmsq = BNRptgetxmsq(p,order,ptveto,Q,Fanom,h,
     &              tBa0,tBb0,tBa1,tBb1,tBa2,tBb2,msq0,msq1,msq2)
     

! This section is a repetition of the above for now;
! could be done more elegantly
      if (central .and. doMultitaucut) then
          scetreweight(:) = 0._dp
          if (xmsq /= 0._dp) then
              origtaucut = taucut
              do m=1,size(tcutarray)
                  taucut = tcutarray(m)
                  if (dynamictau) then
                    ptvetovar=getdynamictau(p)
                  else
                    ptvetovar=taucut
                  endif

                  lnmuonpt=log(scale/ptvetovar)

                  if (kcase==kggfus0) then
                    call Ffill(0,lnmuonpt,R,Fanom)
                    call hfill(0,lnmuonpt,h)
                  else
                    call Ffill(1,lnmuonpt,R,Fanom)
                    call hfill(1,lnmuonpt,h)
                  endif

                  lnmuonpt=log(facscale/ptvetovar)

                  if (order >= 0) then
                    call fdist(ih1,xx(1),facscale,tBa0,1)
                    call fdist(ih2,xx(2),facscale,tBb0,2)
                  endif
                  if (order >= 1) then
                    call BNRptbeam1(ih1,z1,xx(1),lnmuonpt,tBa1,1)
                    call BNRptbeam1(ih2,z2,xx(2),lnmuonpt,tBb1,2)
                  endif
                  if (order >= 2) then
                    call BNRptbeam2(ih1,z1,xx(1),lnmuonpt,R,tBa2,1)
                    call BNRptbeam2(ih2,z2,xx(2),lnmuonpt,R,tBb2,2)
                  endif

                  scetreweight(m) = BNRptgetxmsq(p,order,ptvetovar,Q,Fanom,h,
     &              tBa0,tBb0,tBa1,tBb1,tBa2,tBb2,msq0,msq1,msq2)

              enddo
              taucut = origtaucut
              scetreweight(:) = scetreweight(:) / xmsq
          endif
      endif

      return
      end


