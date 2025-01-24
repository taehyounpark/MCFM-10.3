!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine resummation_BNRptlumxmsq(p,xx,z1,z2,ptveto,Q,R,msq,order,xmsq)
!     Calculate Luminosity in appropriate order
!     hard_routine is the hard routine for the appropriate order
!     p is the four-momentum
!     xx(2) are the Born x values
!     z1,z2 are integration variables
!     ptveto is the value of the ptcut
!     Q^2=x1*x2*S
!     R is jet cut variable      
! note that, for this suite of routines
!  mu = facscale
! muh = scale
      use LHAPDF
      use SCET
c Wrapper routine for new implementation of BNR resummation
      implicit none
      include 'types.f'
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

      integer:: order,m,ii
      real(dp),intent(in)::ptveto,R,Q
      real(dp),intent(out):: xmsq
      real(dp):: p(mxpart,4),xx(2),
     & tBa0(-5:5),tBb0(-5:5),
     & tBa1(-5:5),tBb1(-5:5),
     & tBa2(-5:5),tBb2(-5:5),
     & z1,z2,lnmuonpt,
     & msq(-nf:nf,-nf:nf),
     & resummation_BNRptgetxmsq,Fanom(3),h(3),ptvetovar

      if (first) then
        first=.false.
        call Gammafill(nflav)
      endif

! Special case for ipsgen = 2: gg contributions only at fixed order
! (with extension to NLL also anticipated)
      if (ipsgen == 2) then
        call fdist(ih1,xx(1),facscale,tBa0,1)
        call fdist(ih2,xx(2),facscale,tBb0,2)
        xmsq = tBa0(0)*tBb0(0)*msq(0,0)
        lnmuonpt=log(facscale/ptveto)
!        return    ! stick to fixed order for now
        call Ffill(0,lnmuonpt,R,Fanom)
        h(:)=0._dp
        xmsq = resummation_BNRptgetxmsq(p,0,0,ptveto,Q,Fanom,h,
     &              tBa0,tBb0,tBa1,tBb1,tBa2,tBb2,msq)
        return
      endif

! In the resummed calculation all logs are at the scale ~ ptveto (facscale)
      lnmuonpt=log(facscale/ptveto)

c to include gg->h
      if (kcase==kggfus0) then
        ii=0
      else
        ii=1
      endif
      call Ffill(ii,lnmuonpt,R,Fanom)
      call hfill(ii,lnmuonpt,h)
      
      tBa1(:) = 0._dp
      tBb1(:) = 0._dp
      tBa2(:) = 0._dp
      tBb2(:) = 0._dp

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

      xmsq = resummation_BNRptgetxmsq(p,order,ii,ptveto,Q,Fanom,h,
     &              tBa0,tBb0,tBa1,tBb1,tBa2,tBb2,msq)
     
      return
      end


