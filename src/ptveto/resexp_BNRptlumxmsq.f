!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine resexp_BNRptlumxmsq(p,xx,z1,z2,ptveto,Q,R,order,msq0,msq1,msq2,xmsq)
!     Calculate Luminosity in appropriate order
!     p is the four-momentum
!     msq0, msq1, msq2 are the matrix elements at 0-, 1- and 2-loop order (msq2 should be zero)
!     xx(2) are the Born x values
!     z1,z2 are integration variables
!     ptveto is the value of the ptcut
!     Q^2=x1*x2*S
!     R is jet cut variable
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

      integer:: order,m
      real(dp),intent(in)::ptveto,R,Q
      real(dp),intent(out):: xmsq
      real(dp):: p(mxpart,4),xx(2),
     & tBa0(-5:5),tBb0(-5:5),
     & tBa1(-5:5),tBb1(-5:5),
     & tBa2(-5:5),tBb2(-5:5),
     & z1,z2,lnmuonpt,
     & msq0(-nf:nf,-nf:nf),msq1(-nf:nf,-nf:nf),msq2(-nf:nf,-nf:nf),
     & BNRptgetxmsq,Fanom(3),h(3)

      if (first) then
        first=.false.
        call Gammafill(nflav)
      endif

      lnmuonpt=log(scale/ptveto)

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


      return
      end

