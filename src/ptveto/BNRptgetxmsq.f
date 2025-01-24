!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function BNRptgetxmsq(p,order,ptveto,Q,Fanom,h,
     &              tBa0,tBb0,tBa1,tBb1,tBa2,tBb2,msq0,msq1,msq2)
      use SCET
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'kprocess.f'
      include 'kpart.f'
      include 'nflav.f'
      include 'facscale.f'
      include 'scale.f'
      real(dp):: BNRptgetxmsq,ptgetxmsq
      real(dp), intent(in):: p(mxpart,4),ptveto,Q,
     & msq0(-nf:nf,-nf:nf),msq1(-nf:nf,-nf:nf),msq2(-nf:nf,-nf:nf),
     & Fanom(3),h(3),
     & tBa0(-5:5),tBb0(-5:5),
     & tBa1(-5:5),tBb1(-5:5),
     & tBa2(-5:5),tBb2(-5:5)
      integer:: j,k,order
!      real(dp):: collanom,hfac
      real(dp):: hard(2),xlum,lnQonpt,collanomexp(2),hfacexp(2),beta0

!      collanom=-2*(
!     & ason4pi*(Fanom(1)+ason4pi*(Fanom(2)+ason4pi*Fanom(3))))
!      collanom=(Q/ptveto)**collanom
!      hfac=2*ason4pi*(h(1)+ason4pi*(h(2)+ason4pi*h(3)))
!      hfac=exp(hfac)

      lnQonpt=log(Q/ptveto)
      collanomexp(1) = -2*Fanom(1)*lnQonpt
      collanomexp(2) = 2*Fanom(1)**2*lnQonpt**2-2*Fanom(2)*lnQonpt
      hfacexp(1) = 2*h(1)
      hfacexp(2) = 2*h(1)**2+2*h(2)

      beta0=(11*CA-4*TF*nflav)/3d0

      ptgetxmsq=zip
      do j=-nflav,nflav
      do k=-nflav,nflav

c WARNING: this is meant to catch gg contributions that only enter for
c the first time at NNLO, i.e. everything so far except for gg->h process
      if ((j == 0) .and. (k == 0) .and. (order == 2) .and. (kcase /= kggfus0)) then
        ptgetxmsq = ptgetxmsq + tBa0(j)*tBb0(k)*ason4pi**2*msq2(j,k)
        cycle
      endif

      if (msq0(j,k) == zip) cycle

      hard(1)=msq1(j,k)/msq0(j,k)
      hard(2)=msq2(j,k)/msq0(j,k)

      if (coeffonly) then
         xlum=zip
      else
         xlum=tBa0(j)*tBb0(k)
      endif
      if ( (order == 1) .or.
     &      ((order == 2) .and. (coeffonly .eqv. .false.)) ) then
         xlum=xlum+ason4pi*(
     &    +hard(1)*tBa0(j)*tBb0(k)+tBa0(j)*tBb1(k)+tBa1(j)*tBb0(k)
     &    +(hfacexp(1)+collanomexp(1))*tBa0(j)*tBb0(k))
      endif
      if (order > 1) then
         xlum=xlum+ason4pi**2*(
     &    +hard(2)*tBa0(j)*tBb0(k)
     &    +hard(1)*(tBa0(j)*tBb1(k)+tBa1(j)*tBb0(k))
     &    +tBa0(j)*tBb2(k)+tBa2(j)*tBb0(k)
     &    +tBa1(j)*tBb1(k)
     &    +(hfacexp(1)+collanomexp(1))*(tBa1(j)*tBb0(k)+tBa0(j)*tBb1(k))
     &    +(hfacexp(1)+collanomexp(1))*hard(1)*tBa0(j)*tBb0(k)
     &    +(hfacexp(1)*collanomexp(1))*tBa0(j)*tBb0(k)
     &    +(hfacexp(2)+collanomexp(2))*tBa0(j)*tBb0(k)
! Additional term to account for evaluating beam function corrections
! using as(scale) rather than as(facscale)
     &    +2*beta0*log(scale/facscale)*(
     &    +tBa0(j)*tBb1(k)+tBa1(j)*tBb0(k)))
      endif

      ptgetxmsq=ptgetxmsq+xlum*msq0(j,k)

      enddo
      enddo
      BNRptgetxmsq=ptgetxmsq
      return
      end
