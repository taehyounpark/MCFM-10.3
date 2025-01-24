!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function resummation_BNRptgetxmsq(p,order,ii,ptveto,Q,Fanomin,hin,
     &              tBa0,tBb0,tBa1,tBb1,tBa2,tBb2,msq)
! note that, for this suite of routines
!  mu = facscale
! muh = scale
      use SCET
      use LHAPDF, only: getalphas
      use qtResummation_params, only : scalevar_rapidity_mult, scalevar_rapidity_mult_higgs,
     &                                 scalevar_rapidity_i
      use ptveto, only: gghsinglestep
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'kprocess.f'
      include 'kpart.f'
      include 'facscale.f'
      include 'scale.f'
      include 'masses.f'
      include 'nfl.f'
      include 'taucut.f'
      real(dp):: resummation_BNRptgetxmsq,Ct1,Ct2,Ctsq,hardevol
      real(dp), intent(in):: p(mxpart,4),ptveto,Q,
     & msq(-nf:nf,-nf:nf),
     & Fanomin(3),hin(3),
     & tBa0(-5:5),tBb0(-5:5),
     & tBa1(-5:5),tBb1(-5:5),
     & tBa2(-5:5),tBb2(-5:5)
      integer:: j,k,order,ii
      real(dp):: collanom,hfac,xlum,lnQonpt,Fanom(3),h(3) !,collanomexp(2),hfacexp(2)
      real(dp):: alphasmu,ason4pi,logmt2omu2,alphasptveto,alphasmut,beta0,beta1,beta2
      real(dp):: collanomexp(2),logR

! Compute alphas at the resummation scale (facscale, close to ptveto)
      alphasmu=getalphas(facscale)
      ason4pi=alphasmu/fourpi
      
! only fill entries up to the order required
!      order = 0  NLL
!      order = 1  NNLL
!      order = 2  N3LL

      Fanom(:)=0._dp
      h(:)=0._dp
! keep one extra power since large log
      do k=1,min(order+1,3)
        Fanom(k)=Fanomin(k)
      enddo
! usual expansion
      if (order > 0) then
        do k=1,min(order,3)
          h(k)=hin(k)
        enddo
      endif
!      write(6,*) 'Fanom',Fanom
!      write(6,*) 'h    ',h

! rapidity scale variation according to Eq. (67) of 1307.0025
      if (scalevar_rapidity_i > 0) then
! much smaller range for gluon-fusion Higgs process
          if (ii == 0) then
            logR = log(scalevar_rapidity_mult_higgs(scalevar_rapidity_i))
          else
            logR = log(scalevar_rapidity_mult(scalevar_rapidity_i))
          endif
      else
          logR = 0._dp
      endif
      collanomexp(1) = -2*Fanom(1)*logR
      collanomexp(2) = 2*Fanom(1)**2*logR**2-2*Fanom(2)*logR

!      logR=0._dp
!      collanomexp(1)=0._dp
!      collanomexp(2)=0._dp

      collanom=-2*(
     & ason4pi*(Fanom(1)+ason4pi*(Fanom(2)+ason4pi*Fanom(3))))
      collanom=(Q/ptveto/exp(logR))**collanom
      hfac=2*ason4pi*(h(1)+ason4pi*(h(2)+ason4pi*h(3)))
      hfac=exp(hfac)

!      lnQonpt=log(Q/ptveto)
!      collanomexp(1) = -2*Fanom(1)*lnQonpt
!      collanomexp(2) = 2*Fanom(1)**2*lnQonpt**2-2*Fanom(2)*lnQonpt
!      hfacexp(1) = 2*h(1)
!      hfacexp(2) = 2*h(1)**2+2*h(2)

      resummation_BNRptgetxmsq=zip
      do j=-nf,nf
      do k=-nf,nf

      if (msq(j,k) == zip) cycle

      if (order == 0) then
        xlum=tBa0(j)*tBb0(k)
      elseif (order == 1) then
        xlum=tBa0(j)*tBb0(k)
     &      +ason4pi*(tBa0(j)*tBb1(k)+tBa1(j)*tBb0(k)
     &               +collanomexp(1)*tBa0(j)*tBb0(k))
      else
        xlum=tBa0(j)*tBb0(k)
     &      +ason4pi*(tBa0(j)*tBb1(k)+tBa1(j)*tBb0(k)
     &               +collanomexp(1)*tBa0(j)*tBb0(k))
     &      +ason4pi**2*(tBa0(j)*tBb2(k)+tBa2(j)*tBb0(k)
     &                  +tBa1(j)*tBb1(k)
     &                  +collanomexp(1)*(tBa1(j)*tBb0(k)+tBa0(j)*tBb1(k))
     &                  +collanomexp(2)*tBa0(j)*tBb0(k))
      endif

      resummation_BNRptgetxmsq=resummation_BNRptgetxmsq+collanom*hfac*xlum*msq(j,k)

      enddo
      enddo

! Note: (ii == 0 )should always correspond to (kcase == kggfus0);
!       if not then there could be problems
      if ((ii == 0) .and. (gghsinglestep .eqv. .false.)) then
! LO resummation_BNRptgetxmsq element has been computed with scale mu, not muf
        resummation_BNRptgetxmsq=resummation_BNRptgetxmsq
     &   *(alphasmu/getalphas(scale))**2
! and the previous lines restore the correct choice
! c.f. a similar correction factor in Eq. (19) of arXiv:1307.0025
!        write(6,*) 'alpha-s rescaling:',(alphasmu/getalphas(scale))**2
      endif

! Now perform hard evolution
      resummation_BNRptgetxmsq=resummation_BNRptgetxmsq*hardevol(ii,order,Q**2,facscale,scale)
!      write(6,*) 'hardevol:',hardevol(ii,order,Q**2,facscale,scale)

!      if (scalevar_rapidity_i > 0) then
!!           resummation_BNRptgetxmsq=resummation_BNRptgetxmsq
!!     &      *scalevar_rapidity_mult_higgs(scalevar_rapidity_i)**(-(ason4pi**2*Fanom(2)))
!           resummation_BNRptgetxmsq=resummation_BNRptgetxmsq
!     &      *scalevar_rapidity_mult_higgs(scalevar_rapidity_i)**(-(ason4pi*Fanom(1)))
!      endif

c--- debug - print value of Hbar
!      if (kcase == kggfus0) then
!        write(6,*) 'Hbar',(alphasmu/getalphas(jetptveto))**2*Ctsq*msq(0,0)*hardevol(ii,order,Q**2,facscale,scale)*collanom*hfac
!      else
!        write(6,*) 'Hbar',msq(1,-1)*hardevol(ii,order,Q**2,facscale,scale)*collanom*hfac
!        write(6,*) 'Hbar',msq(2,-2)*hardevol(ii,order,Q**2,facscale,scale)*collanom*hfac
!      endif
!      stop
c---c debug

      return
      end
