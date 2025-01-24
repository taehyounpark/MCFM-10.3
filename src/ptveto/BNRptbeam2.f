!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c#######################################################################
c#### next-to-leading order quark beam functions in BNR scheme         #
c#### normalized to (as/(4 pi))^2 which has been extracted             #
c#### ih is beam type                                                  #
c#### zin is integration variable                                      #
c#### xb is target value of x                                          #
c#### lnmuonpt=ln(mu/pt)                                               #
c#### R is jet radius                                                  #
c#### beam(-5:5) is output                                             #
c#### ibeam single out either first or second beam particle            #
c#######################################################################
      subroutine BNRptbeam2(ih,zin,xb,lnmuonpt,R,beam,ibeam)
      use LHAPDF
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'constants.f'
      include 'facscale.f'
      include 'tiny.f'
      include 'transitionlabels.f'
      include 'distributions.f'
c     dmin=-2,dmax=1,rglr=-2,delt=-1,plus=0,lpls=1
      real(dp), intent(in) ::zin,xb,lnmuonpt,R
      real(dp), intent(out) :: beam(-5:5)
      real(dp)::I2(0:6,dmin:dmax),ccI2(-nf:nf,-nf:nf,dmin:dmax),
     & fx(-5:5),fx0(-5:5),z,jaco,omz,omxb

      integer :: ih,i,j,ibeam


c--- changing variables so z integral ranges from xb to 1
      omxb=1._dp-xb
      z = omxb*zin+xb
      omz=1._dp-z
      jaco = abs(omxb)

c alternative without rescaling
c      z = zin
c      omz=1._dp-z
c      if (z < xb) then
c        beam(:)=zip
c        return
c      endif

c--- catch numerical problems
      if (z > one-tiny) then
        beam(:)=zip
        return
      endif

c calculate parton distribution
      call fdist(ih,xb,facscale,fx0,ibeam)     ! at xb
      call fdist(ih,xb/z,facscale,fx,ibeam)   ! at xb/z


      call I2bar(z,lnmuonpt,R,I2)


c---------------------------------------------
c-----Translate to ij transition notation ----
c---------------------------------------------
      do i=-nf,nf
      do j=-nf,nf
      if ((i == 0) .and. (j == 0)) then
      ccI2(i,j,:)=I2(gg,:)
      elseif ((i == 0) .and. (j  /=  0)) then
      ccI2(i,j,:)=I2(gq,:)
      elseif ((i  /=  0) .and. (j == 0)) then
      ccI2(i,j,:)=I2(qg,:)
      elseif (i  ==  j) then
      ccI2(i,j,:)=I2(qq,:)
      elseif (i  ==  -j) then
      ccI2(i,j,:)=I2(qbq,:)
      else
!     exploits the fact that qbpq=qpq at this order
      ccI2(i,j,:)=I2(qpq,:)
      endif
      enddo
      enddo

c--------------------------------------------------------
c----- Calculation of quark and gluon beam functions ----
c--------------------------------------------------------


      do i=0,nf
      beam(i)=0
      beam(-i)=0
      do j=-nf,nf
c delta function
      beam(i)=beam(i)+ccI2(i,j,delt)*fx0(j)
c plus piece
      beam(i)=beam(i)
     &  +ccI2(i,j,plus)*((fx(j)/z-fx0(j))/omz*jaco
     &  +fx0(j)*log(omxb))
c log plus piece
      beam(i)=beam(i)
     &  +ccI2(i,j,lpls)*((fx(j)/z-fx0(j))*log(omz)/omz*jaco
     &  +fx0(j)*0.5_dp*log(omxb)**2)
c no distribution piece
      beam(i)=beam(i)+ccI2(i,j,rglr)*fx(j)/z*jaco

      if (i > 0) then
c delta function
      beam(-i)=beam(-i)+ccI2(-i,j,delt)*fx0(j)
c plus piece
      beam(-i)=beam(-i)
     &  +ccI2(-i,j,plus)*((fx(j)/z-fx0(j))/omz*jaco
     &  +fx0(j)*log(omxb))
c log plus piece
      beam(-i)=beam(-i)
     &  +ccI2(-i,j,lpls)*((fx(j)/z-fx0(j))*log(omz)/omz*jaco
     &  +fx0(j)*0.5_dp*log(omxb)**2)
c no distribution piece
      beam(-i)=beam(-i)+ccI2(-i,j,rglr)*fx(j)/z*jaco
      endif
      enddo
      enddo

c alternative without rescaling
c      beam(:,:)=beam(:,:)/jaco

      return
      end


