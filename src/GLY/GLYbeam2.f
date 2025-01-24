!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c#######################################################################
c#### next-to-leading order quark beam functions in qt scheme          #
c#### normalized to (as/(4 pi))^2 which has been extracted             #
c#######################################################################
      subroutine GLYbeam2(ih,zin,xb,beam,ibeam)
        use LHAPDF
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'constants.f'
      include 'transitionlabels.f'
      include 'facscale.f'
      include 'scale.f'
      include 'tiny.f'
      include 'distributions.f'
c     dmin=-2,dmax=1,rglr=-2,delt=-1,plus=0,lpls=1
      real(dp), intent(in) ::zin,xb
      real(dp), intent(out) :: beam(-5:5,0:4)
      real(dp)::cI2(0:6,dmin:dmax,0:4),
     & ccI2(-nf:nf,-nf:nf,dmin:dmax,0:4),
     & fx(-5:5),fx0(-5:5),z,jaco,omz,omxb
      integer :: ih,i,j,ibeam

c--- changing variables so z integral ranges from xb to 1
      omxb=1._dp-xb
      z = omxb*zin+xb
      omz=1._dp-z
      jaco = abs(omxb)

c--- catch numerical problems
      if (z > one-tiny) then
        beam(:,:)=zip
        return
      endif

c calculate parton distribution
      call fdist(ih,xb,facscale,fx0,ibeam)     ! at xb
      call fdist(ih,xb/z,facscale,fx,ibeam)   ! at xb/z

c Lw = log(nu/w) = log(mu/w) - log(mu/nu) = log(mu/w) - Lnu
c where w = x*sqrts
      call GLYI2(z,cI2)

c---------------------------------------------
c-----Translate to ij transition notation ----
c---------------------------------------------
      do i=-nf,nf
      do j=-nf,nf
      if ((i == 0) .and. (j == 0)) then
      ccI2(i,j,:,:)=cI2(gg,:,:)
      elseif ((i == 0) .and. (j  /=  0)) then
      ccI2(i,j,:,:)=cI2(gq,:,:)
      elseif ((i  /=  0) .and. (j == 0)) then
      ccI2(i,j,:,:)=cI2(qg,:,:)
      elseif (i  ==  j) then
      ccI2(i,j,:,:)=cI2(qq,:,:)
      elseif (-i  ==  j) then
      ccI2(i,j,:,:)=cI2(qbq,:,:)
      else
      ccI2(i,j,:,:)=cI2(qpq,:,:)
      endif
      enddo
      enddo

c--------------------------------------------------------
c----- Calculation of quark and gluon beam functions ----
c--------------------------------------------------------

c JC: 11/23 -- extended this to i=0 to try to account for gluon case in the same way

      do i=0,nf
      beam(i,:)=0
      beam(-i,:)=0
      do j=-nf,nf
      beam(i,:)=beam(i,:)
     &  +ccI2(i,j,delt,:)*fx0(j)                       ! delta function
      beam(i,:)=beam(i,:)
     &  +ccI2(i,j,plus,:)*((fx(j)/z-fx0(j))/omz*jaco
     &  +fx0(j)*log(omxb))                                 ! plus piece
      beam(i,:)=beam(i,:)
     &  +ccI2(i,j,lpls,:)*((fx(j)/z-fx0(j))*log(omz)/omz*jaco
     &  +fx0(j)*0.5_dp*log(omxb)**2)                        ! log plus piece
      beam(i,:)=beam(i,:)
     &  +ccI2(i,j,rglr,:)*fx(j)/z*jaco                      ! no distribution piece

      if (i > 0) then
      beam(-i,:)=beam(-i,:)
     &  +ccI2(-i,j,delt,:)*fx0(j)                       ! delta function
      beam(-i,:)=beam(-i,:)
     &  +ccI2(-i,j,plus,:)*((fx(j)/z-fx0(j))/omz*jaco
     &  +fx0(j)*log(omxb))                                 ! plus piece
      beam(-i,:)=beam(-i,:)
     &  +ccI2(-i,j,lpls,:)*((fx(j)/z-fx0(j))*log(omz)/omz*jaco
     &  +fx0(j)*0.5_dp*log(omxb)**2)                        ! log plus piece
      beam(-i,:)=beam(-i,:)
     &  +ccI2(-i,j,rglr,:)*fx(j)/z*jaco                      ! no distribution piece
      endif
      enddo
      enddo

      return
      end


