!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c#######################################################################
c#### next-to-leading order quark beam functions in qt scheme          #
c#### normalized to (as/(4 pi))^2 which has been extracted             #
c#######################################################################
      subroutine ptbeam2(ih,zin,xb,Lperp,beam,ibeam)
        use LHAPDF
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'constants.f'
      include 'transitionlabels.f'
c     gg=0,qqV=1,qbV=2,qqS=3,qg=4,gq=5,qqDS=6
      include 'facscale.f'
      include 'tiny.f'
      include 'Lw.f'
      include 'Lnu.f'
      include 'scale.f'
      include 'energy.f'
      include 'distributions.f'
      include 'Rcut.f'
c     dmin=-2,dmax=1,rglr=-2,delt=-1,plus=0,lpls=1
      real(dp), intent(in) ::zin,xb,Lperp
      real(dp), intent(out) :: beam(-5:5)
      real(dp)::cI2(0:6,dmin:dmax,0:2),
     & ccI2(-nf:nf,-nf:nf,dmin:dmax),
     & fx(-5:5),fx0(-5:5),z,jaco,omz,omxb,LQ,
     & I2p(0:6,dmin:dmax,0:2),
     & I2p0(0:6,dmin:dmax),I2p1(0:6,dmin:dmax),I2pm1(0:6,dmin:dmax),
     & I2perparr(0:6,dmin:dmax)

      integer,parameter::qqV=1,qbV=2,qqS=3
      integer :: ih,i,j,ibeam

      logical, parameter :: check = .false.

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

c LQ = log(nu/w) = log(mu/w) - log(mu/nu) = log(mu/w) - Lnu
c where w = x*sqrts
      LQ = log(scale/(xb*sqrts))-Lnu

      call ptI2(z,Rcut,Lperp,LQ,I2perparr)

! Cross-check
      if (check) then
        call I2perp(z, 0._dp,LQ,I2p0)
        call I2perp(z, 1._dp,LQ,I2p1)
        call I2perp(z,-1._dp,LQ,I2pm1)
        I2p(:,:,0)=I2p0(:,:)
        I2p(:,:,2)=(I2p1(:,:)+I2pm1(:,:)-2*I2p(:,:,0))/2._dp
        I2p(:,:,1)=I2p1(:,:)-I2p(:,:,2)-I2p(:,:,0)
! Lperp is log(mu/pt) but Lb is log(musq/bTsq)
        I2p(:,:,2)=I2p(:,:,2)/4._dp
        I2p(:,:,1)=I2p(:,:,1)/2._dp
        I2p(:,:,0)=I2p(:,:,0)
        Lw=LQ
        call calI2(z,cI2)
        write(6,*) '============================================='
        write(6,*) 'cI2(gg,delt,2)',cI2(gg,delt,2)/I2p(gg,delt,2)
        write(6,*) 'cI2(gg,delt,1)',cI2(gg,delt,1)/I2p(gg,delt,1)
        write(6,*) 'cI2(gg,delt,0)',cI2(gg,delt,0)/I2p(gg,delt,0)

        write(6,*) 'cI2(gg,rglr,2)',cI2(gg,rglr,2)/I2p(gg,rglr,2)
        write(6,*) 'cI2(gg,rglr,1)',cI2(gg,rglr,1)/I2p(gg,rglr,1)
        write(6,*) 'cI2(gg,rglr,0)',cI2(gg,rglr,0)/I2p(gg,rglr,0)

        write(6,*) 'cI2(gg,plus,2)',cI2(gg,plus,2)/I2p(gg,plus,2)
        write(6,*) 'cI2(gg,plus,1)',cI2(gg,plus,1)/I2p(gg,plus,1)
        write(6,*) 'cI2(gg,plus,0)',cI2(gg,plus,0)/I2p(gg,plus,0)

        write(6,*) 'cI2(gg,lpls,2)',cI2(gg,lpls,2)/I2p(gg,lpls,2)

        write(6,*) '============================================='
        write(6,*) 'cI2(gq,rglr,2)',cI2(gq,rglr,2)/I2p(gq,rglr,2)
        write(6,*) 'cI2(gq,rglr,1)',cI2(gq,rglr,1)/I2p(gq,rglr,1)
        write(6,*) 'cI2(gq,rglr,0)',cI2(gq,rglr,0)/I2p(gq,rglr,0)

        write(6,*) '============================================='
        write(6,*) 'cI2(qg,rglr,2)',cI2(qg,rglr,2)/I2p(qg,rglr,2)
        write(6,*) 'cI2(qg,rglr,1)',cI2(qg,rglr,1)/I2p(qg,rglr,1)
        write(6,*) 'cI2(qg,rglr,0)',cI2(qg,rglr,0)/I2p(qg,rglr,0)

        write(6,*) '============================================='
        write(6,*) 'cI2(qq,delt,2)',(cI2(qqV,delt,2)+cI2(qqS,delt,2))/I2p(qq,delt,2)
        write(6,*) 'cI2(qq,delt,1)',(cI2(qqV,delt,1)+cI2(qqS,delt,1))/I2p(qq,delt,1)
        write(6,*) 'cI2(qq,delt,0)',(cI2(qqV,delt,0)+cI2(qqS,delt,0))/I2p(qq,delt,0)

        write(6,*) 'cI2(qq,rglr,2)',(cI2(qqV,rglr,2)+cI2(qqS,rglr,2))/I2p(qq,rglr,2)
        write(6,*) 'cI2(qq,rglr,1)',(cI2(qqV,rglr,1)+cI2(qqS,rglr,1))/I2p(qq,rglr,1)
        write(6,*) 'cI2(qq,rglr,0)',(cI2(qqV,rglr,0)+cI2(qqS,rglr,0))/I2p(qq,rglr,0)

        write(6,*) 'cI2(qq,plus,2)',(cI2(qqV,plus,2)+cI2(qqS,plus,2))/I2p(qq,plus,2)
        write(6,*) 'cI2(qq,plus,1)',(cI2(qqV,plus,1)+cI2(qqS,plus,1))/I2p(qq,plus,1)
        write(6,*) 'cI2(qq,plus,0)',(cI2(qqV,plus,0)+cI2(qqS,plus,0))/I2p(qq,plus,0)

        write(6,*) 'cI2(qq,lpls,2)',(cI2(qqV,lpls,2)+cI2(qqS,lpls,2))/I2p(qq,lpls,2)

        write(6,*) '============================================='
        write(6,*) 'cI2(qbq,rglr,2)',(cI2(qbV,rglr,2)+cI2(qqS,rglr,2))/I2p(qbq,rglr,2)
        write(6,*) 'cI2(qbq,rglr,1)',(cI2(qbV,rglr,1)+cI2(qqS,rglr,1))/I2p(qbq,rglr,1)
        write(6,*) 'cI2(qbq,rglr,0)',(cI2(qbV,rglr,0)+cI2(qqS,rglr,0))/I2p(qbq,rglr,0)

        write(6,*) '============================================='
        write(6,*) 'cI2(qpq,rglr,2)',cI2(qqS,rglr,2)/I2p(qpq,rglr,2)
        write(6,*) 'cI2(qpq,rglr,1)',cI2(qqS,rglr,1)/I2p(qpq,rglr,1)
        write(6,*) 'cI2(qpq,rglr,0)',cI2(qqS,rglr,0)/I2p(qpq,rglr,0)

        pause
      endif

c---------------------------------------------
c-----Translate to ij transition notation ----
c---------------------------------------------
      do i=-nf,nf
      do j=-nf,nf
      if ((i == 0) .and. (j == 0)) then
      ccI2(i,j,:)=I2perparr(gg,:)
      elseif ((i == 0) .and. (j  /=  0)) then
      ccI2(i,j,:)=I2perparr(gq,:)
      elseif ((i  /=  0) .and. (j == 0)) then
      ccI2(i,j,:)=I2perparr(qg,:)
      elseif (i  ==  j) then
      ccI2(i,j,:)=I2perparr(qq,:)
      elseif (i  ==  -j) then
      ccI2(i,j,:)=I2perparr(qbq,:)
      else
!     exploits the fact that qbpq=qpq at this order
      ccI2(i,j,:)=I2perparr(qpq,:)
      endif
      enddo
      enddo

c--------------------------------------------------------
c----- Calculation of quark and gluon beam functions ----
c--------------------------------------------------------

c JC: 12/22 -- extended this to i=0 to try to account for gluon case in the same way

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


