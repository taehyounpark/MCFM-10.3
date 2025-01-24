!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_w_ew_v(p,msq)
c----Matrix element for W production
c----averaged over initial colours and spins
c For nwz=+1
c     u(-p1)+dbar(-p2)-->W^+(n(p3)+e^+(p4))
c For nwz=-1
c     d(-p1)+ubar(-p2)-->W^-(e^-(p3)+nbar(p4))
c--- with one loop electroweak radiative corrections
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nwz.f'
      include 'masses.f'
      include 'scale.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'epinv.f'
      include 'scheme.f'
!      include 'msq.f'
      include 'zcouple_cms.f'
      include 'ckm.f'
      include 'ischeme.f'
      include 'ewcouple.f'
      include 'scalarselect.f'
      real (dp) :: p(mxpart,4)
      real(dp) :: s,t,u
      real(dp) :: msq(-nf:nf,-nf:nf),ss,delta(-2:0),deltaqqb,deltaqbq,fac,qqb,qbq
      integer:: j,k
      real(dp) :: theta,phi,rho,csig,muk,beta,ssig,pm(mxpart,4),m
      integer nu
c--statement function
      ss(j,k)=two*(p(j,4)*p(k,4)-p(j,1)*p(k,1)
     &            -p(j,2)*p(k,2)-p(j,3)*p(k,3))

      scheme='tH-V'

      fac=gwsq**2*NC

! This process uses scalar integrals with complex values for the masses.
! Since these are not implemented faithfully in the current QCDLoop version,
! this process enforces the use of OneLOop.
      scalarselect=2

!      m=91._dp
!      include 'kinpoint4.f'
!      p(:,:)=pm(:,:)
!      zmass= 91.15348061918277_dp
!      wmass= 80.357973609878_dp
!      zwidth= 2.4942663787728243_dp
!      wwidth= 2.0842989982782_dp
!      mt=173.2_dp
!      twidth=0._dp
!      hmass=125._dp
!      hwidth=0._dp
!      scale=zmass
!      musq=scale**2
!      zaemmz=1d0/137d0
!      zesq=fourpi*zaemmz

!      write(6,*) 'MW',wmass
!      write(6,*) 'MZ',zmass
!      write(6,*) 'WW',wwidth
!      write(6,*) 'WZ',zwidth
      s=+ss(1,2)
      t=+ss(1,3)
      u=+ss(1,4)

      ischeme=3

c---   calculate propagator
      fac=aveqq*fac/((ss(1,2)-wmass**2)**2+(wmass*wwidth)**2)

c---case dbar-u or ubar-d
      call deltaEW(s,t,ischeme,delta)
      deltaqqb=epinv**2*delta(-2)+epinv*delta(-1)+delta(0)
      qqb=fac*ss(1,4)**2*deltaqqb

!      write(6,*)  'ischeme,qqb_w_ew_v:delta(-2)*2*pi/zaemmz',
!     & ischeme,delta(-2)*2*pi/zaemmz
!      write(6,*)  'ischeme,qqb_w_ew_v:delta(-1)*2*pi/zaemmz',
!     & ischeme,delta(-1)*2*pi/zaemmz
!      write(6,*)  'ischeme,qqb_w_ew_v:delta( 0)*2*pi/zaemmz',
!     & ischeme,delta(0)*2*pi/zaemmz

      call deltaEW(s,u,ischeme,delta)
      deltaqbq=epinv**2*delta(-2)+epinv*delta(-1)+delta(0)
      qbq=fac*ss(1,3)**2*deltaqbq

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=zip
          if ((j > 0) .and. (k < 0)) then
            msq(j,k)=Vsq(j,k)*qqb
          elseif ((j < 0) .and. (k > 0)) then
            msq(j,k)=Vsq(j,k)*qbq
          endif
      enddo
      enddo

      return
      end
