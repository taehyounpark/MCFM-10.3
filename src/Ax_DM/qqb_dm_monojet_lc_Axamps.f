!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_dm_monojet_lc_Axamps(p,i1,i2,i3,i4,i5,amp)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'dm_params.f'
      include 'zprods_decl.f'
      include 'scale.f'
      include 'epinv.f'
c----- LEADING COLOR
c----- fills amplitude for q g qb chi,chib
      complex(dp):: amp(2,2,2,2),Vamp(2,2,2,2),Famp(2,2,2,2)
      complex(dp):: amp_tree(2,2,2,2)
      real(dp):: p(mxpart,4),q(mxpart,4)
      integer:: i1,i2,i3,i4,i5
      integer:: h1,h2,h3,h4
      complex(dp):: L0,Lsm1,L1
      complex(dp):: lnrat
      real(dp):: s(mxpart,mxpart)
      complex(dp):: vfac(2,2)
      real(dp):: bp,beta,s34



      amp(:,:,:,:)=czip
      if(xmass>1d-8) then
c--------- generate massless phase space
      call gen_masslessvecs(p,q,i4,i5)
c--------- generate spinors
      call spinoru(5,q,za,zb)

      else
c-------- massless dm can use usual spinoru
         call spinoru(5,p,za,zb)

      endif
 !     s(:,:)=czip
 !     do j=1,5
 !        do k=1,5
 !           s(j,k)=za(j,k)*zb(k,j)
 !        enddo
 !     enddo

      call dotem(5,p,s)
cbp = 1/2(1+beta)
c beta = sqrt(1-4xmass**2/s34)
      s34=real(za(i4,i5)*zb(i5,i4))
      beta=sqrt(1d0-4d0*xmass**2/s34)
      bp=0.5d0*(one+beta)


c------ basis integrals
c      l12=lnrat(musq,-s(i1,i2))
c      l23=lnrat(musq,-s(i2,i3))

c----- use same notation as Born
c----- split into Vamp and Famp with amp=Vamp+Famp

c------ Vamp piece universial
c      vfac=
c     & -(epinv**2+epinv*l12+0.5d0*l12**2)
c     & -(epinv**2+epinv*l23+0.5d0*l23**2)
c     & -2d0*(epinv+l23)-4d0
c     & +0.5d0*(epinv+l23)+1d0

      vfac(2,2)=  -3d0 - 2*epinv**2 - epinv*lnrat(musq,-s(i1,i2)) -
     &  0.5*lnrat(musq,-s(i1,i2))**2 -
     &  epinv*lnrat(musq,-s(i2,i3)) -
     &  0.5*lnrat(musq,-s(i2,i3))**2 -
     &  1.5*(epinv + lnrat(musq,-s(i2,i3)))
      vfac(2,1)= -3d0 - 2*epinv**2 - epinv*lnrat(musq,-s(i1,i2)) -
     &  0.5*lnrat(musq,-s(i1,i2))**2 -
     &  1.5*(epinv + lnrat(musq,-s(i1,i2))) -
     &  epinv*lnrat(musq,-s(i2,i3)) -
     &  0.5*lnrat(musq,-s(i2,i3))**2
      vfac(1,2)=-3d0 - 2*epinv**2 - epinv*lnrat(musq,-s(i1,i2)) -
     &  0.5*lnrat(musq,-s(i1,i2))**2 -
     &  1.5*(epinv + lnrat(musq,-s(i1,i2))) -
     &  epinv*lnrat(musq,-s(i2,i3)) -
     &     0.5*lnrat(musq,-s(i2,i3))**2
      vfac(1,1)=-3d0 - 2*epinv**2 - epinv*lnrat(musq,-s(i1,i2)) -
     &  0.5*lnrat(musq,-s(i1,i2))**2 -
     &  epinv*lnrat(musq,-s(i2,i3)) -
     &  0.5*lnrat(musq,-s(i2,i3))**2 -
     &  1.5*(epinv + lnrat(musq,-s(i2,i3)))

      call qqb_dm_monojet_Axamps(p,i1,i2,i3,i4,i5,amp_tree)
      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2

                  Vamp(h1,h2,h3,h4)=amp_tree(h1,h2,h3,h4)*vfac(h1,h2)
               enddo
            enddo
         enddo
      enddo

c      write(6,*) 'amp_tree ',amp_tree(2,2,1,2)/s(i4,i5)
c      call A51(i1,i2,i3,i4,i5,za,zb)
c      write(6,*) ' my imp ', Vamp(2,2,1,2)/s(i4,i5)



c      pause
c------ F pieces ==

c===== helicity conserving amplitudes
      Famp(1,1,2,1)= (L1(-s(i2,i3),-s(i4,i5))*za(i1,i2)*za(i1,i5)*
     &     zb(i3,i1)**2*
     &     zb(i4,i2))/(2.*s(i4,i5)*zb(i2,i1)*zb(i3,i2)) -
     &  (bp*L1(-s(i2,i3),-s(i4,i5))*za(i1,i2)*za(i1,i5)*
     &     zb(i3,i1)**2*zb(i4,i2))/(s(i4,i5)*zb(i2,i1)*zb(i3,i2))
     &   + (Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*
     &     za(i2,i5)*zb(i4,i3))/zb(i2,i1) -
     &  (2*bp*Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*
     &     za(i2,i5)*zb(i4,i3))/zb(i2,i1) +
     &  (L0(-s(i2,i3),-s(i4,i5))*za(i1,i5)*zb(i3,i1)*zb(i4,i3))/
     &   (zb(i2,i1)*zb(i3,i2)) -
     &  (2*bp*L0(-s(i2,i3),-s(i4,i5))*za(i1,i5)*zb(i3,i1)*
     &     zb(i4,i3))/(zb(i2,i1)*zb(i3,i2)) +
     &  (Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i1,i5)*
     &     zb(i3,i1)*zb(i4,i3))/(zb(i2,i1)*zb(i3,i2)) -
     &  (2*bp*Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*
     &     za(i1,i5)*zb(i3,i1)*zb(i4,i3))/(zb(i2,i1)*zb(i3,i2)) +
     &  (L1(-s(i2,i3),-s(i4,i5))*za(i1,i3)*za(i1,i5)*zb(i3,i1)**2*
     &     zb(i4,i3))/(2.*s(i4,i5)*zb(i2,i1)*zb(i3,i2)) -
     &  (bp*L1(-s(i2,i3),-s(i4,i5))*za(i1,i3)*za(i1,i5)*
     &     zb(i3,i1)**2*zb(i4,i3))/(s(i4,i5)*zb(i2,i1)*zb(i3,i2))
      Famp(1,2,2,1)= (Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*
     &     za(i1,i5)*
     &     zb(i4,i2))/za(i2,i3) -
     &  (2*bp*Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*za(i1,i5)*
     &     zb(i4,i2))/za(i2,i3) +
     &  (L0(-s(i1,i2),-s(i4,i5))*za(i1,i3)*za(i1,i5)*zb(i4,i3))/
     &   (za(i1,i2)*za(i2,i3)) -
     &  (2*bp*L0(-s(i1,i2),-s(i4,i5))*za(i1,i3)*za(i1,i5)*zb(i4,i3))/
     &   (za(i1,i2)*za(i2,i3)) +
     &  (Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*za(i1,i3)*
     &     za(i1,i5)*zb(i4,i3))/(za(i1,i2)*za(i2,i3)) -
     &  (2*bp*Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*za(i1,i3)*
     &     za(i1,i5)*zb(i4,i3))/(za(i1,i2)*za(i2,i3)) +
     &  (L1(-s(i1,i2),-s(i4,i5))*za(i1,i3)**2*za(i1,i5)*zb(i3,i1)*
     &     zb(i4,i3))/(2.*s(i4,i5)*za(i1,i2)*za(i2,i3)) -
     &  (bp*L1(-s(i1,i2),-s(i4,i5))*za(i1,i3)**2*za(i1,i5)*zb(i3,i1)*
     &     zb(i4,i3))/(s(i4,i5)*za(i1,i2)*za(i2,i3)) +
     &  (L1(-s(i1,i2),-s(i4,i5))*za(i1,i3)**2*za(i2,i5)*zb(i3,i2)*
     &     zb(i4,i3))/(2.*s(i4,i5)*za(i1,i2)*za(i2,i3)) -
     &  (bp*L1(-s(i1,i2),-s(i4,i5))*za(i1,i3)**2*za(i2,i5)*zb(i3,i2)*
     &     zb(i4,i3))/(s(i4,i5)*za(i1,i2)*za(i2,i3))
      Famp(2,1,2,1)=-((Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*
     &     za(i2,i5)*
     &       zb(i4,i1))/zb(i3,i2)) +
     &  (2*bp*Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*za(i2,i5)*
     &     zb(i4,i1))/zb(i3,i2) -
     &  (L0(-s(i1,i2),-s(i4,i5))*za(i3,i5)*zb(i3,i1)*zb(i4,i1))/
     &   (zb(i2,i1)*zb(i3,i2)) +
     &  (2*bp*L0(-s(i1,i2),-s(i4,i5))*za(i3,i5)*zb(i3,i1)*zb(i4,i1))/
     &   (zb(i2,i1)*zb(i3,i2)) -
     &  (Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*za(i3,i5)*
     &     zb(i3,i1)*zb(i4,i1))/(zb(i2,i1)*zb(i3,i2)) +
     &  (2*bp*Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*za(i3,i5)*
     &     zb(i3,i1)*zb(i4,i1))/(zb(i2,i1)*zb(i3,i2)) -
     &  (L1(-s(i1,i2),-s(i4,i5))*za(i1,i3)*za(i3,i5)*zb(i3,i1)**2*
     &     zb(i4,i1))/(2.*s(i4,i5)*zb(i2,i1)*zb(i3,i2)) +
     &  (bp*L1(-s(i1,i2),-s(i4,i5))*za(i1,i3)*za(i3,i5)*zb(i3,i1)**2*
     &     zb(i4,i1))/(s(i4,i5)*zb(i2,i1)*zb(i3,i2)) -
     &  (L1(-s(i1,i2),-s(i4,i5))*za(i2,i3)*za(i3,i5)*zb(i3,i1)**2*
     &     zb(i4,i2))/(2.*s(i4,i5)*zb(i2,i1)*zb(i3,i2)) +
     &  (bp*L1(-s(i1,i2),-s(i4,i5))*za(i2,i3)*za(i3,i5)*zb(i3,i1)**2*
     &     zb(i4,i2))/(s(i4,i5)*zb(i2,i1)*zb(i3,i2))
      Famp(2,2,2,1)= -((L0(-s(i2,i3),-s(i4,i5))*za(i1,i3)*za(i3,i5)*
     &     zb(i4,i1))/
     &     (za(i1,i2)*za(i2,i3))) +
     &  (2*bp*L0(-s(i2,i3),-s(i4,i5))*za(i1,i3)*za(i3,i5)*zb(i4,i1))/
     &   (za(i1,i2)*za(i2,i3)) -
     &  (Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i1,i3)*
     &     za(i3,i5)*zb(i4,i1))/(za(i1,i2)*za(i2,i3)) +
     &  (2*bp*Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i1,i3)*
     &     za(i3,i5)*zb(i4,i1))/(za(i1,i2)*za(i2,i3)) -
     &  (L1(-s(i2,i3),-s(i4,i5))*za(i1,i3)**2*za(i2,i5)*zb(i2,i1)*
     &     zb(i4,i1))/(2.*s(i4,i5)*za(i1,i2)*za(i2,i3)) +
     &  (bp*L1(-s(i2,i3),-s(i4,i5))*za(i1,i3)**2*za(i2,i5)*zb(i2,i1)*
     &     zb(i4,i1))/(s(i4,i5)*za(i1,i2)*za(i2,i3)) -
     &  (L1(-s(i2,i3),-s(i4,i5))*za(i1,i3)**2*za(i3,i5)*zb(i3,i1)*
     &     zb(i4,i1))/(2.*s(i4,i5)*za(i1,i2)*za(i2,i3)) +
     &  (bp*L1(-s(i2,i3),-s(i4,i5))*za(i1,i3)**2*za(i3,i5)*zb(i3,i1)*
     &     zb(i4,i1))/(s(i4,i5)*za(i1,i2)*za(i2,i3)) -
     &  (Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i3,i5)*
     &     zb(i4,i2))/za(i1,i2) +
     &  (2*bp*Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i3,i5)*
     &     zb(i4,i2))/za(i1,i2)

      Famp(1,1,1,2)=-(L1(-s(i2,i3),-s(i4,i5))*za(i1,i2)*za(i1,i4)*
     &     zb(i3,i1)**2*
     &      zb(i5,i2))/(2.*s(i4,i5)*zb(i2,i1)*zb(i3,i2)) +
     &  (bp*L1(-s(i2,i3),-s(i4,i5))*za(i1,i2)*za(i1,i4)*zb(i3,i1)**2*
     &     zb(i5,i2))/(s(i4,i5)*zb(i2,i1)*zb(i3,i2)) -
     &  (Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i2,i4)*
     &     zb(i5,i3))/zb(i2,i1) +
     &  (2*bp*Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i2,i4)*
     &     zb(i5,i3))/zb(i2,i1) -
     &  (L0(-s(i2,i3),-s(i4,i5))*za(i1,i4)*zb(i3,i1)*zb(i5,i3))/
     &   (zb(i2,i1)*zb(i3,i2)) +
     &  (2*bp*L0(-s(i2,i3),-s(i4,i5))*za(i1,i4)*zb(i3,i1)*zb(i5,i3))/
     &   (zb(i2,i1)*zb(i3,i2)) -
     &  (Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i1,i4)*
     &     zb(i3,i1)*zb(i5,i3))/(zb(i2,i1)*zb(i3,i2)) +
     &  (2*bp*Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i1,i4)*
     &     zb(i3,i1)*zb(i5,i3))/(zb(i2,i1)*zb(i3,i2)) -
     &  (L1(-s(i2,i3),-s(i4,i5))*za(i1,i3)*za(i1,i4)*zb(i3,i1)**2*
     &     zb(i5,i3))/(2.*s(i4,i5)*zb(i2,i1)*zb(i3,i2)) +
     &  (bp*L1(-s(i2,i3),-s(i4,i5))*za(i1,i3)*za(i1,i4)*zb(i3,i1)**2*
     &     zb(i5,i3))/(s(i4,i5)*zb(i2,i1)*zb(i3,i2))
      Famp(1,2,1,2)= -((Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*
     &     za(i1,i4)*
     &       zb(i5,i2))/za(i2,i3)) +
     &  (2*bp*Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*za(i1,i4)*
     &     zb(i5,i2))/za(i2,i3) -
     &  (L0(-s(i1,i2),-s(i4,i5))*za(i1,i3)*za(i1,i4)*zb(i5,i3))/
     &   (za(i1,i2)*za(i2,i3)) +
     &  (2*bp*L0(-s(i1,i2),-s(i4,i5))*za(i1,i3)*za(i1,i4)*zb(i5,i3))/
     &   (za(i1,i2)*za(i2,i3)) -
     &  (Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*za(i1,i3)*
     &     za(i1,i4)*zb(i5,i3))/(za(i1,i2)*za(i2,i3)) +
     &  (2*bp*Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*za(i1,i3)*
     &     za(i1,i4)*zb(i5,i3))/(za(i1,i2)*za(i2,i3)) -
     &  (L1(-s(i1,i2),-s(i4,i5))*za(i1,i3)**2*za(i1,i4)*zb(i3,i1)*
     &     zb(i5,i3))/(2.*s(i4,i5)*za(i1,i2)*za(i2,i3)) +
     &  (bp*L1(-s(i1,i2),-s(i4,i5))*za(i1,i3)**2*za(i1,i4)*zb(i3,i1)*
     &     zb(i5,i3))/(s(i4,i5)*za(i1,i2)*za(i2,i3)) -
     &  (L1(-s(i1,i2),-s(i4,i5))*za(i1,i3)**2*za(i2,i4)*zb(i3,i2)*
     &     zb(i5,i3))/(2.*s(i4,i5)*za(i1,i2)*za(i2,i3)) +
     &  (bp*L1(-s(i1,i2),-s(i4,i5))*za(i1,i3)**2*za(i2,i4)*zb(i3,i2)*
     &     zb(i5,i3))/(s(i4,i5)*za(i1,i2)*za(i2,i3))
      Famp(2,1,1,2)= (Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))
     &     *za(i2,i4)*
     &     zb(i5,i1))/zb(i3,i2) -
     &  (2*bp*Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*za(i2,i4)*
     &     zb(i5,i1))/zb(i3,i2) +
     &  (L0(-s(i1,i2),-s(i4,i5))*za(i3,i4)*zb(i3,i1)*zb(i5,i1))/
     &   (zb(i2,i1)*zb(i3,i2)) -
     &  (2*bp*L0(-s(i1,i2),-s(i4,i5))*za(i3,i4)*zb(i3,i1)*zb(i5,i1))/
     &   (zb(i2,i1)*zb(i3,i2)) +
     &  (Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*za(i3,i4)*
     &     zb(i3,i1)*zb(i5,i1))/(zb(i2,i1)*zb(i3,i2)) -
     &  (2*bp*Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*za(i3,i4)*
     &     zb(i3,i1)*zb(i5,i1))/(zb(i2,i1)*zb(i3,i2)) +
     &  (L1(-s(i1,i2),-s(i4,i5))*za(i1,i3)*za(i3,i4)*zb(i3,i1)**2*
     &     zb(i5,i1))/(2.*s(i4,i5)*zb(i2,i1)*zb(i3,i2)) -
     &  (bp*L1(-s(i1,i2),-s(i4,i5))*za(i1,i3)*za(i3,i4)*zb(i3,i1)**2*
     &     zb(i5,i1))/(s(i4,i5)*zb(i2,i1)*zb(i3,i2)) +
     &  (L1(-s(i1,i2),-s(i4,i5))*za(i2,i3)*za(i3,i4)*zb(i3,i1)**2*
     &     zb(i5,i2))/(2.*s(i4,i5)*zb(i2,i1)*zb(i3,i2)) -
     &  (bp*L1(-s(i1,i2),-s(i4,i5))*za(i2,i3)*za(i3,i4)*zb(i3,i1)**2*
     &     zb(i5,i2))/(s(i4,i5)*zb(i2,i1)*zb(i3,i2))
      Famp(2,2,1,2)=(L0(-s(i2,i3),-s(i4,i5))*za(i1,i3)*za(i3,i4)*
     &     zb(i5,i1))/
     &   (za(i1,i2)*za(i2,i3)) -
     &  (2*bp*L0(-s(i2,i3),-s(i4,i5))*za(i1,i3)*za(i3,i4)*zb(i5,i1))/
     &   (za(i1,i2)*za(i2,i3)) +
     &  (Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i1,i3)*
     &     za(i3,i4)*zb(i5,i1))/(za(i1,i2)*za(i2,i3)) -
     &  (2*bp*Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i1,i3)*
     &     za(i3,i4)*zb(i5,i1))/(za(i1,i2)*za(i2,i3)) +
     &  (L1(-s(i2,i3),-s(i4,i5))*za(i1,i3)**2*za(i2,i4)*zb(i2,i1)*
     &     zb(i5,i1))/(2.*s(i4,i5)*za(i1,i2)*za(i2,i3)) -
     &  (bp*L1(-s(i2,i3),-s(i4,i5))*za(i1,i3)**2*za(i2,i4)*zb(i2,i1)*
     &     zb(i5,i1))/(s(i4,i5)*za(i1,i2)*za(i2,i3)) +
     &  (L1(-s(i2,i3),-s(i4,i5))*za(i1,i3)**2*za(i3,i4)*zb(i3,i1)*
     &     zb(i5,i1))/(2.*s(i4,i5)*za(i1,i2)*za(i2,i3)) -
     &  (bp*L1(-s(i2,i3),-s(i4,i5))*za(i1,i3)**2*za(i3,i4)*zb(i3,i1)*
     &     zb(i5,i1))/(s(i4,i5)*za(i1,i2)*za(i2,i3)) +
     &  (Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i3,i4)*
     &     zb(i5,i2))/za(i1,i2) -
     &  (2*bp*Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i3,i4)*
     &     zb(i5,i2))/za(i1,i2)

c-------------- helicity violating

      Famp(1,1,1,1)= (xmass*L1(-s(i2,i3),-s(i4,i5))*za(i1,i2)*za(i1,i4)*
     &     zb(i3,i1)**2*
     &     zb(i4,i2))/(2.*s(i4,i5)*zb(i2,i1)*zb(i3,i2)*zb(i5,i4)) +
     &  (xmass*Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i2,i4)*
     &     zb(i4,i3))/(zb(i2,i1)*zb(i5,i4)) +
     &  (xmass*L0(-s(i2,i3),-s(i4,i5))*za(i1,i4)*zb(i3,i1)*zb(i4,i3))/
     &   (zb(i2,i1)*zb(i3,i2)*zb(i5,i4)) +
     &  (xmass*Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i1,i4)*
     &     zb(i3,i1)*zb(i4,i3))/(zb(i2,i1)*zb(i3,i2)*zb(i5,i4)) +
     &  (xmass*L1(-s(i2,i3),-s(i4,i5))*za(i1,i3)*za(i1,i4)*zb(i3,i1)**2*
     &     zb(i4,i3))/(2.*s(i4,i5)*zb(i2,i1)*zb(i3,i2)*zb(i5,i4)) +
     &  (xmass*L1(-s(i2,i3),-s(i4,i5))*za(i1,i2)*za(i1,i5)*zb(i3,i1)**2*
     &     zb(i5,i2))/(2.*s(i4,i5)*zb(i2,i1)*zb(i3,i2)*zb(i5,i4)) +
     &  (xmass*Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i2,i5)*
     &     zb(i5,i3))/(zb(i2,i1)*zb(i5,i4)) +
     &  (xmass*L0(-s(i2,i3),-s(i4,i5))*za(i1,i5)*zb(i3,i1)*zb(i5,i3))/
     &   (zb(i2,i1)*zb(i3,i2)*zb(i5,i4)) +
     &  (xmass*Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i1,i5)*
     &     zb(i3,i1)*zb(i5,i3))/(zb(i2,i1)*zb(i3,i2)*zb(i5,i4)) +
     &  (xmass*L1(-s(i2,i3),-s(i4,i5))*za(i1,i3)*za(i1,i5)*zb(i3,i1)**2*
     &     zb(i5,i3))/(2.*s(i4,i5)*zb(i2,i1)*zb(i3,i2)*zb(i5,i4))
      Famp(1,2,1,1)=(xmass*Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))
     & *za(i1,i4)*
     &     zb(i4,i2))/(za(i2,i3)*zb(i5,i4)) +
     &  (xmass*L0(-s(i1,i2),-s(i4,i5))*za(i1,i3)*za(i1,i4)*zb(i4,i3))/
     &   (za(i1,i2)*za(i2,i3)*zb(i5,i4)) +
     &  (xmass*Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*za(i1,i3)*
     &     za(i1,i4)*zb(i4,i3))/(za(i1,i2)*za(i2,i3)*zb(i5,i4)) +
     &  (xmass*L1(-s(i1,i2),-s(i4,i5))*za(i1,i3)**2*za(i1,i4)*zb(i3,i1)*
     &     zb(i4,i3))/(2.*s(i4,i5)*za(i1,i2)*za(i2,i3)*zb(i5,i4)) +
     &  (xmass*L1(-s(i1,i2),-s(i4,i5))*za(i1,i3)**2*za(i2,i4)*zb(i3,i2)*
     &     zb(i4,i3))/(2.*s(i4,i5)*za(i1,i2)*za(i2,i3)*zb(i5,i4)) +
     &  (xmass*Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*za(i1,i5)*
     &     zb(i5,i2))/(za(i2,i3)*zb(i5,i4)) +
     &  (xmass*L0(-s(i1,i2),-s(i4,i5))*za(i1,i3)*za(i1,i5)*zb(i5,i3))/
     &   (za(i1,i2)*za(i2,i3)*zb(i5,i4)) +
     &  (xmass*Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*za(i1,i3)*
     &     za(i1,i5)*zb(i5,i3))/(za(i1,i2)*za(i2,i3)*zb(i5,i4)) +
     &  (xmass*L1(-s(i1,i2),-s(i4,i5))*za(i1,i3)**2*za(i1,i5)*zb(i3,i1)*
     &     zb(i5,i3))/(2.*s(i4,i5)*za(i1,i2)*za(i2,i3)*zb(i5,i4)) +
     &  (xmass*L1(-s(i1,i2),-s(i4,i5))*za(i1,i3)**2*za(i2,i5)*zb(i3,i2)*
     &     zb(i5,i3))/(2.*s(i4,i5)*za(i1,i2)*za(i2,i3)*zb(i5,i4))
      Famp(2,1,1,1)=-((xmass*
     &     Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*za(i2,i4)*
     &       zb(i4,i1))/(zb(i3,i2)*zb(i5,i4))) -
     &  (xmass*L0(-s(i1,i2),-s(i4,i5))*za(i3,i4)*zb(i3,i1)*zb(i4,i1))/
     &   (zb(i2,i1)*zb(i3,i2)*zb(i5,i4)) -
     &  (xmass*Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*za(i3,i4)*
     &     zb(i3,i1)*zb(i4,i1))/(zb(i2,i1)*zb(i3,i2)*zb(i5,i4)) -
     &  (xmass*L1(-s(i1,i2),-s(i4,i5))*za(i1,i3)*za(i3,i4)*zb(i3,i1)**2*
     &     zb(i4,i1))/(2.*s(i4,i5)*zb(i2,i1)*zb(i3,i2)*zb(i5,i4)) -
     &  (xmass*L1(-s(i1,i2),-s(i4,i5))*za(i2,i3)*za(i3,i4)*zb(i3,i1)**2*
     &     zb(i4,i2))/(2.*s(i4,i5)*zb(i2,i1)*zb(i3,i2)*zb(i5,i4)) -
     &  (xmass*Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*za(i2,i5)*
     &     zb(i5,i1))/(zb(i3,i2)*zb(i5,i4)) -
     &  (xmass*L0(-s(i1,i2),-s(i4,i5))*za(i3,i5)*zb(i3,i1)*zb(i5,i1))/
     &   (zb(i2,i1)*zb(i3,i2)*zb(i5,i4)) -
     &  (xmass*Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*za(i3,i5)*
     &     zb(i3,i1)*zb(i5,i1))/(zb(i2,i1)*zb(i3,i2)*zb(i5,i4)) -
     &  (xmass*L1(-s(i1,i2),-s(i4,i5))*za(i1,i3)*za(i3,i5)*zb(i3,i1)**2*
     &     zb(i5,i1))/(2.*s(i4,i5)*zb(i2,i1)*zb(i3,i2)*zb(i5,i4)) -
     &  (xmass*L1(-s(i1,i2),-s(i4,i5))*za(i2,i3)*za(i3,i5)*zb(i3,i1)**2*
     &     zb(i5,i2))/(2.*s(i4,i5)*zb(i2,i1)*zb(i3,i2)*zb(i5,i4))
      Famp(2,2,1,1)=-((xmass*L0(-s(i2,i3),-s(i4,i5))*za(i1,i3)*za(i3,i4)
     &     *zb(i4,i1))/
     &     (za(i1,i2)*za(i2,i3)*zb(i5,i4))) -
     &  (xmass*Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i1,i3)*
     &     za(i3,i4)*zb(i4,i1))/(za(i1,i2)*za(i2,i3)*zb(i5,i4)) -
     &  (xmass*L1(-s(i2,i3),-s(i4,i5))*za(i1,i3)**2*za(i2,i4)*zb(i2,i1)*
     &     zb(i4,i1))/(2.*s(i4,i5)*za(i1,i2)*za(i2,i3)*zb(i5,i4)) -
     &  (xmass*L1(-s(i2,i3),-s(i4,i5))*za(i1,i3)**2*za(i3,i4)*zb(i3,i1)*
     &     zb(i4,i1))/(2.*s(i4,i5)*za(i1,i2)*za(i2,i3)*zb(i5,i4)) -
     &  (xmass*Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i3,i4)*
     &     zb(i4,i2))/(za(i1,i2)*zb(i5,i4)) -
     &  (xmass*L0(-s(i2,i3),-s(i4,i5))*za(i1,i3)*za(i3,i5)*zb(i5,i1))/
     &   (za(i1,i2)*za(i2,i3)*zb(i5,i4)) -
     &  (xmass*Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i1,i3)*
     &     za(i3,i5)*zb(i5,i1))/(za(i1,i2)*za(i2,i3)*zb(i5,i4)) -
     &  (xmass*L1(-s(i2,i3),-s(i4,i5))*za(i1,i3)**2*za(i2,i5)*zb(i2,i1)*
     &     zb(i5,i1))/(2.*s(i4,i5)*za(i1,i2)*za(i2,i3)*zb(i5,i4)) -
     &  (xmass*L1(-s(i2,i3),-s(i4,i5))*za(i1,i3)**2*za(i3,i5)*zb(i3,i1)*
     &     zb(i5,i1))/(2.*s(i4,i5)*za(i1,i2)*za(i2,i3)*zb(i5,i4)) -
     &  (xmass*Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i3,i5)*
     &     zb(i5,i2))/(za(i1,i2)*zb(i5,i4))

      Famp(1,1,2,2)=(xmass*L1(-s(i2,i3),-s(i4,i5))*za(i1,i2)*za(i1,i4)*
     &     zb(i3,i1)**2*
     &     zb(i4,i2))/(2.*s(i4,i5)*za(i4,i5)*zb(i2,i1)*zb(i3,i2)) +
     &  (xmass*Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i2,i4)*
     &     zb(i4,i3))/(za(i4,i5)*zb(i2,i1)) +
     &  (xmass*L0(-s(i2,i3),-s(i4,i5))*za(i1,i4)*zb(i3,i1)*zb(i4,i3))/
     &   (za(i4,i5)*zb(i2,i1)*zb(i3,i2)) +
     &  (xmass*Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i1,i4)*
     &     zb(i3,i1)*zb(i4,i3))/(za(i4,i5)*zb(i2,i1)*zb(i3,i2)) +
     &  (xmass*L1(-s(i2,i3),-s(i4,i5))*za(i1,i3)*za(i1,i4)*zb(i3,i1)**2*
     &     zb(i4,i3))/(2.*s(i4,i5)*za(i4,i5)*zb(i2,i1)*zb(i3,i2)) +
     &  (xmass*L1(-s(i2,i3),-s(i4,i5))*za(i1,i2)*za(i1,i5)*zb(i3,i1)**2*
     &     zb(i5,i2))/(2.*s(i4,i5)*za(i4,i5)*zb(i2,i1)*zb(i3,i2)) +
     &  (xmass*Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i2,i5)*
     &     zb(i5,i3))/(za(i4,i5)*zb(i2,i1)) +
     &  (xmass*L0(-s(i2,i3),-s(i4,i5))*za(i1,i5)*zb(i3,i1)*zb(i5,i3))/
     &   (za(i4,i5)*zb(i2,i1)*zb(i3,i2)) +
     &  (xmass*Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i1,i5)*
     &     zb(i3,i1)*zb(i5,i3))/(za(i4,i5)*zb(i2,i1)*zb(i3,i2)) +
     &  (xmass*L1(-s(i2,i3),-s(i4,i5))*za(i1,i3)*za(i1,i5)*zb(i3,i1)**2*
     &     zb(i5,i3))/(2.*s(i4,i5)*za(i4,i5)*zb(i2,i1)*zb(i3,i2))
      Famp(1,2,2,2)=(xmass*Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))
     &     * za(i1,i4)*
     &     zb(i4,i2))/(za(i2,i3)*za(i4,i5)) +
     &  (xmass*L0(-s(i1,i2),-s(i4,i5))*za(i1,i3)*za(i1,i4)*zb(i4,i3))/
     &   (za(i1,i2)*za(i2,i3)*za(i4,i5)) +
     &  (xmass*Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*za(i1,i3)*
     &     za(i1,i4)*zb(i4,i3))/(za(i1,i2)*za(i2,i3)*za(i4,i5)) +
     &  (xmass*L1(-s(i1,i2),-s(i4,i5))*za(i1,i3)**2*za(i1,i4)*zb(i3,i1)*
     &     zb(i4,i3))/(2.*s(i4,i5)*za(i1,i2)*za(i2,i3)*za(i4,i5)) +
     &  (xmass*L1(-s(i1,i2),-s(i4,i5))*za(i1,i3)**2*za(i2,i4)*zb(i3,i2)*
     &     zb(i4,i3))/(2.*s(i4,i5)*za(i1,i2)*za(i2,i3)*za(i4,i5)) +
     &  (xmass*Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*za(i1,i5)*
     &     zb(i5,i2))/(za(i2,i3)*za(i4,i5)) +
     &  (xmass*L0(-s(i1,i2),-s(i4,i5))*za(i1,i3)*za(i1,i5)*zb(i5,i3))/
     &   (za(i1,i2)*za(i2,i3)*za(i4,i5)) +
     &  (xmass*Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*za(i1,i3)*
     &     za(i1,i5)*zb(i5,i3))/(za(i1,i2)*za(i2,i3)*za(i4,i5)) +
     &  (xmass*L1(-s(i1,i2),-s(i4,i5))*za(i1,i3)**2*za(i1,i5)*zb(i3,i1)*
     &     zb(i5,i3))/(2.*s(i4,i5)*za(i1,i2)*za(i2,i3)*za(i4,i5)) +
     &  (xmass*L1(-s(i1,i2),-s(i4,i5))*za(i1,i3)**2*za(i2,i5)*zb(i3,i2)*
     &     zb(i5,i3))/(2.*s(i4,i5)*za(i1,i2)*za(i2,i3)*za(i4,i5))
      Famp(2,1,2,2)= -((xmass*
     &     Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*za(i2,i4)*
     &       zb(i4,i1))/(za(i4,i5)*zb(i3,i2))) -
     &  (xmass*L0(-s(i1,i2),-s(i4,i5))*za(i3,i4)*zb(i3,i1)*zb(i4,i1))/
     &   (za(i4,i5)*zb(i2,i1)*zb(i3,i2)) -
     &  (xmass*Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*za(i3,i4)*
     &     zb(i3,i1)*zb(i4,i1))/(za(i4,i5)*zb(i2,i1)*zb(i3,i2)) -
     &  (xmass*L1(-s(i1,i2),-s(i4,i5))*za(i1,i3)*za(i3,i4)*zb(i3,i1)**2*
     &     zb(i4,i1))/(2.*s(i4,i5)*za(i4,i5)*zb(i2,i1)*zb(i3,i2)) -
     &  (xmass*L1(-s(i1,i2),-s(i4,i5))*za(i2,i3)*za(i3,i4)*zb(i3,i1)**2*
     &     zb(i4,i2))/(2.*s(i4,i5)*za(i4,i5)*zb(i2,i1)*zb(i3,i2)) -
     &  (xmass*Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*za(i2,i5)*
     &     zb(i5,i1))/(za(i4,i5)*zb(i3,i2)) -
     &  (xmass*L0(-s(i1,i2),-s(i4,i5))*za(i3,i5)*zb(i3,i1)*zb(i5,i1))/
     &   (za(i4,i5)*zb(i2,i1)*zb(i3,i2)) -
     &  (xmass*Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))*za(i3,i5)*
     &     zb(i3,i1)*zb(i5,i1))/(za(i4,i5)*zb(i2,i1)*zb(i3,i2)) -
     &  (xmass*L1(-s(i1,i2),-s(i4,i5))*za(i1,i3)*za(i3,i5)*zb(i3,i1)**2*
     &     zb(i5,i1))/(2.*s(i4,i5)*za(i4,i5)*zb(i2,i1)*zb(i3,i2)) -
     &  (xmass*L1(-s(i1,i2),-s(i4,i5))*za(i2,i3)*za(i3,i5)*zb(i3,i1)**2*
     &     zb(i5,i2))/(2.*s(i4,i5)*za(i4,i5)*zb(i2,i1)*zb(i3,i2))
      Famp(2,2,2,2)=-((xmass*L0(-s(i2,i3),-s(i4,i5))*za(i1,i3)*za(i3,i4)
     &     *zb(i4,i1))/
     &     (za(i1,i2)*za(i2,i3)*za(i4,i5))) -
     &  (xmass*Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i1,i3)*
     &     za(i3,i4)*zb(i4,i1))/(za(i1,i2)*za(i2,i3)*za(i4,i5)) -
     &  (xmass*L1(-s(i2,i3),-s(i4,i5))*za(i1,i3)**2*za(i2,i4)*zb(i2,i1)*
     &     zb(i4,i1))/(2.*s(i4,i5)*za(i1,i2)*za(i2,i3)*za(i4,i5)) -
     &  (xmass*L1(-s(i2,i3),-s(i4,i5))*za(i1,i3)**2*za(i3,i4)*zb(i3,i1)*
     &     zb(i4,i1))/(2.*s(i4,i5)*za(i1,i2)*za(i2,i3)*za(i4,i5)) -
     &  (xmass*Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i3,i4)*
     &     zb(i4,i2))/(za(i1,i2)*za(i4,i5)) -
     &  (xmass*L0(-s(i2,i3),-s(i4,i5))*za(i1,i3)*za(i3,i5)*zb(i5,i1))/
     &   (za(i1,i2)*za(i2,i3)*za(i4,i5)) -
     &  (xmass*Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i1,i3)*
     &     za(i3,i5)*zb(i5,i1))/(za(i1,i2)*za(i2,i3)*za(i4,i5)) -
     &  (xmass*L1(-s(i2,i3),-s(i4,i5))*za(i1,i3)**2*za(i2,i5)*zb(i2,i1)*
     &     zb(i5,i1))/(2.*s(i4,i5)*za(i1,i2)*za(i2,i3)*za(i4,i5)) -
     &  (xmass*L1(-s(i2,i3),-s(i4,i5))*za(i1,i3)**2*za(i3,i5)*zb(i3,i1)*
     &     zb(i5,i1))/(2.*s(i4,i5)*za(i1,i2)*za(i2,i3)*za(i4,i5)) -
     &  (xmass*Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))*za(i3,i5)*
     &     zb(i5,i2))/(za(i1,i2)*za(i4,i5))

c----- bulid final amplitude
      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
                  amp(h1,h2,h3,h4)=Vamp(h1,h2,h3,h4)+Famp(h1,h2,h3,h4)
               enddo
            enddo
         enddo
      enddo
c      write(6,*) 'lc Ax amp '
c      do h1=1,2
c         do h2=1,2
c            do h3=1,2
c               do h4=1,2
c             write(6,*) h1,h2,h3,h4,Vamp(h1,h2,h3,h4),Famp(h1,h2,h3,h4)
c          enddo
c       enddo
c      enddo
c      enddo
c      pause

c      call  qqb_dm_monojet_lc_Vamps(p,i1,i2,i3,i4,i5,amp)

c      call A51(i1,i2,i3,i4,i5,za,zb)

c      write(6,*) 'my V*tree ',Vamp(2,2,1,2)/s(i4,i5)
c      write(6,*) 'my F ',Famp(2,2,1,2)/s(i4,i5)
c3      pause

      return
      end
