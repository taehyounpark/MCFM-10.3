!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c============= helicity amplitude for axial insertion,
      subroutine qqb_dm_monojet_Axamps(p,i1,i2,i3,i4,i5,amp)
      implicit none
      include 'types.f'
      include 'dm_params.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      real(dp):: p(mxpart,4)
c----- fills amplitude for q g qb chi,chib
      complex(dp):: amp(2,2,2,2)
      real(dp):: q(mxpart,4)
      integer:: i1,i2,i3,i4,i5
c------ z1jet amplitude for testing (note order is q qb l1 l2 g)
c------ and gluon helicity is summed over
c------ returns amp squared.
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

cbp = 1/2(1+beta)
c beta = sqrt(1-4xmass**2/s34)
      s34=Dble(za(i4,i5)*zb(i5,i4))
      beta=sqrt(1d0-4d0*xmass**2/s34)
      bp=0.5d0*(one+beta)

c========= Helicity conserving

      amp(1,1,2,1)=-((za(i2,i5)*zb(i4,i3))/zb(i2,i1)) +
     &  (2*bp*za(i2,i5)*zb(i4,i3))/zb(i2,i1) -
     &  (za(i1,i5)*zb(i3,i1)*zb(i4,i3))/
     &   (zb(i2,i1)*zb(i3,i2)) +
     &  (2*bp*za(i1,i5)*zb(i3,i1)*zb(i4,i3))/
     &   (zb(i2,i1)*zb(i3,i2))
      amp(1,2,2,1)=-((za(i1,i5)*zb(i4,i2))/za(i2,i3)) +
     &  (2*bp*za(i1,i5)*zb(i4,i2))/za(i2,i3) -
     &  (za(i1,i3)*za(i1,i5)*zb(i4,i3))/
     &   (za(i1,i2)*za(i2,i3)) +
     &  (2*bp*za(i1,i3)*za(i1,i5)*zb(i4,i3))/
     &   (za(i1,i2)*za(i2,i3))
      amp(2,1,2,1)=(za(i2,i5)*zb(i4,i1))/zb(i3,i2) -
     &  (2*bp*za(i2,i5)*zb(i4,i1))/zb(i3,i2) +
     &  (za(i3,i5)*zb(i3,i1)*zb(i4,i1))/
     &   (zb(i2,i1)*zb(i3,i2)) -
     &  (2*bp*za(i3,i5)*zb(i3,i1)*zb(i4,i1))/
     &   (zb(i2,i1)*zb(i3,i2))
      amp(2,2,2,1)=(za(i1,i3)*za(i3,i5)*zb(i4,i1))/
     &   (za(i1,i2)*za(i2,i3)) -
     &  (2*bp*za(i1,i3)*za(i3,i5)*zb(i4,i1))/
     &   (za(i1,i2)*za(i2,i3)) +
     &  (za(i3,i5)*zb(i4,i2))/za(i1,i2) -
     &  (2*bp*za(i3,i5)*zb(i4,i2))/za(i1,i2)

      amp(1,1,1,2)=(za(i2,i4)*zb(i5,i3))/zb(i2,i1) -
     &  (2*bp*za(i2,i4)*zb(i5,i3))/zb(i2,i1) +
     &  (za(i1,i4)*zb(i3,i1)*zb(i5,i3))/
     &   (zb(i2,i1)*zb(i3,i2)) -
     &  (2*bp*za(i1,i4)*zb(i3,i1)*zb(i5,i3))/
     &   (zb(i2,i1)*zb(i3,i2))
      amp(1,2,1,2)=(za(i1,i4)*zb(i5,i2))/za(i2,i3) -
     &  (2*bp*za(i1,i4)*zb(i5,i2))/za(i2,i3) +
     &  (za(i1,i3)*za(i1,i4)*zb(i5,i3))/
     &   (za(i1,i2)*za(i2,i3)) -
     &  (2*bp*za(i1,i3)*za(i1,i4)*zb(i5,i3))/
     &   (za(i1,i2)*za(i2,i3))
      amp(2,1,1,2)=-((za(i2,i4)*zb(i5,i1))/zb(i3,i2)) +
     &  (2*bp*za(i2,i4)*zb(i5,i1))/zb(i3,i2) -
     &  (za(i3,i4)*zb(i3,i1)*zb(i5,i1))/
     &   (zb(i2,i1)*zb(i3,i2)) +
     &  (2*bp*za(i3,i4)*zb(i3,i1)*zb(i5,i1))/
     &   (zb(i2,i1)*zb(i3,i2))
      amp(2,2,1,2)= -((za(i1,i3)*za(i3,i4)*zb(i5,i1))/
     &     (za(i1,i2)*za(i2,i3))) +
     &  (2*bp*za(i1,i3)*za(i3,i4)*zb(i5,i1))/
     &   (za(i1,i2)*za(i2,i3)) -
     &  (za(i3,i4)*zb(i5,i2))/za(i1,i2) +
     &  (2*bp*za(i3,i4)*zb(i5,i2))/za(i1,i2)


      amp(1,1,1,1)=  -((xmass*za(i2,i4)*zb(i4,i3))/
     &     (zb(i2,i1)*zb(i5,i4))) -
     &  (xmass*za(i1,i4)*zb(i3,i1)*zb(i4,i3))/
     &   (zb(i2,i1)*zb(i3,i2)*zb(i5,i4)) -
     &  (xmass*za(i2,i5)*zb(i5,i3))/
     &   (zb(i2,i1)*zb(i5,i4)) -
     &  (xmass*za(i1,i5)*zb(i3,i1)*zb(i5,i3))/
     &   (zb(i2,i1)*zb(i3,i2)*zb(i5,i4))
      amp(1,2,1,1)=-((xmass*za(i1,i4)*zb(i4,i2))/
     &     (za(i2,i3)*zb(i5,i4))) -
     &  (xmass*za(i1,i3)*za(i1,i4)*zb(i4,i3))/
     &   (za(i1,i2)*za(i2,i3)*zb(i5,i4)) -
     &  (xmass*za(i1,i5)*zb(i5,i2))/
     &   (za(i2,i3)*zb(i5,i4)) -
     &  (xmass*za(i1,i3)*za(i1,i5)*zb(i5,i3))/
     &   (za(i1,i2)*za(i2,i3)*zb(i5,i4))
      amp(2,1,1,1)=  (xmass*za(i2,i4)*zb(i4,i1))/
     &   (zb(i3,i2)*zb(i5,i4)) +
     &  (xmass*za(i3,i4)*zb(i3,i1)*zb(i4,i1))/
     &   (zb(i2,i1)*zb(i3,i2)*zb(i5,i4)) +
     &  (xmass*za(i2,i5)*zb(i5,i1))/
     &   (zb(i3,i2)*zb(i5,i4)) +
     &  (xmass*za(i3,i5)*zb(i3,i1)*zb(i5,i1))/
     &   (zb(i2,i1)*zb(i3,i2)*zb(i5,i4))
      amp(2,2,1,1)=(xmass*za(i1,i3)*za(i3,i4)*zb(i4,i1))/
     &   (za(i1,i2)*za(i2,i3)*zb(i5,i4)) +
     &  (xmass*za(i3,i4)*zb(i4,i2))/
     &   (za(i1,i2)*zb(i5,i4)) +
     &  (xmass*za(i1,i3)*za(i3,i5)*zb(i5,i1))/
     &   (za(i1,i2)*za(i2,i3)*zb(i5,i4)) +
     &  (xmass*za(i3,i5)*zb(i5,i2))/(za(i1,i2)*zb(i5,i4))


      amp(1,1,2,2)=-((xmass*za(i2,i4)*zb(i4,i3))/
     &     (za(i4,i5)*zb(i2,i1))) -
     &  (xmass*za(i1,i4)*zb(i3,i1)*zb(i4,i3))/
     &   (za(i4,i5)*zb(i2,i1)*zb(i3,i2)) -
     &  (xmass*za(i2,i5)*zb(i5,i3))/
     &   (za(i4,i5)*zb(i2,i1)) -
     &  (xmass*za(i1,i5)*zb(i3,i1)*zb(i5,i3))/
     &   (za(i4,i5)*zb(i2,i1)*zb(i3,i2))
      amp(1,2,2,2)=-((xmass*za(i1,i4)*zb(i4,i2))/
     &     (za(i2,i3)*za(i4,i5))) -
     &  (xmass*za(i1,i3)*za(i1,i4)*zb(i4,i3))/
     &   (za(i1,i2)*za(i2,i3)*za(i4,i5)) -
     &  (xmass*za(i1,i5)*zb(i5,i2))/
     &   (za(i2,i3)*za(i4,i5)) -
     &  (xmass*za(i1,i3)*za(i1,i5)*zb(i5,i3))/
     &   (za(i1,i2)*za(i2,i3)*za(i4,i5))
      amp(2,1,2,2)=(xmass*za(i2,i4)*zb(i4,i1))/
     &   (za(i4,i5)*zb(i3,i2)) +
     &  (xmass*za(i3,i4)*zb(i3,i1)*zb(i4,i1))/
     &   (za(i4,i5)*zb(i2,i1)*zb(i3,i2)) +
     &  (xmass*za(i2,i5)*zb(i5,i1))/
     &   (za(i4,i5)*zb(i3,i2)) +
     &  (xmass*za(i3,i5)*zb(i3,i1)*zb(i5,i1))/
     &   (za(i4,i5)*zb(i2,i1)*zb(i3,i2))
      amp(2,2,2,2)=(xmass*za(i1,i3)*za(i3,i4)*zb(i4,i1))/
     &   (za(i1,i2)*za(i2,i3)*za(i4,i5)) +
     &  (xmass*za(i3,i4)*zb(i4,i2))/
     &   (za(i1,i2)*za(i4,i5)) +
     &  (xmass*za(i1,i3)*za(i3,i5)*zb(i5,i1))/
     &   (za(i1,i2)*za(i2,i3)*za(i4,i5)) +
     &  (xmass*za(i3,i5)*zb(i5,i2))/(za(i1,i2)*za(i4,i5))

c$$$      amp_test_n(1,1,2,1)=-((za(i2,i4)*zb(i5,i3))/zb(i2,i1)) +
c$$$     &  (2*bp*za(i2,i4)*zb(i5,i3))/zb(i2,i1) -
c$$$     &  (za(i1,i4)*zb(i3,i1)*zb(i5,i3))/(zb(i2,i1)*zb(i3,i2)) +
c$$$     &  (2*bp*za(i1,i4)*zb(i3,i1)*zb(i5,i3))/(zb(i2,i1)*zb(i3,i2))
c$$$      amp_test_n(1,2,2,1)=-((za(i1,i4)*zb(i5,i2))/za(i2,i3)) +
c$$$     &  (2*bp*za(i1,i4)*zb(i5,i2))/za(i2,i3) -
c$$$     &  (za(i1,i3)*za(i1,i4)*zb(i5,i3))/(za(i1,i2)*za(i2,i3)) +
c$$$     &  (2*bp*za(i1,i3)*za(i1,i4)*zb(i5,i3))/(za(i1,i2)*za(i2,i3))
c$$$      amp_test_n(2,1,2,1)=(za(i2,i4)*zb(i5,i1))/zb(i3,i2) -
c$$$     &  (2*bp*za(i2,i4)*zb(i5,i1))/zb(i3,i2) +
c$$$     &  (za(i3,i4)*zb(i3,i1)*zb(i5,i1))/(zb(i2,i1)*zb(i3,i2)) -
c$$$     &  (2*bp*za(i3,i4)*zb(i3,i1)*zb(i5,i1))/(zb(i2,i1)*zb(i3,i2))
c$$$      amp_test_n(2,2,2,1)=(za(i1,i3)*za(i3,i4)*zb(i5,i1))/
c$$$     &   (za(i1,i2)*za(i2,i3)) -
c$$$     &  (2*bp*za(i1,i3)*za(i3,i4)*zb(i5,i1))/
c$$$     &   (za(i1,i2)*za(i2,i3)) +
c$$$     &  (za(i3,i4)*zb(i5,i2))/za(i1,i2) -
c$$$     &  (2*bp*za(i3,i4)*zb(i5,i2))/za(i1,i2)
c$$$      amp_test_n(1,1,1,2)= (za(i2,i5)*zb(i4,i3))/zb(i2,i1) -
c$$$     &  (2*bp*za(i2,i5)*zb(i4,i3))/zb(i2,i1) +
c$$$     &  (za(i1,i5)*zb(i3,i1)*zb(i4,i3))/(zb(i2,i1)*zb(i3,i2)) -
c$$$     &  (2*bp*za(i1,i5)*zb(i3,i1)*zb(i4,i3))/(zb(i2,i1)*zb(i3,i2))
c$$$      amp_test_n(1,2,1,2)=(za(i1,i5)*zb(i4,i2))/za(i2,i3) -
c$$$     &  (2*bp*za(i1,i5)*zb(i4,i2))/za(i2,i3) +
c$$$     &  (za(i1,i3)*za(i1,i5)*zb(i4,i3))/(za(i1,i2)*za(i2,i3)) -
c$$$     &  (2*bp*za(i1,i3)*za(i1,i5)*zb(i4,i3))/(za(i1,i2)*za(i2,i3))
c$$$      amp_test_n(2,1,1,2)=  -((za(i2,i5)*zb(i4,i1))/zb(i3,i2)) +
c$$$     &  (2*bp*za(i2,i5)*zb(i4,i1))/zb(i3,i2) -
c$$$     &  (za(i3,i5)*zb(i3,i1)*zb(i4,i1))/(zb(i2,i1)*zb(i3,i2)) +
c$$$     &  (2*bp*za(i3,i5)*zb(i3,i1)*zb(i4,i1))/(zb(i2,i1)*zb(i3,i2))
c$$$      amp_test_n(2,2,1,2)=-((za(i1,i3)*za(i3,i5)*zb(i4,i1))/
c$$$     &     (za(i1,i2)*za(i2,i3))) +
c$$$     &  (2*bp*za(i1,i3)*za(i3,i5)*zb(i4,i1))/
c$$$     &   (za(i1,i2)*za(i2,i3)) -
c$$$     &  (za(i3,i5)*zb(i4,i2))/za(i1,i2) +
c$$$     &  (2*bp*za(i3,i5)*zb(i4,i2))/za(i1,i2)
c$$$
c$$$!========= Helicity violating
c$$$
c$$$      amp_test_n(1,1,1,1)=-((xmass*za(i2,i4)*zb(i4,i3))/
c$$$     &     (zb(i2,i1)*zb(i5,i4))) -
c$$$     &  (xmass*za(i1,i4)*zb(i3,i1)*zb(i4,i3))/
c$$$     &   (zb(i2,i1)*zb(i3,i2)*zb(i5,i4)) -
c$$$     &  (xmass*za(i2,i5)*zb(i5,i3))/
c$$$     &   (zb(i2,i1)*zb(i5,i4)) -
c$$$     &  (xmass*za(i1,i5)*zb(i3,i1)*zb(i5,i3))/
c$$$     &   (zb(i2,i1)*zb(i3,i2)*zb(i5,i4))
c$$$      amp_test_n(1,2,1,1)=-((xmass*za(i1,i4)*zb(i4,i2))/
c$$$     &     (za(i2,i3)*zb(i5,i4))) -
c$$$     &  (xmass*za(i1,i3)*za(i1,i4)*zb(i4,i3))/
c$$$     &   (za(i1,i2)*za(i2,i3)*zb(i5,i4)) -
c$$$     &  (xmass*za(i1,i5)*zb(i5,i2))/
c$$$     &   (za(i2,i3)*zb(i5,i4)) -
c$$$     &  (xmass*za(i1,i3)*za(i1,i5)*zb(i5,i3))/
c$$$     &   (za(i1,i2)*za(i2,i3)*zb(i5,i4))
c$$$      amp_test_n(2,1,1,1)=(xmass*za(i2,i4)*zb(i4,i1))/
c$$$     &   (zb(i3,i2)*zb(i5,i4)) +
c$$$     &  (xmass*za(i3,i4)*zb(i3,i1)*zb(i4,i1))/
c$$$     &   (zb(i2,i1)*zb(i3,i2)*zb(i5,i4)) +
c$$$     &  (xmass*za(i2,i5)*zb(i5,i1))/
c$$$     &   (zb(i3,i2)*zb(i5,i4)) +
c$$$     &  (xmass*za(i3,i5)*zb(i3,i1)*zb(i5,i1))/
c$$$     &   (zb(i2,i1)*zb(i3,i2)*zb(i5,i4))
c$$$      amp_test_n(2,2,1,1)=(xmass*za(i1,i3)*za(i3,i4)*zb(i4,i1))/
c$$$     &   (za(i1,i2)*za(i2,i3)*zb(i5,i4)) +
c$$$     &  (xmass*za(i3,i4)*zb(i4,i2))/
c$$$     &   (za(i1,i2)*zb(i5,i4)) +
c$$$     &  (xmass*za(i1,i3)*za(i3,i5)*zb(i5,i1))/
c$$$     &   (za(i1,i2)*za(i2,i3)*zb(i5,i4)) +
c$$$     &  (xmass*za(i3,i5)*zb(i5,i2))/
c$$$     &   (za(i1,i2)*zb(i5,i4))
c$$$
c$$$      amp_test_n(1,1,2,2)=  -((xmass*za(i2,i4)*zb(i4,i3))/
c$$$     &     (za(i4,i5)*zb(i2,i1))) -
c$$$     &  (xmass*za(i1,i4)*zb(i3,i1)*zb(i4,i3))/
c$$$     &   (za(i4,i5)*zb(i2,i1)*zb(i3,i2)) -
c$$$     &  (xmass*za(i2,i5)*zb(i5,i3))/
c$$$     &   (za(i4,i5)*zb(i2,i1)) -
c$$$     &  (xmass*za(i1,i5)*zb(i3,i1)*zb(i5,i3))/
c$$$     &   (za(i4,i5)*zb(i2,i1)*zb(i3,i2))
c$$$      amp_test_n(1,2,2,2)=-((xmass*za(i1,i4)*zb(i4,i2))/
c$$$     &     (za(i2,i3)*za(i4,i5))) -
c$$$     &  (xmass*za(i1,i3)*za(i1,i4)*zb(i4,i3))/
c$$$     &   (za(i1,i2)*za(i2,i3)*za(i4,i5)) -
c$$$     &  (xmass*za(i1,i5)*zb(i5,i2))/
c$$$     &   (za(i2,i3)*za(i4,i5)) -
c$$$     &  (xmass*za(i1,i3)*za(i1,i5)*zb(i5,i3))/
c$$$     &   (za(i1,i2)*za(i2,i3)*za(i4,i5))
c$$$      amp_test_n(2,1,2,2)=  (xmass*za(i2,i4)*zb(i4,i1))/
c$$$     &   (za(i4,i5)*zb(i3,i2)) +
c$$$     &  (xmass*za(i3,i4)*zb(i3,i1)*zb(i4,i1))/
c$$$     &   (za(i4,i5)*zb(i2,i1)*zb(i3,i2)) +
c$$$     &  (xmass*za(i2,i5)*zb(i5,i1))/
c$$$     &   (za(i4,i5)*zb(i3,i2)) +
c$$$     &  (xmass*za(i3,i5)*zb(i3,i1)*zb(i5,i1))/
c$$$     &   (za(i4,i5)*zb(i2,i1)*zb(i3,i2))
c$$$      amp_test_n(2,2,2,2)=(xmass*za(i1,i3)*za(i3,i4)*zb(i4,i1))/
c$$$     &   (za(i1,i2)*za(i2,i3)*za(i4,i5)) +
c$$$     &  (xmass*za(i3,i4)*zb(i4,i2))/
c$$$     &   (za(i1,i2)*za(i4,i5)) +
c$$$     &  (xmass*za(i1,i3)*za(i3,i5)*zb(i5,i1))/
c$$$     &   (za(i1,i2)*za(i2,i3)*za(i4,i5)) +
c$$$     &  (xmass*za(i3,i5)*zb(i5,i2))/
c$$$     &   (za(i1,i2)*za(i4,i5))
c$$$
c$$$      do h1=1,2
c$$$         do h2=1,2
c$$$            do h3=1,2
c$$$               do h4=1,2
c$$$                  write(6,*) h1,h2,h3,h4
c$$$                  write(6,*) amp(h1,h2,h3,h4),amp_test_n(h1,h2,h3,h4),
c$$$     &                 amp(h1,h2,h3,h4)/amp_test_n(h1,h2,h3,h4)
c$$$               enddo
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$      pause


      return
      end
