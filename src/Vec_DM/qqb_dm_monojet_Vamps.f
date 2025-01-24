!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_dm_monojet_Vamps(p,i1,i2,i3,i4,i5,amp)
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


c------- helicity conserving amplitudes (can be checked against Z)
      amp(1,1,1,2)=(za(i2,i4)*zb(i5,i3))/zb(i2,i1) +
     &  (za(i1,i4)*zb(i3,i1)*zb(i5,i3))/(zb(i2,i1)*zb(i3,i2))
      amp(1,2,1,2)=(za(i1,i4)*zb(i5,i2))/za(i2,i3) +
     &  (za(i1,i3)*za(i1,i4)*zb(i5,i3))/(za(i1,i2)*za(i2,i3))
      amp(1,1,2,1)=(za(i2,i5)*zb(i4,i3))/zb(i2,i1) +
     &  (za(i1,i5)*zb(i3,i1)*zb(i4,i3))/(zb(i2,i1)*zb(i3,i2))
      amp(1,2,2,1)=(za(i1,i5)*zb(i4,i2))/za(i2,i3) +
     &  (za(i1,i3)*za(i1,i5)*zb(i4,i3))/(za(i1,i2)*za(i2,i3))
      amp(2,1,1,2)= -((za(i2,i4)*zb(i5,i1))/zb(i3,i2)) -
     &  (za(i3,i4)*zb(i3,i1)*zb(i5,i1))/(zb(i2,i1)*zb(i3,i2))
      amp(2,2,1,2)=-((za(i1,i3)*za(i3,i4)*zb(i5,i1))
     &/(za(i1,i2)*za(i2,i3))) - (za(i3,i4)*zb(i5,i2))/za(i1,i2)
      amp(2,1,2,1)=-((za(i2,i5)*zb(i4,i1))/zb(i3,i2)) -
     &  (za(i3,i5)*zb(i3,i1)*zb(i4,i1))/(zb(i2,i1)*zb(i3,i2))
      amp(2,2,2,1)=-((za(i1,i3)*za(i3,i5)*zb(i4,i1))
     &/(za(i1,i2)*za(i2,i3))) - (za(i3,i5)*zb(i4,i2))/za(i1,i2)

c------ helicity violating amplitudes
      amp(1,1,1,1)=-((xmass*za(i2,i4)*zb(i4,i3))/(zb(i2,i1)*zb(i5,i4)))-
     &  (xmass*za(i1,i4)*zb(i3,i1)*zb(i4,i3))/
     &   (zb(i2,i1)*zb(i3,i2)*zb(i5,i4)) +
     &  (xmass*za(i2,i5)*zb(i5,i3))/(zb(i2,i1)*zb(i5,i4)) +
     &  (xmass*za(i1,i5)*zb(i3,i1)*zb(i5,i3))/
     &   (zb(i2,i1)*zb(i3,i2)*zb(i5,i4))
      amp(1,2,1,1)=-((xmass*za(i1,i4)*zb(i4,i2))/(za(i2,i3)*zb(i5,i4)))-
     &  (xmass*za(i1,i3)*za(i1,i4)*zb(i4,i3))/
     &   (za(i1,i2)*za(i2,i3)*zb(i5,i4)) +
     &  (xmass*za(i1,i5)*zb(i5,i2))/(za(i2,i3)*zb(i5,i4)) +
     &  (xmass*za(i1,i3)*za(i1,i5)*zb(i5,i3))/
     &   (za(i1,i2)*za(i2,i3)*zb(i5,i4))
      amp(2,1,1,1)=(xmass*za(i2,i4)*zb(i4,i1))/(zb(i3,i2)*zb(i5,i4)) +
     &  (xmass*za(i3,i4)*zb(i3,i1)*zb(i4,i1))/
     &   (zb(i2,i1)*zb(i3,i2)*zb(i5,i4)) -
     &  (xmass*za(i2,i5)*zb(i5,i1))/(zb(i3,i2)*zb(i5,i4)) -
     &  (xmass*za(i3,i5)*zb(i3,i1)*zb(i5,i1))/
     &   (zb(i2,i1)*zb(i3,i2)*zb(i5,i4))
      amp(2,2,1,1)=(xmass*za(i1,i3)*za(i3,i4)*zb(i4,i1))/
     &   (za(i1,i2)*za(i2,i3)*zb(i5,i4)) +
     &  (xmass*za(i3,i4)*zb(i4,i2))/(za(i1,i2)*zb(i5,i4)) -
     &  (xmass*za(i1,i3)*za(i3,i5)*zb(i5,i1))/
     &   (za(i1,i2)*za(i2,i3)*zb(i5,i4)) -
     &  (xmass*za(i3,i5)*zb(i5,i2))/(za(i1,i2)*zb(i5,i4))
      amp(1,1,2,2)=(xmass*za(i2,i4)*zb(i4,i3))/(za(i4,i5)*zb(i2,i1)) +
     &  (xmass*za(i1,i4)*zb(i3,i1)*zb(i4,i3))/
     &   (za(i4,i5)*zb(i2,i1)*zb(i3,i2)) -
     &  (xmass*za(i2,i5)*zb(i5,i3))/(za(i4,i5)*zb(i2,i1)) -
     &  (xmass*za(i1,i5)*zb(i3,i1)*zb(i5,i3))/
     &   (za(i4,i5)*zb(i2,i1)*zb(i3,i2))
      amp(1,2,2,2)=(xmass*za(i1,i4)*zb(i4,i2))/(za(i2,i3)*za(i4,i5)) +
     &  (xmass*za(i1,i3)*za(i1,i4)*zb(i4,i3))/
     &   (za(i1,i2)*za(i2,i3)*za(i4,i5)) -
     &  (xmass*za(i1,i5)*zb(i5,i2))/(za(i2,i3)*za(i4,i5)) -
     &  (xmass*za(i1,i3)*za(i1,i5)*zb(i5,i3))/
     &   (za(i1,i2)*za(i2,i3)*za(i4,i5))
      amp(2,1,2,2)=-((xmass*za(i2,i4)*zb(i4,i1))/(za(i4,i5)*zb(i3,i2)))-
     &  (xmass*za(i3,i4)*zb(i3,i1)*zb(i4,i1))/
     &   (za(i4,i5)*zb(i2,i1)*zb(i3,i2)) +
     &  (xmass*za(i2,i5)*zb(i5,i1))/(za(i4,i5)*zb(i3,i2)) +
     &  (xmass*za(i3,i5)*zb(i3,i1)*zb(i5,i1))/
     &   (za(i4,i5)*zb(i2,i1)*zb(i3,i2))
      amp(2,2,2,2)= -((xmass*za(i1,i3)*za(i3,i4)*zb(i4,i1))/
     &     (za(i1,i2)*za(i2,i3)*za(i4,i5))) -
     &  (xmass*za(i3,i4)*zb(i4,i2))/(za(i1,i2)*za(i4,i5)) +
     &  (xmass*za(i1,i3)*za(i3,i5)*zb(i5,i1))/
     &   (za(i1,i2)*za(i2,i3)*za(i4,i5)) +
     &  (xmass*za(i3,i5)*zb(i5,i2))/(za(i1,i2)*za(i4,i5))


c------ define square
c      amp2=0d0
c      do h1=1,2
c         do h2=1,2
c            do h3=1,2
c               do h4=1,2
c                  amp2=amp2+abs(dm(L)*amp(1,h2,h3,h4))**2

c               enddo
c            enddo
c         enddo
c      enddo

c      write(6,*) amp2

c---- check against z
c      do h1=1,2
c         do h2=1,2
c            amp_tes(h1,h2)=abs(amp(h1,1,h2,3-h2))**2
c     &           +abs(amp(h1,2,h2,3-h2))**2
c         enddo
c      enddo

c      write(6,*) amp_tes(1,1)**2+amp_tes(1,2)**2
c     &+amp_tes(2,1)**2+amp_tes(2,2)**2

c      write(6,*) i1,i2,i3,i4,i5,amp_tes
c      write(6,*) abs(amp(1,1,1,2))**2+abs(amp(1,2,1,2))**2

c      call zgamps2(i1,i3,i5,i4,i2,za,zb,z1jet_tes)
c      write(6,*) '1, 1',z1jet_tes(1,1)
c      write(6,*) '2, 2',z1jet_tes(2,2)
c      write(6,*) '2, 1',z1jet_tes(2,1)
c      write(6,*) '1, 2',z1jet_tes(1,2)
c      write(6,*) h1,h2
c      write(6,*) ' 1, 1 mine  ',amp_tes(1,1)
c      write(6,*) ' 1, 2 mine  ',amp_tes(1,2)
c      write(6,*) ' 2, 1 mine  ',amp_tes(2,1)
c      write(6,*) ' 2, 2 mine  ',amp_tes(2,2)
c      pause
      return
      end
