!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_dm_gg_Vamps(p,i1,i2,i3,i4,i5,i6,msq1,msq2,msqsl)
      implicit none
      include 'types.f'
      include 'dm_params.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
c----- vector amplitude for
c----- q(i1)+g(i2)+g(i3)+qb(i4)+x(i5)+x(i6)
      real(dp):: p(mxpart,4),q(mxpart,4)
c-----fills amplitude for q g qb chi,chib
      integer:: i1,i2,i3,i4,i5,i6
c------and gluon helicity is summed over
c------returns amp squared.
      integer:: h1,h2,h3,h4,h5
c      complex(dp):: cor_fac
      complex(dp):: amp_hconsx1(2,2,2,2,2),amp_hconsx2(2,2,2,2,2)
      complex(dp):: amp_hconsQED(2,2,2,2,2)
c------ index refers to helicit of qqb (sum over everyting else)
      real(dp):: msq1(2),msq2(2),msqsl(2)


      if(xmass>1d-8) then
c--------- generate massless phase space
      call gen_masslessvecs(p,q,i5,i6)
c--------- generate spinors
      call spinoru(6,q,za,zb)
      else
c-------- massless dm can use usual spinoru
         call spinoru(6,p,za,zb)
      endif

c====== helicity amplitudes
      call dm_gg_helamps(i1,i2,i3,i4,i5,i6,za,zb,amp_hconsx1)
c======= other colour ordering, swap gluons
      call dm_gg_helamps(i1,i3,i2,i4,i5,i6,za,zb,amp_hconsx2)
      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
                  do h5=1,2
               amp_hconsQED(h1,h2,h3,h4,h5)=amp_hconsx1(h1,h2,h3,h4,h5)
     &                +amp_hconsx2(h1,h3,h2,h4,h5)
             enddo
          enddo
       enddo
      enddo
      enddo

      msq1(:)=zero
      msq2(:)=zero
      msqsl(:)=zero
c------ msq
      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
                  do h5=1,2
              msq1(h1)=msq1(h1)+abs(amp_hconsx1(h1,h2,h3,h4,h5))**2
              msq2(h1)=msq2(h1)+abs(amp_hconsx2(h1,h2,h3,h4,h5))**2
              msqsl(h1)=msqsl(h1)+abs(amp_hconsQED(h1,h2,h3,h4,h5))**2
                  enddo
               enddo
            enddo
         enddo
      enddo
c----- include -1/N_c^2 for subleading piece
      do h1=1,2
         msqsl(h1)=-one/9._dp*msqsl(h1)
      enddo

c---- order is q qb -> l l g g
c---- z + 2j test
c      write(6,*) i1,i2,i3,i4,i5,i6
c      call subqcd(i1,i4,i5,i6,i2,i3,za,zb,z2jet_1LL)


c      cor_fac=4._dp/(za(i5,i6)*zb(i6,i5))
c      write(6,*) 'my -1 -1', abs(amp(2,1,1,1,2)*cor_fac)
c      write(6,*) 'my 1 1', abs(amp(2,2,2,1,2)*cor_fac)
c      write(6,*) 'my -1 1', abs(amp(2,1,2,1,2)*cor_fac)
c      write(6,*) 'my 1 -1', abs(amp(2,2,1,1,2)*cor_fac)

c      write(6,*) 'test -1 -1 ',abs(z2jet_1LL(-1,-1))
c      write(6,*) 'test 1 1 ',abs(z2jet_1LL(1,1))
c      write(6,*) 'test -1 1 ',abs(z2jet_1LL(-1,1))
c      write(6,*) 'test 1 -1 ',abs(z2jet_1LL(1,-1))
c      pause
      return
      end


      subroutine dm_gg_helamps(i1,i2,i3,i4,i5,i6,za,zb,amp)
      implicit none
      include 'types.f'
      include 'dm_params.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
c----- vector amplitude for
c----- q(i1)+g(i2)+g(i3)+qb(i4)+x(i5)+x(i6)
c      real(dp):: p(mxpart,4)
c-----fills amplitude for q g qb chi,chib
      complex(dp):: amp(2,2,2,2,2)
c      real(dp):: q(mxpart,4)
      integer:: i1,i2,i3,i4,i5,i6
c------and gluon helicity is summed over
c------returns amp squared.
c      complex(dp):: z2jet_1LL(-1:1,-1:1)
c      real(dp):: amp2,amp_tes(2,2)
      integer:: h1,h2
c      complex(dp):: cor_fac
      do h1=1,6
         do h2=1,6
            s(h1,h2)=Dble(za(h1,h2)*zb(h2,h1))
         enddo
      enddo

      amp(:,:,:,:,:)=czip
c------- helicity conserving amplitudes
      amp(2,1,1,1,2)=(za(i2,i5)*zb(i6,i1))/(zb(i3,i2)*zb(i4,i3)) +
     &  (za(i3,i5)*zb(i3,i1)*zb(i6,i1))/
     &   (zb(i2,i1)*zb(i3,i2)*zb(i4,i3)) +
     &  (za(i4,i5)*zb(i4,i1)*zb(i6,i1))/(zb(i2,i1)*zb(i3,i2)*zb(i4,i3))
      amp(2,2,2,1,2)=(za(i1,i4)*za(i4,i5)*zb(i6,i1))/
     &   (za(i1,i2)*za(i2,i3)*za(i3,i4)) +
     &  (za(i2,i4)*za(i4,i5)*zb(i6,i2))/
     &   (za(i1,i2)*za(i2,i3)*za(i3,i4)) +
     &  (za(i4,i5)*zb(i6,i3))/(za(i1,i2)*za(i2,i3))
      amp(2,1,2,1,2)=-((za(i2,i4)*za(i4,i5)*zb(i3,i1)*zb(i6,i1))/
     &     (s(i2,i3)*za(i3,i4)*zb(i2,i1))) +
     &  (za(i1,i2)*za(i4,i5)*zb(i3,i1)**2*zb(i6,i1))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*zb(i2,i1)) -
     &  (za(i2,i4)**2*za(i2,i5)*zb(i3,i2)*zb(i6,i1))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i3,i4)) +
     &  (za(i2,i4)**2*za(i4,i5)*zb(i4,i3)*zb(i6,i1))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i3,i4)) -
     &  (za(i2,i3)*za(i4,i5)*zb(i3,i1)**2*zb(i6,i3))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*zb(i2,i1))
      amp(2,2,1,1,2)=(za(i1,i3)**2*za(i4,i5)*zb(i2,i1)*zb(i6,i1))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i1,i2)) -
     &  (za(i1,i3)*za(i3,i5)*zb(i3,i2)*zb(i6,i1))/
     &   (s(i2,i3)*za(i1,i2)*zb(i4,i3)) -
     &  (za(i1,i3)*za(i4,i5)*zb(i4,i2)*zb(i6,i1))/
     &   (s(i2,i3)*za(i1,i2)*zb(i4,i3)) +
     &  (za(i3,i4)*za(i3,i5)*zb(i3,i2)*zb(i4,i2)*zb(i6,i1))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i4,i3)) +
     &  (za(i3,i4)*za(i4,i5)*zb(i4,i2)**2*zb(i6,i1))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i4,i3)) +
     &  (za(i1,i3)*za(i2,i3)*za(i4,i5)*zb(i2,i1)*zb(i6,i2))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i1,i2)) -
     &  (za(i2,i3)*za(i3,i5)*zb(i3,i2)*zb(i6,i2))/
     &   (s(i2,i3)*za(i1,i2)*zb(i4,i3)) -
     &  (za(i2,i3)*za(i4,i5)*zb(i4,i2)*zb(i6,i2))/
     &   (s(i2,i3)*za(i1,i2)*zb(i4,i3))
c------ current flip : quarks
      amp(1,1,1,1,2)=(za(i3,i5)*zb(i6,i4))/(zb(i2,i1)*zb(i3,i2)) +
     &  (za(i1,i5)*zb(i4,i1)*zb(i6,i4))/
     &   (zb(i2,i1)*zb(i3,i2)*zb(i4,i3)) +
     &  (za(i2,i5)*zb(i4,i2)*zb(i6,i4))/(zb(i2,i1)*zb(i3,i2)*zb(i4,i3))
      amp(1,1,2,1,2)=-((za(i1,i5)*za(i2,i3)*zb(i3,i1)*zb(i6,i3))/
     &     (s(i2,i3)*za(i3,i4)*zb(i2,i1))) -
     &  (za(i2,i3)*za(i2,i5)*zb(i3,i2)*zb(i6,i3))/
     &   (s(i2,i3)*za(i3,i4)*zb(i2,i1)) +
     &  (za(i1,i5)*za(i2,i3)*za(i2,i4)*zb(i4,i3)*zb(i6,i3))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i3,i4)) -
     &  (za(i1,i5)*za(i2,i4)*zb(i3,i1)*zb(i6,i4))/
     &   (s(i2,i3)*za(i3,i4)*zb(i2,i1)) +
     &  (za(i1,i2)*za(i1,i5)*zb(i3,i1)**2*zb(i6,i4))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*zb(i2,i1)) -
     &  (za(i2,i4)*za(i2,i5)*zb(i3,i2)*zb(i6,i4))/
     &   (s(i2,i3)*za(i3,i4)*zb(i2,i1)) +
     &  (za(i1,i2)*za(i2,i5)*zb(i3,i1)*zb(i3,i2)*zb(i6,i4))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*zb(i2,i1)) +
     &  (za(i1,i5)*za(i2,i4)**2*zb(i4,i3)*zb(i6,i4))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i3,i4))
      amp(1,2,1,1,2)= -((za(i1,i5)*za(i2,i3)*zb(i4,i2)**2*zb(i6,i2))/
     &     (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i4,i3))) +
     &  (za(i1,i3)**2*za(i1,i5)*zb(i2,i1)*zb(i6,i4))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i1,i2)) -
     &  (za(i1,i3)**2*za(i3,i5)*zb(i3,i2)*zb(i6,i4))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i1,i2)) -
     &  (za(i1,i3)*za(i1,i5)*zb(i4,i2)*zb(i6,i4))/
     &   (s(i2,i3)*za(i1,i2)*zb(i4,i3)) +
     &  (za(i1,i5)*za(i3,i4)*zb(i4,i2)**2*zb(i6,i4))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i4,i3))
      amp(1,2,2,1,2)=(za(i1,i5)*zb(i6,i2))/(za(i2,i3)*za(i3,i4)) +
     &  (za(i1,i3)*za(i1,i5)*zb(i6,i3))/
     &   (za(i1,i2)*za(i2,i3)*za(i3,i4)) +
     &  (za(i1,i4)*za(i1,i5)*zb(i6,i4))/(za(i1,i2)*za(i2,i3)*za(i3,i4))

c----- current flip : dm
      amp(2,1,1,2,1)= (za(i2,i6)*zb(i5,i1))/(zb(i3,i2)*zb(i4,i3)) +
     &  (za(i3,i6)*zb(i3,i1)*zb(i5,i1))/
     &   (zb(i2,i1)*zb(i3,i2)*zb(i4,i3)) +
     &  (za(i4,i6)*zb(i4,i1)*zb(i5,i1))/(zb(i2,i1)*zb(i3,i2)*zb(i4,i3))
      amp(2,1,2,2,1)=-((za(i2,i4)*za(i4,i6)*zb(i3,i1)*zb(i5,i1))/
     &     (s(i2,i3)*za(i3,i4)*zb(i2,i1))) +
     &  (za(i1,i2)*za(i4,i6)*zb(i3,i1)**2*zb(i5,i1))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*zb(i2,i1)) -
     &  (za(i2,i4)**2*za(i2,i6)*zb(i3,i2)*zb(i5,i1))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i3,i4)) +
     &  (za(i2,i4)**2*za(i4,i6)*zb(i4,i3)*zb(i5,i1))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i3,i4)) -
     &  (za(i2,i3)*za(i4,i6)*zb(i3,i1)**2*zb(i5,i3))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*zb(i2,i1))
      amp(2,2,1,2,1)= (za(i1,i3)**2*za(i4,i6)*zb(i2,i1)*zb(i5,i1))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i1,i2)) -
     &  (za(i1,i3)*za(i3,i6)*zb(i3,i2)*zb(i5,i1))/
     &   (s(i2,i3)*za(i1,i2)*zb(i4,i3)) -
     &  (za(i1,i3)*za(i4,i6)*zb(i4,i2)*zb(i5,i1))/
     &   (s(i2,i3)*za(i1,i2)*zb(i4,i3)) +
     &  (za(i3,i4)*za(i3,i6)*zb(i3,i2)*zb(i4,i2)*zb(i5,i1))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i4,i3)) +
     &  (za(i3,i4)*za(i4,i6)*zb(i4,i2)**2*zb(i5,i1))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i4,i3)) +
     &  (za(i1,i3)*za(i2,i3)*za(i4,i6)*zb(i2,i1)*zb(i5,i2))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i1,i2)) -
     &  (za(i2,i3)*za(i3,i6)*zb(i3,i2)*zb(i5,i2))/
     &   (s(i2,i3)*za(i1,i2)*zb(i4,i3)) -
     &  (za(i2,i3)*za(i4,i6)*zb(i4,i2)*zb(i5,i2))/
     &   (s(i2,i3)*za(i1,i2)*zb(i4,i3))
      amp(2,2,2,2,1)=(za(i1,i4)*za(i4,i6)*zb(i5,i1))/
     &   (za(i1,i2)*za(i2,i3)*za(i3,i4)) +
     &  (za(i2,i4)*za(i4,i6)*zb(i5,i2))/
     &   (za(i1,i2)*za(i2,i3)*za(i3,i4)) +
     &  (za(i4,i6)*zb(i5,i3))/(za(i1,i2)*za(i2,i3))

c------ current flip : dm and quark
      amp(1,1,1,2,1)=(za(i3,i6)*zb(i5,i4))/(zb(i2,i1)*zb(i3,i2)) +
     &  (za(i1,i6)*zb(i4,i1)*zb(i5,i4))/
     &   (zb(i2,i1)*zb(i3,i2)*zb(i4,i3)) +
     &  (za(i2,i6)*zb(i4,i2)*zb(i5,i4))/(zb(i2,i1)*zb(i3,i2)*zb(i4,i3))
      amp(1,1,2,2,1)=-((za(i1,i6)*za(i2,i3)*zb(i3,i1)*zb(i5,i3))/
     &     (s(i2,i3)*za(i3,i4)*zb(i2,i1))) -
     &  (za(i2,i3)*za(i2,i6)*zb(i3,i2)*zb(i5,i3))/
     &   (s(i2,i3)*za(i3,i4)*zb(i2,i1)) +
     &  (za(i1,i6)*za(i2,i3)*za(i2,i4)*zb(i4,i3)*zb(i5,i3))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i3,i4)) -
     &  (za(i1,i6)*za(i2,i4)*zb(i3,i1)*zb(i5,i4))/
     &   (s(i2,i3)*za(i3,i4)*zb(i2,i1)) +
     &  (za(i1,i2)*za(i1,i6)*zb(i3,i1)**2*zb(i5,i4))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*zb(i2,i1)) -
     &  (za(i2,i4)*za(i2,i6)*zb(i3,i2)*zb(i5,i4))/
     &   (s(i2,i3)*za(i3,i4)*zb(i2,i1)) +
     &  (za(i1,i2)*za(i2,i6)*zb(i3,i1)*zb(i3,i2)*zb(i5,i4))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*zb(i2,i1)) +
     &  (za(i1,i6)*za(i2,i4)**2*zb(i4,i3)*zb(i5,i4))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i3,i4))
      amp(1,2,1,2,1)=-((za(i1,i6)*za(i2,i3)*zb(i4,i2)**2*zb(i5,i2))/
     &     (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i4,i3))) +
     &  (za(i1,i3)**2*za(i1,i6)*zb(i2,i1)*zb(i5,i4))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i1,i2)) -
     &  (za(i1,i3)**2*za(i3,i6)*zb(i3,i2)*zb(i5,i4))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i1,i2)) -
     &  (za(i1,i3)*za(i1,i6)*zb(i4,i2)*zb(i5,i4))/
     &   (s(i2,i3)*za(i1,i2)*zb(i4,i3)) +
     &  (za(i1,i6)*za(i3,i4)*zb(i4,i2)**2*zb(i5,i4))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i4,i3))
      amp(1,2,2,2,1)=(za(i1,i6)*zb(i5,i2))/(za(i2,i3)*za(i3,i4)) +
     &  (za(i1,i3)*za(i1,i6)*zb(i5,i3))/
     &   (za(i1,i2)*za(i2,i3)*za(i3,i4)) +
     &  (za(i1,i4)*za(i1,i6)*zb(i5,i4))/(za(i1,i2)*za(i2,i3)*za(i3,i4))



c======= HELICITY VIOLATING AMPLITUDES

c------
      amp(1,1,1,1,1)=-((xmass*za(i3,i5)*zb(i5,i4))/
     &     (zb(i2,i1)*zb(i3,i2)*zb(i6,i5))) -
     &  (xmass*za(i1,i5)*zb(i4,i1)*zb(i5,i4))/
     &   (zb(i2,i1)*zb(i3,i2)*zb(i4,i3)*zb(i6,i5)) -
     &  (xmass*za(i2,i5)*zb(i4,i2)*zb(i5,i4))/
     &   (zb(i2,i1)*zb(i3,i2)*zb(i4,i3)*zb(i6,i5)) +
     &  (xmass*za(i3,i6)*zb(i6,i4))/(zb(i2,i1)*zb(i3,i2)*zb(i6,i5)) +
     &  (xmass*za(i1,i6)*zb(i4,i1)*zb(i6,i4))/
     &   (zb(i2,i1)*zb(i3,i2)*zb(i4,i3)*zb(i6,i5)) +
     &  (xmass*za(i2,i6)*zb(i4,i2)*zb(i6,i4))/
     &   (zb(i2,i1)*zb(i3,i2)*zb(i4,i3)*zb(i6,i5))
      amp(1,1,2,1,1)=(xmass*za(i1,i5)*za(i2,i3)*zb(i3,i1)*zb(i5,i3))/
     &   (s(i2,i3)*za(i3,i4)*zb(i2,i1)*zb(i6,i5)) +
     &  (xmass*za(i2,i3)*za(i2,i5)*zb(i3,i2)*zb(i5,i3))/
     &   (s(i2,i3)*za(i3,i4)*zb(i2,i1)*zb(i6,i5)) -
     &  (xmass*za(i1,i5)*za(i2,i3)*za(i2,i4)*zb(i4,i3)*zb(i5,i3))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i3,i4)*zb(i6,i5))
     &    + (xmass*za(i1,i5)*za(i2,i4)*zb(i3,i1)*zb(i5,i4))/
     &   (s(i2,i3)*za(i3,i4)*zb(i2,i1)*zb(i6,i5)) -
     &  (xmass*za(i1,i2)*za(i1,i5)*zb(i3,i1)**2*zb(i5,i4))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*zb(i2,i1)*zb(i6,i5))
     &    + (xmass*za(i2,i4)*za(i2,i5)*zb(i3,i2)*zb(i5,i4))/
     &   (s(i2,i3)*za(i3,i4)*zb(i2,i1)*zb(i6,i5)) -
     &  (xmass*za(i1,i2)*za(i2,i5)*zb(i3,i1)*zb(i3,i2)*zb(i5,i4))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*zb(i2,i1)*zb(i6,i5))
     &    - (xmass*za(i1,i5)*za(i2,i4)**2*zb(i4,i3)*zb(i5,i4))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i3,i4)*zb(i6,i5))
     &    - (xmass*za(i1,i6)*za(i2,i3)*zb(i3,i1)*zb(i6,i3))/
     &   (s(i2,i3)*za(i3,i4)*zb(i2,i1)*zb(i6,i5)) -
     &  (xmass*za(i2,i3)*za(i2,i6)*zb(i3,i2)*zb(i6,i3))/
     &   (s(i2,i3)*za(i3,i4)*zb(i2,i1)*zb(i6,i5)) +
     &  (xmass*za(i1,i6)*za(i2,i3)*za(i2,i4)*zb(i4,i3)*zb(i6,i3))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i3,i4)*zb(i6,i5))
     &    - (xmass*za(i1,i6)*za(i2,i4)*zb(i3,i1)*zb(i6,i4))/
     &   (s(i2,i3)*za(i3,i4)*zb(i2,i1)*zb(i6,i5)) +
     &  (xmass*za(i1,i2)*za(i1,i6)*zb(i3,i1)**2*zb(i6,i4))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*zb(i2,i1)*zb(i6,i5))
     &    - (xmass*za(i2,i4)*za(i2,i6)*zb(i3,i2)*zb(i6,i4))/
     &   (s(i2,i3)*za(i3,i4)*zb(i2,i1)*zb(i6,i5)) +
     &  (xmass*za(i1,i2)*za(i2,i6)*zb(i3,i1)*zb(i3,i2)*zb(i6,i4))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*zb(i2,i1)*zb(i6,i5))
     &    + (xmass*za(i1,i6)*za(i2,i4)**2*zb(i4,i3)*zb(i6,i4))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i3,i4)*zb(i6,i5))
      amp(1,2,1,1,1)=(xmass*za(i1,i5)*za(i2,i3)*zb(i4,i2)**2*zb(i5,i2))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i4,i3)*zb(i6,i5))
     &    - (xmass*za(i1,i3)**2*za(i1,i5)*zb(i2,i1)*zb(i5,i4))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i1,i2)*zb(i6,i5))
     &    + (xmass*za(i1,i3)**2*za(i3,i5)*zb(i3,i2)*zb(i5,i4))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i1,i2)*zb(i6,i5))
     &    + (xmass*za(i1,i3)*za(i1,i5)*zb(i4,i2)*zb(i5,i4))/
     &   (s(i2,i3)*za(i1,i2)*zb(i4,i3)*zb(i6,i5)) -
     &  (xmass*za(i1,i5)*za(i3,i4)*zb(i4,i2)**2*zb(i5,i4))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i4,i3)*zb(i6,i5))
     &    - (xmass*za(i1,i6)*za(i2,i3)*zb(i4,i2)**2*zb(i6,i2))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i4,i3)*zb(i6,i5))
     &    + (xmass*za(i1,i3)**2*za(i1,i6)*zb(i2,i1)*zb(i6,i4))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i1,i2)*zb(i6,i5))
     &    - (xmass*za(i1,i3)**2*za(i3,i6)*zb(i3,i2)*zb(i6,i4))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i1,i2)*zb(i6,i5))
     &    - (xmass*za(i1,i3)*za(i1,i6)*zb(i4,i2)*zb(i6,i4))/
     &   (s(i2,i3)*za(i1,i2)*zb(i4,i3)*zb(i6,i5)) +
     &  (xmass*za(i1,i6)*za(i3,i4)*zb(i4,i2)**2*zb(i6,i4))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i4,i3)*zb(i6,i5))
      amp(2,1,1,1,1)=-((xmass*za(i2,i5)*zb(i5,i1))/
     &     (zb(i3,i2)*zb(i4,i3)*zb(i6,i5))) -
     &  (xmass*za(i3,i5)*zb(i3,i1)*zb(i5,i1))/
     &   (zb(i2,i1)*zb(i3,i2)*zb(i4,i3)*zb(i6,i5)) -
     &  (xmass*za(i4,i5)*zb(i4,i1)*zb(i5,i1))/
     &   (zb(i2,i1)*zb(i3,i2)*zb(i4,i3)*zb(i6,i5)) +
     &  (xmass*za(i2,i6)*zb(i6,i1))/(zb(i3,i2)*zb(i4,i3)*zb(i6,i5)) +
     &  (xmass*za(i3,i6)*zb(i3,i1)*zb(i6,i1))/
     &   (zb(i2,i1)*zb(i3,i2)*zb(i4,i3)*zb(i6,i5)) +
     &  (xmass*za(i4,i6)*zb(i4,i1)*zb(i6,i1))/
     &   (zb(i2,i1)*zb(i3,i2)*zb(i4,i3)*zb(i6,i5))
      amp(1,2,2,1,1)= -((xmass*za(i1,i5)*zb(i5,i2))/
     &     (za(i2,i3)*za(i3,i4)*zb(i6,i5))) -
     &  (xmass*za(i1,i3)*za(i1,i5)*zb(i5,i3))/
     &   (za(i1,i2)*za(i2,i3)*za(i3,i4)*zb(i6,i5)) -
     &  (xmass*za(i1,i4)*za(i1,i5)*zb(i5,i4))/
     &   (za(i1,i2)*za(i2,i3)*za(i3,i4)*zb(i6,i5)) +
     &  (xmass*za(i1,i6)*zb(i6,i2))/(za(i2,i3)*za(i3,i4)*zb(i6,i5)) +
     &  (xmass*za(i1,i3)*za(i1,i6)*zb(i6,i3))/
     &   (za(i1,i2)*za(i2,i3)*za(i3,i4)*zb(i6,i5)) +
     &  (xmass*za(i1,i4)*za(i1,i6)*zb(i6,i4))/
     &   (za(i1,i2)*za(i2,i3)*za(i3,i4)*zb(i6,i5))
      amp(2,1,2,1,1)=(xmass*za(i2,i4)*za(i4,i5)*zb(i3,i1)*zb(i5,i1))/
     &   (s(i2,i3)*za(i3,i4)*zb(i2,i1)*zb(i6,i5)) -
     &  (xmass*za(i1,i2)*za(i4,i5)*zb(i3,i1)**2*zb(i5,i1))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*zb(i2,i1)*zb(i6,i5))
     &    + (xmass*za(i2,i4)**2*za(i2,i5)*zb(i3,i2)*zb(i5,i1))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i3,i4)*zb(i6,i5))
     &    - (xmass*za(i2,i4)**2*za(i4,i5)*zb(i4,i3)*zb(i5,i1))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i3,i4)*zb(i6,i5))
     &    + (xmass*za(i2,i3)*za(i4,i5)*zb(i3,i1)**2*zb(i5,i3))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*zb(i2,i1)*zb(i6,i5))
     &    - (xmass*za(i2,i4)*za(i4,i6)*zb(i3,i1)*zb(i6,i1))/
     &   (s(i2,i3)*za(i3,i4)*zb(i2,i1)*zb(i6,i5)) +
     &  (xmass*za(i1,i2)*za(i4,i6)*zb(i3,i1)**2*zb(i6,i1))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*zb(i2,i1)*zb(i6,i5))
     &    - (xmass*za(i2,i4)**2*za(i2,i6)*zb(i3,i2)*zb(i6,i1))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i3,i4)*zb(i6,i5))
     &    + (xmass*za(i2,i4)**2*za(i4,i6)*zb(i4,i3)*zb(i6,i1))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i3,i4)*zb(i6,i5))
     &    - (xmass*za(i2,i3)*za(i4,i6)*zb(i3,i1)**2*zb(i6,i3))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*zb(i2,i1)*zb(i6,i5))
      amp(2,2,1,1,1)=-((xmass*za(i1,i3)**2*za(i4,i5)*zb(i2,i1)*
     &     zb(i5,i1))/
     &     (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i1,i2)*
     &       zb(i6,i5))) + (xmass*za(i1,i3)*za(i3,i5)*zb(i3,i2)*
     &     zb(i5,i1))/(s(i2,i3)*za(i1,i2)*zb(i4,i3)*zb(i6,i5)) +
     &  (xmass*za(i1,i3)*za(i4,i5)*zb(i4,i2)*zb(i5,i1))/
     &   (s(i2,i3)*za(i1,i2)*zb(i4,i3)*zb(i6,i5)) -
     &  (xmass*za(i3,i4)*za(i3,i5)*zb(i3,i2)*zb(i4,i2)*zb(i5,i1))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i4,i3)*zb(i6,i5))
     &    - (xmass*za(i3,i4)*za(i4,i5)*zb(i4,i2)**2*zb(i5,i1))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i4,i3)*zb(i6,i5))
     &    - (xmass*za(i1,i3)*za(i2,i3)*za(i4,i5)*zb(i2,i1)*zb(i5,i2))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i1,i2)*zb(i6,i5))
     &    + (xmass*za(i2,i3)*za(i3,i5)*zb(i3,i2)*zb(i5,i2))/
     &   (s(i2,i3)*za(i1,i2)*zb(i4,i3)*zb(i6,i5)) +
     &  (xmass*za(i2,i3)*za(i4,i5)*zb(i4,i2)*zb(i5,i2))/
     &   (s(i2,i3)*za(i1,i2)*zb(i4,i3)*zb(i6,i5)) +
     &  (xmass*za(i1,i3)**2*za(i4,i6)*zb(i2,i1)*zb(i6,i1))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i1,i2)*zb(i6,i5))
     &    - (xmass*za(i1,i3)*za(i3,i6)*zb(i3,i2)*zb(i6,i1))/
     &   (s(i2,i3)*za(i1,i2)*zb(i4,i3)*zb(i6,i5)) -
     &  (xmass*za(i1,i3)*za(i4,i6)*zb(i4,i2)*zb(i6,i1))/
     &   (s(i2,i3)*za(i1,i2)*zb(i4,i3)*zb(i6,i5)) +
     &  (xmass*za(i3,i4)*za(i3,i6)*zb(i3,i2)*zb(i4,i2)*zb(i6,i1))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i4,i3)*zb(i6,i5))
     &    + (xmass*za(i3,i4)*za(i4,i6)*zb(i4,i2)**2*zb(i6,i1))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i4,i3)*zb(i6,i5))
     &    + (xmass*za(i1,i3)*za(i2,i3)*za(i4,i6)*zb(i2,i1)*zb(i6,i2))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i1,i2)*zb(i6,i5))
     &    - (xmass*za(i2,i3)*za(i3,i6)*zb(i3,i2)*zb(i6,i2))/
     &   (s(i2,i3)*za(i1,i2)*zb(i4,i3)*zb(i6,i5)) -
     &  (xmass*za(i2,i3)*za(i4,i6)*zb(i4,i2)*zb(i6,i2))/
     &   (s(i2,i3)*za(i1,i2)*zb(i4,i3)*zb(i6,i5))
      amp(2,2,2,1,1)=-((xmass*za(i1,i4)*za(i4,i5)*zb(i5,i1))/
     &     (za(i1,i2)*za(i2,i3)*za(i3,i4)*zb(i6,i5))) -
     &  (xmass*za(i2,i4)*za(i4,i5)*zb(i5,i2))/
     &   (za(i1,i2)*za(i2,i3)*za(i3,i4)*zb(i6,i5)) -
     &  (xmass*za(i4,i5)*zb(i5,i3))/(za(i1,i2)*za(i2,i3)*zb(i6,i5)) +
     &  (xmass*za(i1,i4)*za(i4,i6)*zb(i6,i1))/
     &   (za(i1,i2)*za(i2,i3)*za(i3,i4)*zb(i6,i5)) +
     &  (xmass*za(i2,i4)*za(i4,i6)*zb(i6,i2))/
     &   (za(i1,i2)*za(i2,i3)*za(i3,i4)*zb(i6,i5)) +
     &  (xmass*za(i4,i6)*zb(i6,i3))/(za(i1,i2)*za(i2,i3)*zb(i6,i5))

c----- spin flip on dm
      amp(1,1,1,2,2)=(xmass*za(i3,i5)*zb(i5,i4))/(za(i5,i6)*
     &     zb(i2,i1)*zb(i3,i2)) +
     &  (xmass*za(i1,i5)*zb(i4,i1)*zb(i5,i4))/
     &   (za(i5,i6)*zb(i2,i1)*zb(i3,i2)*zb(i4,i3)) +
     &  (xmass*za(i2,i5)*zb(i4,i2)*zb(i5,i4))/
     &   (za(i5,i6)*zb(i2,i1)*zb(i3,i2)*zb(i4,i3)) -
     &  (xmass*za(i3,i6)*zb(i6,i4))/(za(i5,i6)*zb(i2,i1)*zb(i3,i2)) -
     &  (xmass*za(i1,i6)*zb(i4,i1)*zb(i6,i4))/
     &   (za(i5,i6)*zb(i2,i1)*zb(i3,i2)*zb(i4,i3)) -
     &  (xmass*za(i2,i6)*zb(i4,i2)*zb(i6,i4))/
     &   (za(i5,i6)*zb(i2,i1)*zb(i3,i2)*zb(i4,i3))
      amp(1,1,2,2,2)=-((xmass*za(i1,i5)*za(i2,i3)*zb(i3,i1)*zb(i5,i3))/
     &     (s(i2,i3)*za(i3,i4)*za(i5,i6)*zb(i2,i1))) -
     &  (xmass*za(i2,i3)*za(i2,i5)*zb(i3,i2)*zb(i5,i3))/
     &   (s(i2,i3)*za(i3,i4)*za(i5,i6)*zb(i2,i1)) +
     &  (xmass*za(i1,i5)*za(i2,i3)*za(i2,i4)*zb(i4,i3)*zb(i5,i3))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i3,i4)*za(i5,i6))
     &    - (xmass*za(i1,i5)*za(i2,i4)*zb(i3,i1)*zb(i5,i4))/
     &   (s(i2,i3)*za(i3,i4)*za(i5,i6)*zb(i2,i1)) +
     &  (xmass*za(i1,i2)*za(i1,i5)*zb(i3,i1)**2*zb(i5,i4))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i5,i6)*zb(i2,i1))
     &    - (xmass*za(i2,i4)*za(i2,i5)*zb(i3,i2)*zb(i5,i4))/
     &   (s(i2,i3)*za(i3,i4)*za(i5,i6)*zb(i2,i1)) +
     &  (xmass*za(i1,i2)*za(i2,i5)*zb(i3,i1)*zb(i3,i2)*zb(i5,i4))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i5,i6)*zb(i2,i1))
     &    + (xmass*za(i1,i5)*za(i2,i4)**2*zb(i4,i3)*zb(i5,i4))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i3,i4)*za(i5,i6))
     &    + (xmass*za(i1,i6)*za(i2,i3)*zb(i3,i1)*zb(i6,i3))/
     &   (s(i2,i3)*za(i3,i4)*za(i5,i6)*zb(i2,i1)) +
     &  (xmass*za(i2,i3)*za(i2,i6)*zb(i3,i2)*zb(i6,i3))/
     &   (s(i2,i3)*za(i3,i4)*za(i5,i6)*zb(i2,i1)) -
     &  (xmass*za(i1,i6)*za(i2,i3)*za(i2,i4)*zb(i4,i3)*zb(i6,i3))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i3,i4)*za(i5,i6))
     &    + (xmass*za(i1,i6)*za(i2,i4)*zb(i3,i1)*zb(i6,i4))/
     &   (s(i2,i3)*za(i3,i4)*za(i5,i6)*zb(i2,i1)) -
     &  (xmass*za(i1,i2)*za(i1,i6)*zb(i3,i1)**2*zb(i6,i4))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i5,i6)*zb(i2,i1))
     &    + (xmass*za(i2,i4)*za(i2,i6)*zb(i3,i2)*zb(i6,i4))/
     &   (s(i2,i3)*za(i3,i4)*za(i5,i6)*zb(i2,i1)) -
     &  (xmass*za(i1,i2)*za(i2,i6)*zb(i3,i1)*zb(i3,i2)*zb(i6,i4))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i5,i6)*zb(i2,i1))
     &    - (xmass*za(i1,i6)*za(i2,i4)**2*zb(i4,i3)*zb(i6,i4))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i3,i4)*za(i5,i6))
      amp(1,2,1,2,2)=-((xmass*za(i1,i5)*za(i2,i3)*zb(i4,i2)**2*
     &     zb(i5,i2))/
     &     (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6)*
     &       zb(i4,i3))) + (xmass*za(i1,i3)**2*za(i1,i5)*zb(i2,i1)*
     &     zb(i5,i4))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i1,i2)*za(i5,i6))
     &    - (xmass*za(i1,i3)**2*za(i3,i5)*zb(i3,i2)*zb(i5,i4))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i1,i2)*za(i5,i6))
     &    - (xmass*za(i1,i3)*za(i1,i5)*zb(i4,i2)*zb(i5,i4))/
     &   (s(i2,i3)*za(i1,i2)*za(i5,i6)*zb(i4,i3)) +
     &  (xmass*za(i1,i5)*za(i3,i4)*zb(i4,i2)**2*zb(i5,i4))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6)*zb(i4,i3))
     &    + (xmass*za(i1,i6)*za(i2,i3)*zb(i4,i2)**2*zb(i6,i2))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6)*zb(i4,i3))
     &    - (xmass*za(i1,i3)**2*za(i1,i6)*zb(i2,i1)*zb(i6,i4))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i1,i2)*za(i5,i6))
     &    + (xmass*za(i1,i3)**2*za(i3,i6)*zb(i3,i2)*zb(i6,i4))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i1,i2)*za(i5,i6))
     &    + (xmass*za(i1,i3)*za(i1,i6)*zb(i4,i2)*zb(i6,i4))/
     &   (s(i2,i3)*za(i1,i2)*za(i5,i6)*zb(i4,i3)) -
     &  (xmass*za(i1,i6)*za(i3,i4)*zb(i4,i2)**2*zb(i6,i4))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6)*zb(i4,i3))
      amp(2,1,1,2,2)= (xmass*za(i2,i5)*zb(i5,i1))/(za(i5,i6)*zb(i3,i2)*
     &     zb(i4,i3)) +
     &  (xmass*za(i3,i5)*zb(i3,i1)*zb(i5,i1))/
     &   (za(i5,i6)*zb(i2,i1)*zb(i3,i2)*zb(i4,i3)) +
     &  (xmass*za(i4,i5)*zb(i4,i1)*zb(i5,i1))/
     &   (za(i5,i6)*zb(i2,i1)*zb(i3,i2)*zb(i4,i3)) -
     &  (xmass*za(i2,i6)*zb(i6,i1))/(za(i5,i6)*zb(i3,i2)*zb(i4,i3)) -
     &  (xmass*za(i3,i6)*zb(i3,i1)*zb(i6,i1))/
     &   (za(i5,i6)*zb(i2,i1)*zb(i3,i2)*zb(i4,i3)) -
     &  (xmass*za(i4,i6)*zb(i4,i1)*zb(i6,i1))/
     &   (za(i5,i6)*zb(i2,i1)*zb(i3,i2)*zb(i4,i3))
      amp(1,2,2,2,2)=  (xmass*za(i1,i5)*zb(i5,i2))/
     &   (za(i2,i3)*za(i3,i4)*za(i5,i6)) +
     &  (xmass*za(i1,i3)*za(i1,i5)*zb(i5,i3))/
     &   (za(i1,i2)*za(i2,i3)*za(i3,i4)*za(i5,i6)) +
     &  (xmass*za(i1,i4)*za(i1,i5)*zb(i5,i4))/
     &   (za(i1,i2)*za(i2,i3)*za(i3,i4)*za(i5,i6)) -
     &  (xmass*za(i1,i6)*zb(i6,i2))/
     &   (za(i2,i3)*za(i3,i4)*za(i5,i6)) -
     &  (xmass*za(i1,i3)*za(i1,i6)*zb(i6,i3))/
     &   (za(i1,i2)*za(i2,i3)*za(i3,i4)*za(i5,i6)) -
     &  (xmass*za(i1,i4)*za(i1,i6)*zb(i6,i4))/
     &   (za(i1,i2)*za(i2,i3)*za(i3,i4)*za(i5,i6))
      amp(2,1,2,2,2)=-((xmass*za(i2,i4)*za(i4,i5)*zb(i3,i1)*zb(i5,i1))/
     &     (s(i2,i3)*za(i3,i4)*za(i5,i6)*zb(i2,i1))) +
     &  (xmass*za(i1,i2)*za(i4,i5)*zb(i3,i1)**2*zb(i5,i1))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i5,i6)*zb(i2,i1))
     &    - (xmass*za(i2,i4)**2*za(i2,i5)*zb(i3,i2)*zb(i5,i1))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i3,i4)*za(i5,i6))
     &    + (xmass*za(i2,i4)**2*za(i4,i5)*zb(i4,i3)*zb(i5,i1))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i3,i4)*za(i5,i6))
     &    - (xmass*za(i2,i3)*za(i4,i5)*zb(i3,i1)**2*zb(i5,i3))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i5,i6)*zb(i2,i1))
     &    + (xmass*za(i2,i4)*za(i4,i6)*zb(i3,i1)*zb(i6,i1))/
     &   (s(i2,i3)*za(i3,i4)*za(i5,i6)*zb(i2,i1)) -
     &  (xmass*za(i1,i2)*za(i4,i6)*zb(i3,i1)**2*zb(i6,i1))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i5,i6)*zb(i2,i1))
     &    + (xmass*za(i2,i4)**2*za(i2,i6)*zb(i3,i2)*zb(i6,i1))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i3,i4)*za(i5,i6))
     &    - (xmass*za(i2,i4)**2*za(i4,i6)*zb(i4,i3)*zb(i6,i1))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i3,i4)*za(i5,i6))
     &    + (xmass*za(i2,i3)*za(i4,i6)*zb(i3,i1)**2*zb(i6,i3))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i5,i6)*zb(i2,i1))
      amp(2,2,1,2,2)=(xmass*za(i1,i3)**2*za(i4,i5)*zb(i2,i1)*zb(i5,i1))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i1,i2)*za(i5,i6))
     &    - (xmass*za(i1,i3)*za(i3,i5)*zb(i3,i2)*zb(i5,i1))/
     &   (s(i2,i3)*za(i1,i2)*za(i5,i6)*zb(i4,i3)) -
     &  (xmass*za(i1,i3)*za(i4,i5)*zb(i4,i2)*zb(i5,i1))/
     &   (s(i2,i3)*za(i1,i2)*za(i5,i6)*zb(i4,i3)) +
     &  (xmass*za(i3,i4)*za(i3,i5)*zb(i3,i2)*zb(i4,i2)*zb(i5,i1))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6)*zb(i4,i3))
     &    + (xmass*za(i3,i4)*za(i4,i5)*zb(i4,i2)**2*zb(i5,i1))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6)*zb(i4,i3))
     &    + (xmass*za(i1,i3)*za(i2,i3)*za(i4,i5)*zb(i2,i1)*zb(i5,i2))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i1,i2)*za(i5,i6))
     &    - (xmass*za(i2,i3)*za(i3,i5)*zb(i3,i2)*zb(i5,i2))/
     &   (s(i2,i3)*za(i1,i2)*za(i5,i6)*zb(i4,i3)) -
     &  (xmass*za(i2,i3)*za(i4,i5)*zb(i4,i2)*zb(i5,i2))/
     &   (s(i2,i3)*za(i1,i2)*za(i5,i6)*zb(i4,i3)) -
     &  (xmass*za(i1,i3)**2*za(i4,i6)*zb(i2,i1)*zb(i6,i1))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i1,i2)*za(i5,i6))
     &    + (xmass*za(i1,i3)*za(i3,i6)*zb(i3,i2)*zb(i6,i1))/
     &   (s(i2,i3)*za(i1,i2)*za(i5,i6)*zb(i4,i3)) +
     &  (xmass*za(i1,i3)*za(i4,i6)*zb(i4,i2)*zb(i6,i1))/
     &   (s(i2,i3)*za(i1,i2)*za(i5,i6)*zb(i4,i3)) -
     &  (xmass*za(i3,i4)*za(i3,i6)*zb(i3,i2)*zb(i4,i2)*zb(i6,i1))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6)*zb(i4,i3))
     &    - (xmass*za(i3,i4)*za(i4,i6)*zb(i4,i2)**2*zb(i6,i1))/
     &   (s(i2,i3)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6)*zb(i4,i3))
     &    - (xmass*za(i1,i3)*za(i2,i3)*za(i4,i6)*zb(i2,i1)*zb(i6,i2))/
     &   (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i1,i2)*za(i5,i6))
     &    + (xmass*za(i2,i3)*za(i3,i6)*zb(i3,i2)*zb(i6,i2))/
     &   (s(i2,i3)*za(i1,i2)*za(i5,i6)*zb(i4,i3)) +
     &  (xmass*za(i2,i3)*za(i4,i6)*zb(i4,i2)*zb(i6,i2))/
     &   (s(i2,i3)*za(i1,i2)*za(i5,i6)*zb(i4,i3))
      amp(2,2,2,2,2)=(xmass*za(i1,i4)*za(i4,i5)*zb(i5,i1))/
     &   (za(i1,i2)*za(i2,i3)*za(i3,i4)*za(i5,i6)) +
     &  (xmass*za(i2,i4)*za(i4,i5)*zb(i5,i2))/
     &   (za(i1,i2)*za(i2,i3)*za(i3,i4)*za(i5,i6)) +
     &  (xmass*za(i4,i5)*zb(i5,i3))/(za(i1,i2)*za(i2,i3)*za(i5,i6)) -
     &  (xmass*za(i1,i4)*za(i4,i6)*zb(i6,i1))/
     &   (za(i1,i2)*za(i2,i3)*za(i3,i4)*za(i5,i6)) -
     &  (xmass*za(i2,i4)*za(i4,i6)*zb(i6,i2))/
     &   (za(i1,i2)*za(i2,i3)*za(i3,i4)*za(i5,i6)) -
     &  (xmass*za(i4,i6)*zb(i6,i3))/(za(i1,i2)*za(i2,i3)*za(i5,i6))

      return
      end



