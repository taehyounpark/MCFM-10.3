!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_dm_qqb(p,i1,i2,i5,i6,amp_a,amp_b)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'dm_params.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      include 'first.f'
c------ amplitudes=> ME's for q(i1)+Qb(i5)+Q(i6)+qb(i2)+x(3)+x~(4)
      real(dp):: p(mxpart,4),q(mxpart,4)
      integer:: i1,i2,i5,i6
c      real(dp):: msqA,msqB
      complex(dp):: amp_a(2,2,2,2),amp_b(2,2,2,2)
      integer:: h1,h2,h3,h4
c      complex(dp):: aqqb_zbb_new,cor
      complex(dp):: qqb_dm_qqb_VLR,qqb_dm_qqb_VRL,qqb_dm_qqb_VLL
     & ,qqb_dm_qqb_VRR
      complex(dp):: qqb_dm_qqb_Ax_VLR,qqb_dm_qqb_Ax_VRL
     & ,qqb_dm_qqb_Ax_VLL,qqb_dm_qqb_Ax_VRR
      complex(dp):: amp_pa(2,2),amp_pb(2,2),amp_dec(2,2)
      real(dp):: bp,beta,s34

c      complex(dp):: tc1(2,2,2,2),Rc2(2,2,2,2)

      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
                  amp_a(h1,h2,h3,h4)=czip
                  amp_b(h1,h2,h3,h4)=czip
               enddo
            enddo
         enddo
      enddo

c      ampA(:,:,:,:)=czip

      if(xmass>1d-8) then
c---------generate massless phase space
         call gen_masslessvecs(p,q,3,4)
c---------generate spinors
         call spinoru(6,q,za,zb)
      else
c--------massless dm can use usual spinoru
         call spinoru(6,p,za,zb)
      endif

      do h1=1,6
         do h2=1,6
            s(h1,h2)=Dble(za(h1,h2)*zb(h2,h1))
         enddo
      enddo



      if(dm_mediator=='vector') then
c----- Z => qqb L L (LR,RL,LL,RR)
         amp_a(1,1,1,2)=qqb_dm_qqb_VLR(i1,i6,i5,i2,3,4)
         amp_a(1,1,2,1)=qqb_dm_qqb_VRL(i1,i6,i5,i2,3,4)
         amp_a(1,1,1,1)=qqb_dm_qqb_VLL(i1,i6,i5,i2,3,4)
         amp_a(1,1,2,2)=qqb_dm_qqb_VRR(i1,i6,i5,i2,3,4)
c-----Z => bb L L (LR,RL,LL,RR)
         amp_b(1,1,1,2)=qqb_dm_qqb_VLR(i6,i1,i2,i5,3,4)
         amp_b(1,1,2,1)=qqb_dm_qqb_VRL(i6,i1,i2,i5,3,4)
         amp_b(1,1,1,1)=qqb_dm_qqb_VLL(i6,i1,i2,i5,3,4)
         amp_b(1,1,2,2)=qqb_dm_qqb_VRR(i6,i1,i2,i5,3,4)
c-----Z => qqb R R (LR,RL,LL,RR)
         amp_a(2,2,1,2)=-qqb_dm_qqb_VLR(i2,i5,i6,i1,3,4)
         amp_a(2,2,2,1)=-qqb_dm_qqb_VRL(i2,i5,i6,i1,3,4)
         amp_a(2,2,1,1)=-qqb_dm_qqb_VLL(i2,i5,i6,i1,3,4)
         amp_a(2,2,2,2)=-qqb_dm_qqb_VRR(i2,i5,i6,i1,3,4)
c-----Z => bb R R (LR,RL,LL,RR)
         amp_b(2,2,1,2)=-qqb_dm_qqb_VLR(i5,i2,i1,i6,3,4)
         amp_b(2,2,2,1)=-qqb_dm_qqb_VRL(i5,i2,i1,i6,3,4)
         amp_b(2,2,1,1)=-qqb_dm_qqb_VLL(i5,i2,i1,i6,3,4)
         amp_b(2,2,2,2)=-qqb_dm_qqb_VRR(i5,i2,i1,i6,3,4)
c-----Z => qqb L R (LR,RL,LL,RR)
         amp_a(1,2,1,2)=qqb_dm_qqb_VLR(i1,i5,i6,i2,3,4)
         amp_a(1,2,2,1)=qqb_dm_qqb_VRL(i1,i5,i6,i2,3,4)
         amp_a(1,2,1,1)=qqb_dm_qqb_VLL(i1,i5,i6,i2,3,4)
         amp_a(1,2,2,2)=qqb_dm_qqb_VRR(i1,i5,i6,i2,3,4)
c-----Z => bb L R (LR,RL,LL,RR)
         amp_b(1,2,1,2)=-qqb_dm_qqb_VLR(i5,i1,i2,i6,3,4)
         amp_b(1,2,2,1)=-qqb_dm_qqb_VRL(i5,i1,i2,i6,3,4)
         amp_b(1,2,1,1)=-qqb_dm_qqb_VLL(i5,i1,i2,i6,3,4)
         amp_b(1,2,2,2)=-qqb_dm_qqb_VRR(i5,i1,i2,i6,3,4)
c-----Z => qqb R L (LR,RL,LL,RR)
         amp_a(2,1,1,2)=-qqb_dm_qqb_VLR(i2,i6,i5,i1,3,4)
         amp_a(2,1,2,1)=-qqb_dm_qqb_VRL(i2,i6,i5,i1,3,4)
         amp_a(2,1,1,1)=-qqb_dm_qqb_VLL(i2,i6,i5,i1,3,4)
         amp_a(2,1,2,2)=-qqb_dm_qqb_VRR(i2,i6,i5,i1,3,4)
c-----Z => bb R L (LR,RL,LL,RR)
         amp_b(2,1,1,2)=qqb_dm_qqb_VLR(i6,i2,i1,i5,3,4)
         amp_b(2,1,2,1)=qqb_dm_qqb_VRL(i6,i2,i1,i5,3,4)
         amp_b(2,1,1,1)=qqb_dm_qqb_VLL(i6,i2,i1,i5,3,4)
         amp_b(2,1,2,2)=qqb_dm_qqb_VRR(i6,i2,i1,i5,3,4)
      elseif(dm_mediator=='axvect') then
         !----- Z => qqb L L (LR,RL,LL,RR)
         amp_a(1,1,1,2)=qqb_dm_qqb_Ax_VLR(i1,i6,i5,i2,3,4)
         amp_a(1,1,2,1)=qqb_dm_qqb_Ax_VRL(i1,i6,i5,i2,3,4)
         amp_a(1,1,1,1)=qqb_dm_qqb_Ax_VLL(i1,i6,i5,i2,3,4)
         amp_a(1,1,2,2)=qqb_dm_qqb_Ax_VRR(i1,i6,i5,i2,3,4)
c-----Z => bb L L (LR,RL,LL,RR)
         amp_b(1,1,1,2)=qqb_dm_qqb_Ax_VLR(i6,i1,i2,i5,3,4)
         amp_b(1,1,2,1)=qqb_dm_qqb_Ax_VRL(i6,i1,i2,i5,3,4)
         amp_b(1,1,1,1)=qqb_dm_qqb_Ax_VLL(i6,i1,i2,i5,3,4)
         amp_b(1,1,2,2)=qqb_dm_qqb_Ax_VRR(i6,i1,i2,i5,3,4)
c-----Z => qqb R R (LR,RL,LL,RR)
         amp_a(2,2,1,2)=-qqb_dm_qqb_Ax_VLR(i2,i5,i6,i1,3,4)
         amp_a(2,2,2,1)=-qqb_dm_qqb_Ax_VRL(i2,i5,i6,i1,3,4)
         amp_a(2,2,1,1)=-qqb_dm_qqb_Ax_VLL(i2,i5,i6,i1,3,4)
         amp_a(2,2,2,2)=-qqb_dm_qqb_Ax_VRR(i2,i5,i6,i1,3,4)
c-----Z => bb R R (LR,RL,LL,RR)
         amp_b(2,2,1,2)=-qqb_dm_qqb_Ax_VLR(i5,i2,i1,i6,3,4)
         amp_b(2,2,2,1)=-qqb_dm_qqb_Ax_VRL(i5,i2,i1,i6,3,4)
         amp_b(2,2,1,1)=-qqb_dm_qqb_Ax_VLL(i5,i2,i1,i6,3,4)
         amp_b(2,2,2,2)=-qqb_dm_qqb_Ax_VRR(i5,i2,i1,i6,3,4)
c-----Z => qqb L R (LR,RL,LL,RR)
         amp_a(1,2,1,2)=qqb_dm_qqb_Ax_VLR(i1,i5,i6,i2,3,4)
         amp_a(1,2,2,1)=qqb_dm_qqb_Ax_VRL(i1,i5,i6,i2,3,4)
         amp_a(1,2,1,1)=qqb_dm_qqb_Ax_VLL(i1,i5,i6,i2,3,4)
         amp_a(1,2,2,2)=qqb_dm_qqb_Ax_VRR(i1,i5,i6,i2,3,4)
c-----Z => bb L R (LR,RL,LL,RR)
         amp_b(1,2,1,2)=-qqb_dm_qqb_Ax_VLR(i5,i1,i2,i6,3,4)
         amp_b(1,2,2,1)=-qqb_dm_qqb_Ax_VRL(i5,i1,i2,i6,3,4)
         amp_b(1,2,1,1)=-qqb_dm_qqb_Ax_VLL(i5,i1,i2,i6,3,4)
         amp_b(1,2,2,2)=-qqb_dm_qqb_Ax_VRR(i5,i1,i2,i6,3,4)
c-----Z => qqb R L (LR,RL,LL,RR)
         amp_a(2,1,1,2)=-qqb_dm_qqb_Ax_VLR(i2,i6,i5,i1,3,4)
         amp_a(2,1,2,1)=-qqb_dm_qqb_Ax_VRL(i2,i6,i5,i1,3,4)
         amp_a(2,1,1,1)=-qqb_dm_qqb_Ax_VLL(i2,i6,i5,i1,3,4)
         amp_a(2,1,2,2)=-qqb_dm_qqb_Ax_VRR(i2,i6,i5,i1,3,4)
c-----Z => bb R L (LR,RL,LL,RR)
         amp_b(2,1,1,2)=qqb_dm_qqb_Ax_VLR(i6,i2,i1,i5,3,4)
         amp_b(2,1,2,1)=qqb_dm_qqb_Ax_VRL(i6,i2,i1,i5,3,4)
         amp_b(2,1,1,1)=qqb_dm_qqb_Ax_VLL(i6,i2,i1,i5,3,4)
         amp_b(2,1,2,2)=qqb_dm_qqb_Ax_VRR(i6,i2,i1,i5,3,4)
      elseif(dm_mediator=='scalar') then
         call qqb_dm_qqb_Samps(p,i1,i5,i6,i2,za,zb,amp_pa)
         call qqb_dm_qqb_Samps(p,i5,i1,i2,i6,za,zb,amp_pb)
         s34=s(3,4)
         beta=sqrt(1d0-4d0*xmass**2/s34)
         bp=0.5d0*(one+beta)
         call dm_scal_decay(3,4,za,zb,bp,amp_dec)
         do h1=1,2
            do h2=1,2
               do h3=1,2
                  do h4=1,2
            amp_a(h1,h2,h3,h4)=amp_pa(h1,h2)*amp_dec(h3,h4)
            amp_b(h1,h2,h3,h4)=amp_pb(h1,h2)*amp_dec(h3,h4)
         enddo
      enddo
      enddo
      enddo
      elseif(dm_mediator=='pseudo') then
         if(first) then
            first = .false.
            !call check_dmAxC
            error stop
         endif

         call qqb_dm_qqb_Samps(p,i1,i5,i6,i2,za,zb,amp_pa)
         call qqb_dm_qqb_Samps(p,i5,i1,i2,i6,za,zb,amp_pb)
         s34=s(3,4)
         beta=sqrt(1d0-4d0*xmass**2/s34)
         bp=0.5d0*(one+beta)

c         write(6,*) 'In QQ ',3,4,bp,s34,xmass
         call dm_Pscal_decay(3,4,za,zb,bp,amp_dec)
         do h1=1,2
            do h2=1,2
               do h3=1,2
                  do h4=1,2
            amp_a(h1,h2,h3,h4)=amp_pa(h1,h2)*amp_dec(h3,h4)
            amp_b(h1,h2,h3,h4)=amp_pb(h1,h2)*amp_dec(h3,h4)
         enddo
      enddo
      enddo
      enddo
      endif
      return
      end


c--------- function for helicity conserving amplitudes
      function qqb_dm_qqb_VLR(i1,i2,i3,i4,i5,i6)
      implicit none
      include 'types.f'
      complex(dp):: qqb_dm_qqb_VLR

      include 'mxpart.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6
c------ copy exisitng MCFM structure
c--- This corresponds to A++(1,2,3,4) of eq. (12.3) in BDK
c    The notation of BDK calculates the following amplitude

c     q3(L)----<----------q2            q3(L)------<--------q2
c                 0                             0
c                 0                             0
c                 0                             0
c     q1(R)------<--------q4            q1(R)------<--------q4
c             )                                         )
c            (                                         (
c             )                                         )
c     l5(L)-------<-------l6            l5(L)-------<-------l6

c     Note that this function has the property
c     Conjg(aqqb_zbb_new(i1,i2,i3,i4,i5,i6))=
c          -aqqb_zbb_new(i4,i3,i2,i1,i6,i5)
c------ note that the default expresion (i.e. 123456) corresponds
c------ to amp_A(1,2,1,2) and amp_B(2,1,1,2)

      qqb_dm_qqb_VLR=-((-((s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i1,i3)*
     &     za(i4,i5)*
     &         zb(i2,i1)*zb(i6,i1)) +
     &      (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*za(i3,i5)*
     &       zb(i3,i2)*zb(i6,i1) +
     &      (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*za(i4,i5)*
     &       zb(i4,i2)*zb(i6,i1) -
     &      (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i2,i3)*za(i4,i5)*
     &       zb(i2,i1)*zb(i6,i2))/
     &    (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     &      (s(i2,i3) + s(i2,i4) + s(i3,i4))))

      return
      end

      function qqb_dm_qqb_VRL(i1,i2,i3,i4,i5,i6)
      implicit none
      include 'types.f'
      complex(dp):: qqb_dm_qqb_VRL

      include 'mxpart.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6
c------ copy exisitng MCFM structure
c--- This corresponds to A++(1,2,3,4) of eq. (12.3) in BDK
c    The notation of BDK calculates the following amplitude

c     q3(L)----<----------q2            q3(L)------<--------q2
c                 0                             0
c                 0                             0
c                 0                             0
c     q1(R)------<--------q4            q1(R)------<--------q4
c             )                                         )
c            (                                         (
c             )                                         )
c     l5(R)-------<-------l6            l5(R)-------<-------l6

c     Note that this function has the property
c     Conjg(aqqb_zbb_new(i1,i2,i3,i4,i5,i6))=
c          -aqqb_zbb_new(i4,i3,i2,i1,i6,i5)
c------ note that the default expresion (i.e. 123456) corresponds
c------ to amp_A(1,2,2,1) and amp_B(2,1,2,1)

      qqb_dm_qqb_VRL=-((-((s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i1,i3)
     &     *za(i4,i6)*
     &         zb(i2,i1)*zb(i5,i1)) +
     &      (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*za(i3,i6)*
     &       zb(i3,i2)*zb(i5,i1) +
     &      (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*za(i4,i6)*
     &       zb(i4,i2)*zb(i5,i1) -
     &      (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i2,i3)*za(i4,i6)*
     &       zb(i2,i1)*zb(i5,i2))/
     &    (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     &      (s(i2,i3) + s(i2,i4) + s(i3,i4))))
      return
      end

      function qqb_dm_qqb_VLL(i1,i2,i3,i4,i5,i6)
      implicit none
      include 'types.f'
      complex(dp):: qqb_dm_qqb_VLL

      include 'mxpart.f'
      include 'dm_params.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6
c------ copy exisitng MCFM structure
c--- This corresponds to A++(1,2,3,4) of eq. (12.3) in BDK
c    The notation of BDK calculates the following amplitude

c     q3(L)----<----------q2            q3(L)------<--------q2
c                 0                             0
c                 0                             0
c                 0                             0
c     q1(R)------<--------q4            q1(R)------<--------q4
c             )                                         )
c            (                                         (
c             )                                         )
c     l5(L)-----X-<-------l6(L)            l5(R)-------<-------l6

c     Note that this function has the property
c     Conjg(aqqb_zbb_new(i1,i2,i3,i4,i5,i6))=
c          -aqqb_zbb_new(i4,i3,i2,i1,i6,i5)
c------ note that the default expresion (i.e. 123456) corresponds
c------ to amp_A(1,2,1,1) and amp_B(2,1,1,1)

      qqb_dm_qqb_VLL=(xmass*(-((s(i2,i3) + s(i2,i4) + s(i3,i4))*
     &     za(i1,i3)*za(i4,i5)*
     &         zb(i2,i1)*zb(i5,i1)) +
     &      (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*za(i3,i5)*
     &       zb(i3,i2)*zb(i5,i1) +
     &      (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*za(i4,i5)*
     &       zb(i4,i2)*zb(i5,i1) -
     &      (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i2,i3)*za(i4,i5)*
     &       zb(i2,i1)*zb(i5,i2) +
     &      (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i1,i3)*za(i4,i6)*
     &       zb(i2,i1)*zb(i6,i1) -
     &      (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*za(i3,i6)*
     &       zb(i3,i2)*zb(i6,i1) -
     &      (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*za(i4,i6)*
     &       zb(i4,i2)*zb(i6,i1) +
     &      (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i2,i3)*za(i4,i6)*
     &       zb(i2,i1)*zb(i6,i2)))/
     &  (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     &    (s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i6,i5))
      return
      end

            function qqb_dm_qqb_VRR(i1,i2,i3,i4,i5,i6)
            implicit none
      include 'types.f'
      complex(dp):: qqb_dm_qqb_VRR

      include 'mxpart.f'
      include 'dm_params.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6
c------ copy exisitng MCFM structure
c--- This corresponds to A++(1,2,3,4) of eq. (12.3) in BDK
c    The notation of BDK calculates the following amplitude

c     q3(L)----<----------q2            q3(L)------<--------q2
c                 0                             0
c                 0                             0
c                 0                             0
c     q1(R)------<--------q4            q1(R)------<--------q4
c             )                                         )
c            (                                         (
c             )                                         )
c     l5(R)-----X-<-------l6(R)            l5(R)-------<-------l6

c     Note that this function has the property
c     Conjg(aqqb_zbb_new(i1,i2,i3,i4,i5,i6))=
c          -aqqb_zbb_new(i4,i3,i2,i1,i6,i5)
c------ note that the default expresion (i.e. 123456) corresponds
c------ to amp_A(1,2,2,2) and amp_B(2,1,2,2)

      qqb_dm_qqb_VRR=-((xmass*(-((s(i2,i3) + s(i2,i4) + s(i3,i4))*
     &za(i1,i3)*
     &           za(i4,i5)*zb(i2,i1)*zb(i5,i1)) +
     &        (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*za(i3,i5)*
     &         zb(i3,i2)*zb(i5,i1) +
     &        (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*za(i4,i5)*
     &         zb(i4,i2)*zb(i5,i1) -
     &        (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i2,i3)*za(i4,i5)*
     &         zb(i2,i1)*zb(i5,i2) +
     &        (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i1,i3)*za(i4,i6)*
     &         zb(i2,i1)*zb(i6,i1) -
     &        (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*za(i3,i6)*
     &         zb(i3,i2)*zb(i6,i1) -
     &        (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*za(i4,i6)*
     &         zb(i4,i2)*zb(i6,i1) +
     &        (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i2,i3)*za(i4,i6)*
     &         zb(i2,i1)*zb(i6,i2)))/
     &    (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     &      (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6)))
      return
      end
