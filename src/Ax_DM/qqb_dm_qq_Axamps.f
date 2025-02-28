!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c--------- function for helicity conserving amplitudes
      function qqb_dm_qqb_Ax_VLR(i1,i2,i3,i4,i5,i6)
      implicit none
      include 'types.f'
      complex(dp):: qqb_dm_qqb_Ax_VLR

      include 'constants.f'
      include 'mxpart.f'
      include 'dm_params.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6
      real(dp):: bp,beta,s34

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

c     bp = 1/2(1+beta)
c     beta = sqrt(1-4xmass**2/s34)
      s34=Dble(za(i6,i5)*zb(i5,i6))
      beta=sqrt(1d0-4d0*xmass**2/s34)
      bp=0.5d0*(one+beta)

      qqb_dm_qqb_Ax_VLR=((-cone + 2*bp)*
     &     (-((s(i2,i3) + s(i2,i4) + s(i3,i4))*
     &         za(i1,i3)*za(i4,i6)*zb(i2,i1)*zb(i5,i1)) +
     &      (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*
     &       za(i3,i6)*zb(i3,i2)*zb(i5,i1) +
     &      (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*
     &       za(i4,i6)*zb(i4,i2)*zb(i5,i1) -
     &      (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i2,i3)*
     &       za(i4,i6)*zb(i2,i1)*zb(i5,i2)))/
     &  (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     &    (s(i2,i3) + s(i2,i4) + s(i3,i4)))

      return
      end

      function qqb_dm_qqb_Ax_VRL(i1,i2,i3,i4,i5,i6)
      implicit none
      include 'types.f'
      complex(dp):: qqb_dm_qqb_Ax_VRL
      include 'constants.f'
      include 'mxpart.f'
      include 'dm_params.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6
      real(dp):: bp,beta,s34

c------ copy existing MCFM structure
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

c     bp = 1/2(1+beta)
c     beta = sqrt(1-4xmass**2/s34)
      s34=Dble(za(i6,i5)*zb(i5,i6))
      beta=sqrt(1d0-4d0*xmass**2/s34)
      bp=0.5d0*(one+beta)

      qqb_dm_qqb_Ax_VRL=-(((-cone + 2*bp)*
     &     (-((s(i2,i3) + s(i2,i4) + s(i3,i4))*
     &           za(i1,i3)*za(i4,i5)*zb(i2,i1)*zb(i6,i1)) +
     &        (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*
     &         za(i3,i5)*zb(i3,i2)*zb(i6,i1) +
     &        (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*
     &         za(i4,i5)*zb(i4,i2)*zb(i6,i1) -
     &        (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i2,i3)*
     &         za(i4,i5)*zb(i2,i1)*zb(i6,i2)))/
     &    (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     &      (s(i2,i3) + s(i2,i4) + s(i3,i4))))
      return
      end

      function qqb_dm_qqb_Ax_VLL(i1,i2,i3,i4,i5,i6)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'dm_params.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      complex(dp):: qqb_dm_qqb_Ax_VLL
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


      qqb_dm_qqb_Ax_VLL= (xmass*(-((s(i2,i3) + s(i2,i4) + s(i3,i4))*
     &         za(i1,i3)*za(i4,i5)*zb(i2,i1)*zb(i5,i1)) +
     &      (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*
     &       za(i3,i5)*zb(i3,i2)*zb(i5,i1) +
     &      (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*
     &       za(i4,i5)*zb(i4,i2)*zb(i5,i1) -
     &      (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i2,i3)*
     &       za(i4,i5)*zb(i2,i1)*zb(i5,i2) -
     &      (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i1,i3)*
     &       za(i4,i6)*zb(i2,i1)*zb(i6,i1) +
     &      (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*
     &       za(i3,i6)*zb(i3,i2)*zb(i6,i1) +
     &      (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*
     &       za(i4,i6)*zb(i4,i2)*zb(i6,i1) -
     &      (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i2,i3)*
     &       za(i4,i6)*zb(i2,i1)*zb(i6,i2)))/
     &  (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     &    (s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i6,i5))

      return
      end

      function qqb_dm_qqb_Ax_VRR(i1,i2,i3,i4,i5,i6)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'dm_params.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      complex(dp):: qqb_dm_qqb_Ax_VRR
      integer:: i1,i2,i3,i4,i5,i6

c------ copy existing MCFM structure
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

      qqb_dm_qqb_Ax_VRR=(xmass*(-((s(i2,i3) + s(i2,i4) + s(i3,i4))*
     &         za(i1,i3)*za(i4,i5)*zb(i2,i1)*zb(i5,i1)) +
     &      (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*
     &       za(i3,i5)*zb(i3,i2)*zb(i5,i1) +
     &      (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*
     &       za(i4,i5)*zb(i4,i2)*zb(i5,i1) -
     &      (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i2,i3)*
     &       za(i4,i5)*zb(i2,i1)*zb(i5,i2) -
     &      (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i1,i3)*
     &       za(i4,i6)*zb(i2,i1)*zb(i6,i1) +
     &      (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*
     &       za(i3,i6)*zb(i3,i2)*zb(i6,i1) +
     &      (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*
     &       za(i4,i6)*zb(i4,i2)*zb(i6,i1) -
     &      (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i2,i3)*
     &       za(i4,i6)*zb(i2,i1)*zb(i6,i2)))/
     &  (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     &    (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6))
      return
      end
