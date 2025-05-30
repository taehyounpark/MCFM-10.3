!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_dm_qq_vecamps(p,i1,i2,i3,i4,i5,i6,amp)
      implicit none
      include 'types.f'
      include 'dm_params.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      real(dp):: p(mxpart,4)
c----- fills amplitude for q qb Q Qb chi,chib
c-----
      complex(dp):: amp(2,2,2,2)
c------ helicity ordering is sign of q, and Q and chi chib
      real(dp):: q(mxpart,4)
      integer:: i1,i2,i3,i4,i5,i6

       amp(:,:,:,:)=czip
c--------- generate massless phase space
      call gen_masslessvecs(p,q,i5,i6)
c--------- generate spinors
      call spinoru(6,q,za,zb)
      call dotem(6,q,s)

      amp(1,1,1,1)=-((xmass*za(i1,i3)*za(i1,i5)*zb(i4,i1)*zb(i5,i2))/
     &     (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*zb(i6,i5))) +
     &  (xmass*za(i1,i5)*za(i2,i3)*zb(i4,i2)*zb(i5,i2))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i6,i5)) -
     &  (xmass*za(i1,i3)*za(i3,i5)*zb(i4,i3)*zb(i5,i2))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*zb(i6,i5)) -
     &  (xmass*za(i1,i5)*za(i3,i4)*zb(i4,i2)*zb(i5,i4))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i6,i5)) +
     &  (xmass*za(i1,i3)*za(i1,i6)*zb(i4,i1)*zb(i6,i2))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*zb(i6,i5)) -
     &  (xmass*za(i1,i6)*za(i2,i3)*zb(i4,i2)*zb(i6,i2))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i6,i5)) +
     &  (xmass*za(i1,i3)*za(i3,i6)*zb(i4,i3)*zb(i6,i2))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*zb(i6,i5)) +
     &  (xmass*za(i1,i6)*za(i3,i4)*zb(i4,i2)*zb(i6,i4))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i6,i5))

      amp(1,1,1,2)=(za(i1,i3)*za(i1,i5)*zb(i4,i1)*zb(i6,i2))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))) -
     &  (za(i1,i5)*za(i2,i3)*zb(i4,i2)*zb(i6,i2))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))) +
     &  (za(i1,i3)*za(i3,i5)*zb(i4,i3)*zb(i6,i2))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))) +
     &  (za(i1,i5)*za(i3,i4)*zb(i4,i2)*zb(i6,i4))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4)))

      amp(1,1,2,1)=(za(i1,i3)*za(i1,i6)*zb(i4,i1)*zb(i5,i2))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))) -
     &  (za(i1,i6)*za(i2,i3)*zb(i4,i2)*zb(i5,i2))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))) +
     &  (za(i1,i3)*za(i3,i6)*zb(i4,i3)*zb(i5,i2))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))) +
     &  (za(i1,i6)*za(i3,i4)*zb(i4,i2)*zb(i5,i4))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4)))

      amp(1,1,2,2)=(xmass*za(i1,i3)*za(i1,i5)*zb(i4,i1)*zb(i5,i2))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*za(i5,i6)) -
     &  (xmass*za(i1,i5)*za(i2,i3)*zb(i4,i2)*zb(i5,i2))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6)) +
     &  (xmass*za(i1,i3)*za(i3,i5)*zb(i4,i3)*zb(i5,i2))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*za(i5,i6)) +
     &  (xmass*za(i1,i5)*za(i3,i4)*zb(i4,i2)*zb(i5,i4))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6)) -
     &  (xmass*za(i1,i3)*za(i1,i6)*zb(i4,i1)*zb(i6,i2))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*za(i5,i6)) +
     &  (xmass*za(i1,i6)*za(i2,i3)*zb(i4,i2)*zb(i6,i2))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6)) -
     &  (xmass*za(i1,i3)*za(i3,i6)*zb(i4,i3)*zb(i6,i2))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*za(i5,i6)) -
     &  (xmass*za(i1,i6)*za(i3,i4)*zb(i4,i2)*zb(i6,i4))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6))

      amp(1,2,1,1)=-((xmass*za(i1,i4)*za(i1,i5)*zb(i3,i1)*zb(i5,i2))/
     &     (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*zb(i6,i5))) +
     &  (xmass*za(i1,i5)*za(i2,i4)*zb(i3,i2)*zb(i5,i2))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i6,i5)) +
     &  (xmass*za(i1,i4)*za(i4,i5)*zb(i4,i3)*zb(i5,i2))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*zb(i6,i5)) +
     &  (xmass*za(i1,i5)*za(i3,i4)*zb(i3,i2)*zb(i5,i3))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i6,i5)) +
     &  (xmass*za(i1,i4)*za(i1,i6)*zb(i3,i1)*zb(i6,i2))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*zb(i6,i5)) -
     &  (xmass*za(i1,i6)*za(i2,i4)*zb(i3,i2)*zb(i6,i2))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i6,i5)) -
     &  (xmass*za(i1,i4)*za(i4,i6)*zb(i4,i3)*zb(i6,i2))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*zb(i6,i5)) -
     &  (xmass*za(i1,i6)*za(i3,i4)*zb(i3,i2)*zb(i6,i3))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i6,i5))

      amp(1,2,1,2)=(za(i1,i4)*za(i1,i5)*zb(i3,i1)*zb(i6,i2))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))) -
     &  (za(i1,i5)*za(i2,i4)*zb(i3,i2)*zb(i6,i2))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))) -
     &  (za(i1,i4)*za(i4,i5)*zb(i4,i3)*zb(i6,i2))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))) -
     &  (za(i1,i5)*za(i3,i4)*zb(i3,i2)*zb(i6,i3))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4)))

      amp(1,2,2,1)=(za(i1,i4)*za(i1,i6)*zb(i3,i1)*zb(i5,i2))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))) -
     &  (za(i1,i6)*za(i2,i4)*zb(i3,i2)*zb(i5,i2))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))) -
     &  (za(i1,i4)*za(i4,i6)*zb(i4,i3)*zb(i5,i2))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))) -
     &  (za(i1,i6)*za(i3,i4)*zb(i3,i2)*zb(i5,i3))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4)))

      amp(1,2,2,2)=(xmass*za(i1,i4)*za(i1,i5)*zb(i3,i1)*zb(i5,i2))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*za(i5,i6)) -
     &  (xmass*za(i1,i5)*za(i2,i4)*zb(i3,i2)*zb(i5,i2))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6)) -
     &  (xmass*za(i1,i4)*za(i4,i5)*zb(i4,i3)*zb(i5,i2))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*za(i5,i6)) -
     &  (xmass*za(i1,i5)*za(i3,i4)*zb(i3,i2)*zb(i5,i3))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6)) -
     &  (xmass*za(i1,i4)*za(i1,i6)*zb(i3,i1)*zb(i6,i2))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*za(i5,i6)) +
     &  (xmass*za(i1,i6)*za(i2,i4)*zb(i3,i2)*zb(i6,i2))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6)) +
     &  (xmass*za(i1,i4)*za(i4,i6)*zb(i4,i3)*zb(i6,i2))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*za(i5,i6)) +
     &  (xmass*za(i1,i6)*za(i3,i4)*zb(i3,i2)*zb(i6,i3))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6))

      amp(2,1,1,1)=(xmass*za(i1,i3)*za(i2,i5)*zb(i4,i1)*zb(i5,i1))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*zb(i6,i5)) -
     &  (xmass*za(i2,i3)*za(i2,i5)*zb(i4,i2)*zb(i5,i1))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i6,i5)) -
     &  (xmass*za(i2,i3)*za(i3,i5)*zb(i4,i3)*zb(i5,i1))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i6,i5)) -
     &  (xmass*za(i2,i5)*za(i3,i4)*zb(i4,i1)*zb(i5,i4))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*zb(i6,i5)) -
     &  (xmass*za(i1,i3)*za(i2,i6)*zb(i4,i1)*zb(i6,i1))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*zb(i6,i5)) +
     &  (xmass*za(i2,i3)*za(i2,i6)*zb(i4,i2)*zb(i6,i1))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i6,i5)) +
     &  (xmass*za(i2,i3)*za(i3,i6)*zb(i4,i3)*zb(i6,i1))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i6,i5)) +
     &  (xmass*za(i2,i6)*za(i3,i4)*zb(i4,i1)*zb(i6,i4))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*zb(i6,i5))

      amp(2,1,1,2)=-((za(i1,i3)*za(i2,i5)*zb(i4,i1)*zb(i6,i1))/
     &     (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4)))) +
     &  (za(i2,i3)*za(i2,i5)*zb(i4,i2)*zb(i6,i1))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))) +
     &  (za(i2,i3)*za(i3,i5)*zb(i4,i3)*zb(i6,i1))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))) +
     &  (za(i2,i5)*za(i3,i4)*zb(i4,i1)*zb(i6,i4))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4)))

      amp(2,1,2,1)=-((za(i1,i3)*za(i2,i6)*zb(i4,i1)*zb(i5,i1))/
     &     (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4)))) +
     &  (za(i2,i3)*za(i2,i6)*zb(i4,i2)*zb(i5,i1))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))) +
     &  (za(i2,i3)*za(i3,i6)*zb(i4,i3)*zb(i5,i1))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))) +
     &  (za(i2,i6)*za(i3,i4)*zb(i4,i1)*zb(i5,i4))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4)))

      amp(2,1,2,2)=-((xmass*za(i1,i3)*za(i2,i5)*zb(i4,i1)*zb(i5,i1))/
     &     (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*za(i5,i6))) +
     &  (xmass*za(i2,i3)*za(i2,i5)*zb(i4,i2)*zb(i5,i1))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6)) +
     &  (xmass*za(i2,i3)*za(i3,i5)*zb(i4,i3)*zb(i5,i1))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6)) +
     &  (xmass*za(i2,i5)*za(i3,i4)*zb(i4,i1)*zb(i5,i4))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*za(i5,i6)) +
     &  (xmass*za(i1,i3)*za(i2,i6)*zb(i4,i1)*zb(i6,i1))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*za(i5,i6)) -
     &  (xmass*za(i2,i3)*za(i2,i6)*zb(i4,i2)*zb(i6,i1))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6)) -
     &  (xmass*za(i2,i3)*za(i3,i6)*zb(i4,i3)*zb(i6,i1))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6)) -
     &  (xmass*za(i2,i6)*za(i3,i4)*zb(i4,i1)*zb(i6,i4))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*za(i5,i6))

      amp(2,2,1,1)=(xmass*za(i1,i4)*za(i2,i5)*zb(i3,i1)*zb(i5,i1))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*zb(i6,i5)) -
     &  (xmass*za(i2,i4)*za(i2,i5)*zb(i3,i2)*zb(i5,i1))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i6,i5)) +
     &  (xmass*za(i2,i4)*za(i4,i5)*zb(i4,i3)*zb(i5,i1))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i6,i5)) +
     &  (xmass*za(i2,i5)*za(i3,i4)*zb(i3,i1)*zb(i5,i3))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*zb(i6,i5)) -
     &  (xmass*za(i1,i4)*za(i2,i6)*zb(i3,i1)*zb(i6,i1))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*zb(i6,i5)) +
     &  (xmass*za(i2,i4)*za(i2,i6)*zb(i3,i2)*zb(i6,i1))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i6,i5)) -
     &  (xmass*za(i2,i4)*za(i4,i6)*zb(i4,i3)*zb(i6,i1))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i6,i5)) -
     &  (xmass*za(i2,i6)*za(i3,i4)*zb(i3,i1)*zb(i6,i3))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*zb(i6,i5))

      amp(2,2,1,2)=-((za(i1,i4)*za(i2,i5)*zb(i3,i1)*zb(i6,i1))/
     &     (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4)))) +
     &  (za(i2,i4)*za(i2,i5)*zb(i3,i2)*zb(i6,i1))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))) -
     &  (za(i2,i4)*za(i4,i5)*zb(i4,i3)*zb(i6,i1))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))) -
     &  (za(i2,i5)*za(i3,i4)*zb(i3,i1)*zb(i6,i3))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4)))

      amp(2,2,2,1)=-((za(i1,i4)*za(i2,i6)*zb(i3,i1)*zb(i5,i1))/
     &     (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4)))) +
     &  (za(i2,i4)*za(i2,i6)*zb(i3,i2)*zb(i5,i1))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))) -
     &  (za(i2,i4)*za(i4,i6)*zb(i4,i3)*zb(i5,i1))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))) -
     &  (za(i2,i6)*za(i3,i4)*zb(i3,i1)*zb(i5,i3))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4)))

      amp(2,2,2,2)=-((xmass*za(i1,i4)*za(i2,i5)*zb(i3,i1)*zb(i5,i1))/
     &     (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*za(i5,i6))) +
     &  (xmass*za(i2,i4)*za(i2,i5)*zb(i3,i2)*zb(i5,i1))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6)) -
     &  (xmass*za(i2,i4)*za(i4,i5)*zb(i4,i3)*zb(i5,i1))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6)) -
     &  (xmass*za(i2,i5)*za(i3,i4)*zb(i3,i1)*zb(i5,i3))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*za(i5,i6)) +
     &  (xmass*za(i1,i4)*za(i2,i6)*zb(i3,i1)*zb(i6,i1))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*za(i5,i6)) -
     &  (xmass*za(i2,i4)*za(i2,i6)*zb(i3,i2)*zb(i6,i1))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6)) +
     &  (xmass*za(i2,i4)*za(i4,i6)*zb(i4,i3)*zb(i6,i1))/
     &   (s(i3,i4)*(s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6)) +
     &  (xmass*za(i2,i6)*za(i3,i4)*zb(i3,i1)*zb(i6,i3))/
     &   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*za(i5,i6))


      return
      end
