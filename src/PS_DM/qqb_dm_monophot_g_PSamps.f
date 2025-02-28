!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_dm_monophot_g_PSamps(p,i3,i2,i1,i4,i5,i6,amp)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'dm_params.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      real(dp):: p(mxpart,4),q(mxpart,4)
      complex(dp):: amp(2,2,2,2,2)
      integer:: i1,i2,i3,i4,i5,i6
      integer:: h1,h2,h3,h4,h5
      complex(dp):: amp_p(2,2,2),amp_dec(2,2)
      real(dp):: s34,beta,bp
      complex(dp):: s1234
c=========== q(i1)+g(i2)+qb(i3)+gamma(i4)+x(i5)+x~(i6)
c---- order is quark helicity,gluon helicity photon helciity dm
      if(xmass>1d-8) then
c---------generate massless phase space
         call gen_masslessvecs(p,q,i5,i6)
c--------- generate spinors
         call spinoru(6,q,za,zb)
      else
c--------massless dm can use usual spinoru
         call spinoru(6,p,za,zb)
      endif

      amp_p(:,:,:)=czip
      amp(:,:,:,:,:)=czip
c========== Helicity conserving amplitudes

      do h1=1,6
         do h2=1,6
            s(h1,h2)=real(za(h1,h2)*zb(h2,h1))
         enddo
      enddo



      s34=real(za(i5,i6)*zb(i6,i5))
      beta=sqrt(1d0-4d0*xmass**2/s34)
      bp=0.5d0*(one+beta)

      s1234= za(i1,i2)*zb(i2,i1) + za(i1,i3)*zb(i3,i1) +
     &  za(i2,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i1) +
     &  za(i2,i4)*zb(i4,i2) + za(i3,i4)*zb(i4,i3)
      call dm_Pscal_decay(i5,i6,za,zb,bp,amp_dec)

      amp_p(1,1,1)= (s1234*zb(i3,i1))/
     &  (zb(i2,i1)*zb(i3,i2)*zb(i4,i1)*
     &    zb(i4,i3))
      amp_p(1,1,2)=-(((s(i1,i2) + s(i1,i3) + s(i2,i3))*
     &       za(i1,i3))/
     &     (za(i1,i4)*za(i3,i4)*zb(i2,i1)*
     &       zb(i3,i2))) +
     &  (za(i1,i2)*
     &     (za(i1,i3)*zb(i4,i1) +
     &       za(i2,i3)*zb(i4,i2)))/
     &   ((s(i1,i2) + s(i1,i4) + s(i2,i4))*
     &     za(i1,i4)*zb(i2,i1)) -
     &  (za(i2,i3)*
     &     (-(za(i1,i2)*zb(i4,i2)) -
     &       za(i1,i3)*zb(i4,i3)))/
     &   ((s(i2,i3) + s(i2,i4) + s(i3,i4))*
     &     za(i3,i4)*zb(i3,i2))
      amp_p(1,2,1)= (za(i1,i4)*
     &     (za(i1,i3)*zb(i2,i1) +
     &       za(i3,i4)*zb(i4,i2)))/
     &   ((s(i1,i2) + s(i1,i4) + s(i2,i4))*
     &     za(i1,i2)*zb(i4,i1)) -
     &  ((s(i1,i3) + s(i1,i4) + s(i3,i4))*
     &     za(i1,i3))/
     &   (za(i1,i2)*za(i2,i3)*zb(i4,i1)*
     &     zb(i4,i3)) +
     &  (za(i3,i4)*
     &     (za(i1,i3)*zb(i3,i2) +
     &       za(i1,i4)*zb(i4,i2)))/
     &   ((s(i2,i3) + s(i2,i4) + s(i3,i4))*
     &     za(i2,i3)*zb(i4,i3))
      amp_p(2,1,1)= zb(i3,i1)**3/
     &  (zb(i2,i1)*zb(i3,i2)*zb(i4,i1)*
     &    zb(i4,i3))


      amp_p(1,2,2)=-(za(i1,i3)**3/
     &    (za(i1,i2)*za(i1,i4)*za(i2,i3)*
     &      za(i3,i4)))
      amp_p(2,1,2)= ((s(i1,i3) + s(i1,i4) + s(i3,i4))*
     &     zb(i3,i1))/
     &   (za(i1,i4)*za(i3,i4)*zb(i2,i1)*
     &     zb(i3,i2)) -
     &  ((za(i2,i3)*zb(i3,i1) +
     &       za(i2,i4)*zb(i4,i1))*zb(i4,i3))/
     &   ((s(i2,i3) + s(i2,i4) + s(i3,i4))*
     &     za(i3,i4)*zb(i3,i2)) -
     &  (zb(i4,i1)*
     &     (za(i1,i2)*zb(i3,i1) +
     &       za(i2,i4)*zb(i4,i3)))/
     &   ((s(i1,i2) + s(i1,i4) + s(i2,i4))*
     &     za(i1,i4)*zb(i2,i1))
      amp_p(2,2,1)=-((zb(i2,i1)*
     &       (za(i1,i4)*zb(i3,i1) +
     &         za(i2,i4)*zb(i3,i2)))/
     &     ((s(i1,i2) + s(i1,i4) + s(i2,i4))*
     &       za(i1,i2)*zb(i4,i1))) +
     &  ((-(za(i2,i4)*zb(i2,i1)) -
     &       za(i3,i4)*zb(i3,i1))*zb(i3,i2))/
     &   ((s(i2,i3) + s(i2,i4) + s(i3,i4))*
     &     za(i2,i3)*zb(i4,i3)) +
     &  ((s(i1,i2) + s(i1,i3) + s(i2,i3))*
     &     zb(i3,i1))/
     &   (za(i1,i2)*za(i2,i3)*zb(i4,i1)*
     &     zb(i4,i3))
      amp_p(2,2,2)=-((s1234*za(i1,i3))/(za(i1,i2)*za(i1,i4)*
     &     za(i2,i3)*za(i3,i4)))


      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
                  do h5=1,2
                     amp(h1,h2,h3,h4,h5)=amp_p(h1,h2,h3)*amp_dec(h4,h5)
                  enddo
               enddo
            enddo
         enddo
      enddo

      return
      end
