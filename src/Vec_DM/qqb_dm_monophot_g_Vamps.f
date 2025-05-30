!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_dm_monophot_g_Vamps(p,i1,i2,i3,i4,i5,i6,amp)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'dm_params.f'
      real(dp):: p(mxpart,4),q(mxpart,4)
      complex(dp):: amp(2,2,2,2,2)
      integer:: i1,i2,i3,i4,i5,i6
      complex(dp):: Vec_monophot_helppC
      complex(dp):: Vec_monophot_helpmC
      complex(dp):: Vec_monophot_helppV1
      complex(dp):: Vec_monophot_helpmV1
      complex(dp):: Vec_monophot_helppV2
      complex(dp):: Vec_monophot_helpmV2

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

      amp(:,:,:,:,:)=czip


c===== helicity conserving same sign gluon/photon
      amp(2,2,2,2,1)=Vec_monophot_helppC(i1,i2,i3,i4,i5,i6,za,zb)
      amp(2,2,2,1,2)=Vec_monophot_helppC(i1,i2,i3,i4,i6,i5,za,zb)
c===== line flip
      amp(1,2,2,2,1)=Vec_monophot_helppC(i3,i2,i1,i4,i5,i6,za,zb)
      amp(1,2,2,1,2)=Vec_monophot_helppC(i3,i2,i1,i4,i6,i5,za,zb)

c===== helicity violaiting same sign gluon/photon
      amp(2,2,2,2,2)=Vec_monophot_helppV1(i1,i2,i3,i4,i5,i6,za,zb)
      amp(1,2,2,2,2)=Vec_monophot_helppV1(i3,i2,i1,i4,i5,i6,za,zb)
      amp(2,2,2,1,1)=Vec_monophot_helppV2(i1,i2,i3,i4,i5,i6,za,zb)
      amp(1,2,2,1,1)=Vec_monophot_helppV2(i3,i2,i1,i4,i5,i6,za,zb)


c====== helicity conserving opposite sign gluon/photon
      amp(1,2,1,2,1)=Vec_monophot_helpmC(i1,i2,i3,i4,i5,i6,za,zb)
      amp(1,2,1,1,2)=Vec_monophot_helpmC(i1,i2,i3,i4,i6,i5,za,zb)
c====== line flip
      amp(2,2,1,2,1)=Vec_monophot_helpmC(i3,i2,i1,i4,i5,i6,za,zb)
      amp(2,2,1,1,2)=Vec_monophot_helpmC(i3,i2,i1,i4,i6,i5,za,zb)

c===== helicity violaiting same sign gluon/photon
      amp(2,2,1,2,2)=Vec_monophot_helpmV1(i3,i2,i1,i4,i5,i6,za,zb)
      amp(1,2,1,2,2)=Vec_monophot_helpmV1(i1,i2,i3,i4,i5,i6,za,zb)
      amp(2,2,1,1,1)=Vec_monophot_helpmV2(i3,i2,i1,i4,i5,i6,za,zb)
      amp(1,2,1,1,1)=Vec_monophot_helpmV2(i1,i2,i3,i4,i5,i6,za,zb)



c======= Conjugate
      amp(1,1,1,1,2)=-conjg(amp(2,2,2,2,1))
      amp(1,1,1,2,1)=-conjg(amp(2,2,2,1,2))
      amp(2,1,1,1,2)=-conjg(amp(1,2,2,2,1))
      amp(2,1,1,2,1)=-conjg(amp(1,2,2,1,2))

      amp(1,1,2,2,1)=-conjg(amp(2,2,1,1,2))
      amp(1,1,2,1,2)=-conjg(amp(2,2,1,2,1))
      amp(2,1,2,2,1)=-conjg(amp(1,2,1,1,2))
      amp(2,1,2,1,2)=-conjg(amp(1,2,1,2,1))

      amp(1,1,1,1,1)=conjg(amp(2,2,2,2,2))
      amp(2,1,1,1,1)=conjg(amp(1,2,2,2,2))
      amp(1,1,1,2,2)=conjg(amp(2,2,2,1,1))
      amp(2,1,1,2,2)=conjg(amp(1,2,2,1,1))

      amp(1,1,2,1,1)=conjg(amp(2,2,1,2,2))
      amp(2,1,2,1,1)=conjg(amp(1,2,1,2,2))
      amp(1,1,2,2,2)=conjg(amp(2,2,1,1,1))
      amp(2,1,2,2,2)=conjg(amp(1,2,1,1,1))

      return
      end


      function Vec_monophot_helppC(i1,i2,i3,i4,i5,i6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: Vec_monophot_helppC
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: i1,i2,i3,i4,i5,i6
c===== helicity amplitude for
c===== q(i1)^+g(i2)^+qb(i3)^-+gamma(i4)^+ * Axial (+,-) Current


      Vec_monophot_helppC=-((za(i1,i3)*za(i3,i6)**2*zb(i6,i5))/
     &    (za(i1,i2)*za(i1,i4)*za(i2,i3)*za(i3,i4)))

      return
      end

      function Vec_monophot_helppV1(i1,i2,i3,i4,i5,i6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: Vec_monophot_helppV1
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'dm_params.f'
      integer:: i1,i2,i3,i4,i5,i6
c===== helicity amplitude for
c===== q(i1)^+g(i2)^+qb(i3)^-+gamma(i4)^+ * Axial (+,+) Current

      Vec_monophot_helppV1=(-2*xmass*za(i1,i3)*za(i3,i5)*za(i3,i6)*
     &    zb(i6,i5))/
     &  (za(i1,i2)*za(i1,i4)*za(i2,i3)*za(i3,i4)*
     &    za(i5,i6))
      return
      end


      function Vec_monophot_helppV2(i1,i2,i3,i4,i5,i6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: Vec_monophot_helppV2
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'dm_params.f'
      integer:: i1,i2,i3,i4,i5,i6
c===== helicity amplitude for
c===== q(i1)^+g(i2)^+qb(i3)^-+gamma(i4)^+ * Axial (-,-) Current

      Vec_monophot_helppV2=(ctwo*xmass*za(i1,i3)*za(i3,i5)*za(i3,i6))/
     &  (za(i1,i2)*za(i1,i4)*za(i2,i3)*za(i3,i4))
      return
      end

      function Vec_monophot_helpmC(i1,i2,i3,i4,i5,i6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: Vec_monophot_helpmC
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'dm_params.f'
      integer:: i1,i2,i3,i4,i5,i6
c===== helicity amplitude for
c===== q(i1)^+g(i2)^+qb(i3)^-+gamma(i4)^+ * Axial (+,+) Current
      real(dp):: bp,beta
      complex(dp):: Ax_monophot_helpmC

      beta=sqrt(1._dp-4._dp*xmass**2/s(i5,i6))
      bp=half*(one+beta)



c====== A terms

      Vec_monophot_helpmC=Ax_monophot_helpmC(i1,i2,i3,i4,i5,i6,za,zb)
      Vec_monophot_helpmC=-Vec_monophot_helpmC/(cone-ctwo*bp)

      return
      end

      function Vec_monophot_helpmV1(i1,i2,i3,i4,i5,i6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: Vec_monophot_helpmV1
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'dm_params.f'
      integer:: i1,i2,i3,i4,i5,i6
c===== helicity amplitude for
c===== q(i1)^+g(i2)^+qb(i3)^-+gamma(i4)^+ * Axial (+,+) Current
      integer:: i,j,k
      real(dp):: st(mxpart,mxpart,mxpart)

      st(:,:,:)=zero

      do i=1,6
         do j=1,6
            do k=1,6
               st(i,j,k)=s(i,j)+s(j,k)+s(i,k)
            enddo
         enddo
      enddo

c====== A terms

      Vec_monophot_helpmV1=
     & (xmass*zb(i3,i2)*(za(i1,i5)*za(i2,i4)*zb(i5,i2) +
     &      za(i1,i5)*za(i3,i4)*zb(i5,i3) -
     &      za(i1,i6)*za(i2,i4)*zb(i6,i2) -
     &      za(i1,i6)*za(i3,i4)*zb(i6,i3)))/
     &  (st(i2,i3,i4)*za(i2,i3)*za(i5,i6)*zb(i4,i3))


      Vec_monophot_helpmV1=Vec_monophot_helpmV1+
     & (xmass*(za(i1,i2)**2*za(i4,i5)*zb(i3,i2)*zb(i5,i2) +
     &      za(i1,i2)*za(i1,i3)*za(i4,i5)*zb(i3,i2)*zb(i5,i3) -
     &      za(i1,i2)**2*za(i4,i6)*zb(i3,i2)*zb(i6,i2) +
     &      za(i1,i2)*za(i1,i4)*za(i5,i6)*zb(i5,i3)*zb(i6,i2) -
     &      za(i1,i2)*za(i1,i3)*za(i4,i6)*zb(i3,i2)*zb(i6,i3) +
     &      za(i1,i2)*za(i1,i4)*za(i5,i6)*zb(i5,i2)*zb(i6,i3) +
     &      2*za(i1,i3)*za(i1,i4)*za(i5,i6)*zb(i5,i3)*zb(i6,i3)))/
     &  (za(i1,i2)*za(i2,i3)*za(i5,i6)*zb(i4,i3)*
     &    (za(i1,i4)*zb(i4,i1) + za(i2,i4)*zb(i4,i2) +
     &      za(i3,i4)*zb(i4,i3)))

c===== B terms
      Vec_monophot_helpmV1=Vec_monophot_helpmV1+
     &  (xmass*za(i1,i4)*(za(i1,i5)*zb(i2,i1)*zb(i5,i3) -
     &      za(i4,i5)*zb(i4,i2)*zb(i5,i3) -
     &      za(i1,i6)*zb(i2,i1)*zb(i6,i3) +
     &      za(i4,i6)*zb(i4,i2)*zb(i6,i3)))/
     &  (st(i1,i2,i4)*za(i1,i2)*za(i5,i6)*zb(i4,i1))


      Vec_monophot_helpmV1=Vec_monophot_helpmV1+
     &(xmass*(-(s(i2,i3)*za(i1,i2)*za(i4,i5)*zb(i5,i2)) -
     &      s(i2,i3)*za(i1,i3)*za(i4,i5)*zb(i5,i3) +
     &      s(i2,i3)*za(i1,i2)*za(i4,i6)*zb(i6,i2) +
     &      2*za(i1,i2)*za(i2,i4)*za(i5,i6)*zb(i5,i2)*zb(i6,i2) +
     &      za(i1,i3)*za(i2,i4)*za(i5,i6)*zb(i5,i3)*zb(i6,i2) +
     &      za(i1,i2)*za(i3,i4)*za(i5,i6)*zb(i5,i3)*zb(i6,i2) +
     &      s(i2,i3)*za(i1,i3)*za(i4,i6)*zb(i6,i3) +
     &      za(i1,i3)*za(i2,i4)*za(i5,i6)*zb(i5,i2)*zb(i6,i3) +
     &      za(i1,i2)*za(i3,i4)*za(i5,i6)*zb(i5,i2)*zb(i6,i3) +
     &      2*za(i1,i3)*za(i3,i4)*za(i5,i6)*zb(i5,i3)*zb(i6,i3)))/
     &  (za(i1,i2)*za(i2,i3)*za(i5,i6)*zb(i4,i1)*
     &    (za(i1,i4)*zb(i4,i1) + za(i2,i4)*zb(i4,i2) +
     &      za(i3,i4)*zb(i4,i3)))

      return
      end
      function Vec_monophot_helpmV2(i1,i2,i3,i4,i5,i6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: Vec_monophot_helpmV2
c===== helicity amplitude for
c====== q(i1)^+g(i2)^+qb(i3)^-+gamma(i4)^+ * Axial (-,-) Current
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: i1,i2,i3,i4,i5,i6
      complex(dp):: Vec_monophot_helpmV1

      Vec_monophot_helpmV2=Vec_monophot_helpmV1(i1,i2,i3,i4,i5,i6,za,zb)
      Vec_monophot_helpmV2=Vec_monophot_helpmV2*za(i5,i6)/zb(i5,i6)

      return
      end
