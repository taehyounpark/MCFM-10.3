!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a7treen(k1,k2,k3,k4,k5,k6,k7,p,n,za,zb,Qid,a7n)
c     The function calculates
c      qbar(p1)+q(p2)-->W^+ + W^- + g(p7)
c                       |     |
c                       |     |--> e^-(p6)+nu(p5)
c                       |
c                       |----> nu(p3)+e^+(p4)
c     with line 7 contracted with n
c     jtype=1 for a-type diagrams
c     jtype=2 for b-type diagrams
c     jtype=3 for singly-resonant diagrams
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'srdiags.f'
      integer i,j,k1,k2,k3,k4,k5,k6,k7,Qid
      real(dp):: n(4),p(mxpart,4),
     & t134,t256,t127,s34,s56
      complex(dp):: a7n(3,2),vecm(mxpart,mxpart)
      complex(dp):: zba2,zba3,iza(7,7),izb(7,7),sra(2),srb(2)
      zba2(k1,k2,k3,k4)=zb(k1,k2)*za(k2,k4)+zb(k1,k3)*za(k3,k4)
      zba3(k1,k2,k3,k4,k5)=
     & zb(k1,k2)*za(k2,k5)+zb(k1,k3)*za(k3,k5)+zb(k1,k4)*za(k4,k5)

      call checkndotp(p,n,k7)

      do i=1,7
      do j=i+1,7
         iza(i,j)=cone/za(i,j)
         izb(i,j)=cone/zb(i,j)
         iza(j,i)=-iza(i,j)
         izb(j,i)=-izb(i,j)
      enddo
      enddo
      s34=s(k3,k4)
      s56=s(k5,k6)
      t127=s(k1,k2)+s(k2,k7)+s(k7,k1)
      t134=s(k1,k3)+s(k3,k4)+s(k4,k1)
      t256=s(k2,k5)+s(k2,k6)+s(k5,k6)
      do i=1,7
      do j=1,7
      call ndveccur(i,j,n,p,vecm)
      enddo
      enddo
      a7n(1,1)= + iza(k1,k7)*izb(k1,k7)*vecm(k1,k1)*t256**(-1)*
     & s34**(-1)*s56**(-1) * (  - zba2(k4,k2,k5,k6)*za(k1,k3)*zb(k2,k5)
     &     )
      a7n(1,1) = a7n(1,1) + iza(k1,k7)*izb(k1,k7)*vecm(k1,k7)*
     & t256**(-1)*s34**(-1)*s56**(-1) * ( zba2(k4,k2,k5,k6)*za(k3,k7)*
     &    zb(k2,k5) )
      a7n(1,1) = a7n(1,1) + iza(k2,k7)*izb(k2,k7)*vecm(k2,k2)*
     & t134**(-1)*s34**(-1)*s56**(-1) * (  - zba2(k4,k1,k3,k6)*za(k1,k3
     &    )*zb(k2,k5) )
      a7n(1,1) = a7n(1,1) + iza(k2,k7)*izb(k2,k7)*vecm(k7,k2)*
     & t134**(-1)*s34**(-1)*s56**(-1) * ( zba2(k4,k1,k3,k6)*za(k1,k3)*
     &    zb(k5,k7) )
      a7n(1,1) = a7n(1,1) + vecm(k1,k2)*t256**(-1)*t134**(-1)*s34**(-1)
     & *s56**(-1) * (  - za(k1,k3)*za(k2,k6)*zb(k1,k4)*zb(k2,k5) )
      a7n(1,1) = a7n(1,1) + vecm(k1,k5)*t256**(-1)*t134**(-1)*s34**(-1)
     & *s56**(-1) * (  - za(k1,k3)*za(k5,k6)*zb(k1,k4)*zb(k2,k5) )
      a7n(1,1) = a7n(1,1) + vecm(k3,k2)*t256**(-1)*t134**(-1)*s34**(-1)
     & *s56**(-1) * (  - za(k1,k3)*za(k2,k6)*zb(k2,k5)*zb(k3,k4) )
      a7n(1,1) = a7n(1,1) + vecm(k3,k5)*t256**(-1)*t134**(-1)*s34**(-1)
     & *s56**(-1) * (  - za(k1,k3)*za(k5,k6)*zb(k2,k5)*zb(k3,k4) )

      a7n(1,2)=czip
      a7n(2,1)= + iza(k1,k7)*izb(k1,k7)*vecm(k1,k1)*s34**(-1)*s56**(-1)
     & *t127**(-1) * ( zba2(k2,k6,k5,k1)*za(k3,k6)*zb(k4,k5) - zba2(k4,
     &    k6,k5,k3)*za(k1,k6)*zb(k2,k5) - zba3(k5,k1,k2,k7,k6)*za(k1,k3
     &    )*zb(k2,k4) )
      a7n(2,1) = a7n(2,1) + iza(k1,k7)*izb(k1,k7)*vecm(k1,k7)*s34**(-1)
     & *s56**(-1)*t127**(-1) * ( zba2(k2,k6,k5,k7)*za(k3,k6)*zb(k4,k5)
     &     + zba2(k4,k6,k5,k3)*za(k6,k7)*zb(k2,k5) + zba3(k5,k1,k2,k7,
     &    k6)*za(k3,k7)*zb(k2,k4) )
      a7n(2,1) = a7n(2,1) + iza(k2,k7)*izb(k2,k7)*vecm(k2,k2)*s34**(-1)
     & *s56**(-1)*t127**(-1) * (  - zba2(k2,k6,k5,k1)*za(k3,k6)*zb(k4,
     &    k5) + zba2(k4,k6,k5,k3)*za(k1,k6)*zb(k2,k5) + zba3(k5,k1,k2,
     &    k7,k6)*za(k1,k3)*zb(k2,k4) )
      a7n(2,1) = a7n(2,1) + iza(k2,k7)*izb(k2,k7)*vecm(k7,k2)*s34**(-1)
     & *s56**(-1)*t127**(-1) * (  - zba2(k4,k6,k5,k3)*za(k1,k6)*zb(k5,
     &    k7) - zba2(k7,k6,k5,k1)*za(k3,k6)*zb(k4,k5) - zba3(k5,k1,k2,
     &    k7,k6)*za(k1,k3)*zb(k4,k7) )

      a7n(2,2)= + iza(k1,k7)*izb(k1,k7)*vecm(k1,k1)*s34**(-1)*s56**(-1)
     & *t127**(-1) * ( zba2(k1,k5,k6,k2)*za(k3,k6)*zb(k4,k5) - zba2(k4,
     &    k6,k5,k3)*za(k2,k6)*zb(k1,k5) - zba3(k5,k1,k2,k7,k6)*za(k2,k3
     &    )*zb(k1,k4) )
      a7n(2,2) = a7n(2,2) + iza(k1,k7)*izb(k1,k7)*vecm(k7,k1)*s34**(-1)
     & *s56**(-1)*t127**(-1) * ( zba2(k4,k6,k5,k3)*za(k2,k6)*zb(k5,k7)
     &     + zba2(k7,k5,k6,k2)*za(k3,k6)*zb(k4,k5) + zba3(k5,k1,k2,k7,
     &    k6)*za(k2,k3)*zb(k4,k7) )
      a7n(2,2) = a7n(2,2) + iza(k2,k7)*izb(k2,k7)*vecm(k2,k2)*s34**(-1)
     & *s56**(-1)*t127**(-1) * (  - zba2(k1,k5,k6,k2)*za(k3,k6)*zb(k4,
     &    k5) + zba2(k4,k6,k5,k3)*za(k2,k6)*zb(k1,k5) + zba3(k5,k1,k2,
     &    k7,k6)*za(k2,k3)*zb(k1,k4) )
      a7n(2,2) = a7n(2,2) + iza(k2,k7)*izb(k2,k7)*vecm(k2,k7)*s34**(-1)
     & *s56**(-1)*t127**(-1) * (  - zba2(k1,k5,k6,k7)*za(k3,k6)*zb(k4,
     &    k5) - zba2(k4,k6,k5,k3)*za(k6,k7)*zb(k1,k5) - zba3(k5,k1,k2,
     &    k7,k6)*za(k3,k7)*zb(k1,k4) )

      if (srdiags) then

        call srn(k1,k2,3,4,5,6,k7,za,zb,vecm,Qid,.false.,sra)
        call srn(k1,k2,6,5,4,3,k7,za,zb,vecm,Qid,.true.,srb)

        a7n(3,1)=-(sra(1)+srb(1))
        a7n(3,2)=-(sra(2)+srb(2))

      else

        a7n(3,:)=czip

      endif

      return
      end
