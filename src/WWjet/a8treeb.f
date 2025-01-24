!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a8treeb(k1,k2,k3,k4,k5,k6,k7,k8,za,zb,a8b)
      implicit none
c     (b)-type amplitude calculated in DKS notation
c     0->u(q1)+ubar(q2)+l(q3)+a(q4)+a'(q5)+l'(q6)+g(p7)+g(p8)
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      real(dp):: s34,s56,s78,t3456,t178,t278
      complex(dp):: a8b(2,2),zbaba222,zbab22,zbab32
      complex(dp):: iza(8,8),izb(8,8),zab2,zba2,zba3
      complex(dp):: zaba22,zaba32,zab3
      integer i,j,k1,k2,k3,k4,k5,k6,k7,k8
c---statement functions
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
      zba2(k1,k2,k3,k4)=zb(k1,k2)*za(k2,k4)+zb(k1,k3)*za(k3,k4)
      zba3(k1,k2,k3,k4,k5)=
     & zb(k1,k2)*za(k2,k5)+zb(k1,k3)*za(k3,k5)+zb(k1,k4)*za(k4,k5)
      zab3(k1,k2,k3,k4,k5)=
     & za(k1,k2)*zb(k2,k5)+za(k1,k3)*zb(k3,k5)+za(k1,k4)*zb(k4,k5)
      zbaba222(k1,k2,k3,k4,k5,k6,k7,k8)=
     & +zba2(k1,k2,k3,k4)*zba2(k4,k6,k7,k8)
     & +zba2(k1,k2,k3,k5)*zba2(k5,k6,k7,k8)
      zbab32(k1,k2,k3,k4,k5,k6)=
     & +zba3(k1,k2,k3,k4,k5)*zb(k5,k1)+zba3(k1,k2,k3,k4,k6)*zb(k6,k1)
      zaba32(k1,k2,k3,k4,k5,k6)=
     & +zab3(k1,k2,k3,k4,k5)*za(k5,k1)+zab3(k1,k2,k3,k4,k6)*za(k6,k1)
      zbab22(k1,k2,k3,k4,k5,k6)=
     & +zba2(k1,k2,k3,k4)*zb(k4,k6)+zba2(k1,k2,k3,k5)*zb(k5,k6)
      zaba22(k1,k2,k3,k4,k5,k6)=
     & +zab2(k1,k2,k3,k4)*za(k4,k6)+zab2(k1,k2,k3,k5)*za(k5,k6)
c---statement functions
      s34=s(k3,k4)
      s56=s(k5,k6)
      s78=s(k7,k8)
      t178=s(k1,k7)+s(k1,k8)+s(k7,k8)
      t278=s(k2,k7)+s(k2,k8)+s(k7,k8)
      t3456=s(k3,k4)+s(k3,k5)+s(k3,k6)+s(k4,k5)+s(k4,k6)+s(k5,k6)
      do i=1,8
      do j=i+1,8
         iza(i,j)=cone/za(i,j)
         izb(i,j)=cone/zb(i,j)
         iza(j,i)=-iza(i,j)
         izb(j,i)=-izb(i,j)
      enddo
      enddo
      a8b(2,2)= + iza(k1,k7)*t278**(-1)*t3456**(-1)*s34**(-1)*s56**(-1)
     &  * (  - zba2(k4,k6,k5,k3)*zba3(k5,k2,k7,k8,k1)*za(k1,k6)*zb(k2,
     &    k8)*zb(k7,k8)*s78**(-1) + zba2(k5,k4,k3,k6)*zba3(k4,k2,k7,k8,
     &    k1)*za(k1,k3)*zb(k2,k8)*zb(k7,k8)*s78**(-1) + za(k3,k6)*zb(k2
     &    ,k8)*zb(k4,k5)*zb(k7,k8)*zaba32(k1,k2,k7,k8,k5,k6)*s78**(-1)
     &     )
      a8b(2,2) = a8b(2,2) + iza(k1,k7)*iza(k1,k8)*iza(k2,k8)*t278**(-1)
     & *t3456**(-1)*s34**(-1)*s56**(-1) * (  - zba2(k4,k6,k5,k3)*zba2(
     &    k7,k2,k8,k1)*zba3(k5,k2,k7,k8,k1)*za(k1,k6) + zba2(k5,k4,k3,
     &    k6)*zba2(k7,k2,k8,k1)*zba3(k4,k2,k7,k8,k1)*za(k1,k3) + zba2(
     &    k7,k2,k8,k1)*za(k3,k6)*zb(k4,k5)*zaba32(k1,k2,k7,k8,k5,k6) )
      a8b(2,2) = a8b(2,2) + iza(k1,k8)*t278**(-1)*t3456**(-1)*s34**(-1)
     & *s56**(-1) * (  - zba2(k4,k6,k5,k3)*zba3(k5,k2,k7,k8,k1)*za(k1,
     &    k6)*zb(k2,k7)*zb(k7,k8)*s78**(-1) + zba2(k5,k4,k3,k6)*zba3(k4
     &    ,k2,k7,k8,k1)*za(k1,k3)*zb(k2,k7)*zb(k7,k8)*s78**(-1) + za(k3
     &    ,k6)*zb(k2,k7)*zb(k4,k5)*zb(k7,k8)*zaba32(k1,k2,k7,k8,k5,k6)*
     &    s78**(-1) )

      a8b(1,1)= + izb(k1,k7)*izb(k2,k7)*izb(k2,k8)*t178**(-1)*
     & t3456**(-1)*s34**(-1)*s56**(-1) * (  - zba2(k2,k1,k7,k8)*zba2(k4
     &    ,k6,k5,k3)*zba3(k2,k1,k7,k8,k6)*zb(k2,k5) + zba2(k2,k1,k7,k8)
     &    *zba2(k5,k4,k3,k6)*zba3(k2,k1,k7,k8,k3)*zb(k2,k4) + zba2(k2,
     &    k1,k7,k8)*za(k3,k6)*zb(k4,k5)*zbab32(k2,k1,k7,k8,k5,k6) )
      a8b(1,1) = a8b(1,1) + izb(k2,k7)*t178**(-1)*t3456**(-1)*s34**(-1)
     & *s56**(-1) * ( zba2(k4,k6,k5,k3)*zba3(k2,k1,k7,k8,k6)*za(k1,k8)*
     &    za(k7,k8)*zb(k2,k5)*s78**(-1) - zba2(k5,k4,k3,k6)*zba3(k2,k1,
     &    k7,k8,k3)*za(k1,k8)*za(k7,k8)*zb(k2,k4)*s78**(-1) - za(k1,k8)
     &    *za(k3,k6)*za(k7,k8)*zb(k4,k5)*zbab32(k2,k1,k7,k8,k5,k6)*
     &    s78**(-1) )
      a8b(1,1) = a8b(1,1) + izb(k2,k8)*t178**(-1)*t3456**(-1)*s34**(-1)
     & *s56**(-1) * ( zba2(k4,k6,k5,k3)*zba3(k2,k1,k7,k8,k6)*za(k1,k7)*
     &    za(k7,k8)*zb(k2,k5)*s78**(-1) - zba2(k5,k4,k3,k6)*zba3(k2,k1,
     &    k7,k8,k3)*za(k1,k7)*za(k7,k8)*zb(k2,k4)*s78**(-1) - za(k1,k7)
     &    *za(k3,k6)*za(k7,k8)*zb(k4,k5)*zbab32(k2,k1,k7,k8,k5,k6)*
     &    s78**(-1) )

      a8b(2,1)= + iza(k1,k7)*iza(k7,k8)*izb(k2,k8)*izb(k7,k8)*
     & t3456**(-1)*s34**(-1)*s56**(-1) * ( zba2(k2,k6,k5,k1)*za(k1,k8)*
     &    za(k3,k6)*zb(k2,k7)*zb(k4,k5) - zba2(k4,k6,k5,k3)*za(k1,k6)*
     &    za(k1,k8)*zb(k2,k5)*zb(k2,k7) + zba2(k5,k4,k3,k6)*za(k1,k3)*
     &    za(k1,k8)*zb(k2,k4)*zb(k2,k7) )
      a8b(2,1) = a8b(2,1) + iza(k1,k7)*iza(k7,k8)*izb(k7,k8)*t178**(-1)
     & *t3456**(-1)*s34**(-1)*s56**(-1) * ( zba2(k4,k6,k5,k3)*zba2(k7,
     &    k1,k8,k6)*za(k1,k8)**2*zb(k2,k5) - zba2(k5,k4,k3,k6)*zba2(k7,
     &    k1,k8,k3)*za(k1,k8)**2*zb(k2,k4) - za(k1,k8)**2*za(k3,k6)*zb(
     &    k4,k5)*zbab22(k7,k1,k8,k5,k6,k2) )
      a8b(2,1) = a8b(2,1) + iza(k7,k8)*izb(k2,k8)*izb(k7,k8)*t278**(-1)
     & *t3456**(-1)*s34**(-1)*s56**(-1) * (  - zba2(k4,k2,k7,k8)*zba2(
     &    k5,k4,k3,k6)*za(k1,k3)*zb(k2,k7)**2 + zba2(k4,k6,k5,k3)*zba2(
     &    k5,k2,k7,k8)*za(k1,k6)*zb(k2,k7)**2 + za(k3,k6)*zb(k2,k7)**2*
     &    zb(k4,k5)*zaba22(k1,k5,k6,k2,k7,k8) )

      a8b(1,2)= + iza(k2,k8)*iza(k7,k8)*izb(k1,k7)*izb(k7,k8)*
     & t3456**(-1)*s34**(-1)*s56**(-1) * ( zba2(k4,k2,k8,k7)*zba2(k5,k4
     &    ,k3,k6)*zba2(k8,k1,k7,k3) - zba2(k4,k6,k5,k3)*zba2(k5,k2,k8,
     &    k7)*zba2(k8,k1,k7,k6) - za(k3,k6)*zb(k4,k5)*zbaba222(k8,k1,k7
     &    ,k5,k6,k2,k8,k7) )
      a8b(1,2) = a8b(1,2) + iza(k2,k8)*iza(k7,k8)*izb(k7,k8)*t278**(-1)
     & *t3456**(-1)*s34**(-1)*s56**(-1) * (  - zba2(k4,k2,k8,k7)*zba2(
     &    k5,k4,k3,k6)*za(k1,k3)*za(k2,k7)*zb(k2,k8) + zba2(k4,k6,k5,k3
     &    )*zba2(k5,k2,k8,k7)*za(k1,k6)*za(k2,k7)*zb(k2,k8) + za(k2,k7)
     &    *za(k3,k6)*zb(k2,k8)*zb(k4,k5)*zaba22(k1,k5,k6,k2,k8,k7) )
      a8b(1,2) = a8b(1,2) + iza(k7,k8)*izb(k1,k7)*izb(k7,k8)*t178**(-1)
     & *t3456**(-1)*s34**(-1)*s56**(-1) * ( zba2(k4,k6,k5,k3)*zba2(k8,
     &    k1,k7,k6)*za(k1,k7)*zb(k1,k8)*zb(k2,k5) - zba2(k5,k4,k3,k6)*
     &    zba2(k8,k1,k7,k3)*za(k1,k7)*zb(k1,k8)*zb(k2,k4) - za(k1,k7)*
     &    za(k3,k6)*zb(k1,k8)*zb(k4,k5)*zbab22(k8,k1,k7,k5,k6,k2) )

c     To fix with our standard notation removing the overall factor of i
      a8b(:,:)=-a8b(:,:)
      return
      end
