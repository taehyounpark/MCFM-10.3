!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a8treea(k1,k2,k3,k4,k5,k6,k7,k8,za,zb,a8a)
      implicit none
c     (a)-type amplitude calculated in DKS notation
c     0->u(q1)+ubar(q2)+l(q3)+a(q4)+a'(q5)+l'(q6)+g(p7)+g(p8)
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      real(dp):: s34,s56,s78,t1567,t178,t278,t156,t234
      complex(dp):: a8a(2,2)
      complex(dp):: iza(8,8),izb(8,8),zba2,zba3
      integer i,j,k1,k2,k3,k4,k5,k6,k7,k8
c---statement functions
      zba2(k1,k2,k3,k4)=zb(k1,k2)*za(k2,k4)+zb(k1,k3)*za(k3,k4)
      zba3(k1,k2,k3,k4,k5)=
     & zb(k1,k2)*za(k2,k5)+zb(k1,k3)*za(k3,k5)+zb(k1,k4)*za(k4,k5)
c---statement functions
      s34=s(k3,k4)
      s56=s(k5,k6)
      s78=s(k7,k8)
      t156=s(k1,k5)+s(k1,k6)+s(k5,k6)
      t234=s(k2,k3)+s(k2,k4)+s(k3,k4)
      t178=s(k1,k7)+s(k7,k8)+s(k1,k8)
      t278=s(k2,k7)+s(k7,k8)+s(k2,k8)
      t1567=s(k1,k5)+s(k1,k6)+s(k1,k7)+s(k5,k6)+s(k5,k7)+s(k6,k7)
      do i=1,8
      do j=i+1,8
         iza(i,j)=cone/za(i,j)
         izb(i,j)=cone/zb(i,j)
         iza(j,i)=-iza(i,j)
         izb(j,i)=-izb(i,j)
      enddo
      enddo
      a8a(2,2)= + iza(k1,k7)*t278**(-1) * ( zba2(k5,k1,k6,k3)*zba3(k4,
     &    k2,k7,k8,k1)*za(k1,k6)*zb(k2,k8)*zb(k7,k8)*s78**(-1)*
     &    s34**(-1)*s56**(-1)*t156**(-1) )
      a8a(2,2) = a8a(2,2) + iza(k1,k7)*iza(k1,k8)*iza(k2,k8)*t278**(-1)
     &  * ( zba2(k5,k1,k6,k3)*zba2(k7,k2,k8,k1)*zba3(k4,k2,k7,k8,k1)*
     &    za(k1,k6)*s34**(-1)*s56**(-1)*t156**(-1) )
      a8a(2,2) = a8a(2,2) + iza(k1,k7)*iza(k1,k8)*iza(k2,k8) * (  -
     &    zba2(k4,k2,k8,k1)*zba3(k7,k1,k5,k6,k3)*za(k1,k6)*za(k6,k1)*
     &    zb(k5,k6)*s34**(-1)*s56**(-1)*t156**(-1)*t1567**(-1) )
      a8a(2,2) = a8a(2,2) + iza(k1,k7)*iza(k1,k8) * ( zba2(k7,k5,k6,k1)
     &    *zba2(k8,k2,k4,k3)*za(k1,k6)*za(k6,k1)*zb(k2,k4)*zb(k5,k6)*
     &    s34**(-1)*s56**(-1)*t156**(-1)*t234**(-1)*t1567**(-1) )
      a8a(2,2) = a8a(2,2) + iza(k1,k7) * ( zba2(k8,k2,k4,k3)*za(k1,k6)*
     &    za(k6,k1)*zb(k2,k4)*zb(k5,k6)*zb(k7,k8)*s78**(-1)*s34**(-1)*
     &    s56**(-1)*t156**(-1)*t234**(-1) )
      a8a(2,2) = a8a(2,2) + iza(k1,k8)*t278**(-1) * ( zba2(k5,k1,k6,k3)
     &    *zba3(k4,k2,k7,k8,k1)*za(k1,k6)*zb(k2,k7)*zb(k7,k8)*s78**(-1)
     &    *s34**(-1)*s56**(-1)*t156**(-1) )
      a8a(2,2) = a8a(2,2) + iza(k1,k8) * ( zba2(k7,k2,k4,k3)*za(k1,k6)*
     &    za(k6,k1)*zb(k2,k4)*zb(k5,k6)*zb(k7,k8)*s78**(-1)*s34**(-1)*
     &    s56**(-1)*t156**(-1)*t234**(-1) )

      a8a(1,1)= + izb(k1,k7)*izb(k2,k7)*izb(k2,k8)*t178**(-1)*s34**(-1)
     & *s56**(-1) * (  - zba2(k2,k1,k7,k8)*zba2(k5,k2,k4,k3)*zba3(k2,k1
     &    ,k7,k8,k6)*zb(k2,k4)*t234**(-1) )
      a8a(1,1) = a8a(1,1) + izb(k1,k7)*izb(k2,k7)*izb(k2,k8)*s34**(-1)*
     & s56**(-1) * (  - zba2(k2,k1,k7,k6)*zba3(k5,k1,k6,k7,k8)*za(k4,k3
     &    )*zb(k2,k4)**2*t234**(-1)*t1567**(-1) )
      a8a(1,1) = a8a(1,1) + izb(k2,k7)*t178**(-1)*s34**(-1)*s56**(-1)
     &  * ( zba2(k5,k2,k4,k3)*zba3(k2,k1,k7,k8,k6)*za(k1,k8)*za(k7,k8)*
     &    zb(k2,k4)*s78**(-1)*t234**(-1) )
      a8a(1,1) = a8a(1,1) + izb(k2,k7)*s34**(-1)*s56**(-1) * ( zba2(k5,
     &    k1,k6,k8)*za(k1,k6)*za(k4,k3)*za(k7,k8)*zb(k2,k4)**2*
     &    s78**(-1)*t156**(-1)*t234**(-1) )
      a8a(1,1) = a8a(1,1) + izb(k2,k7)*izb(k2,k8)*s34**(-1)*s56**(-1)
     &  * (  - zba2(k2,k3,k4,k8)*zba2(k5,k1,k6,k7)*za(k1,k6)*za(k4,k3)*
     &    zb(k2,k4)**2*t156**(-1)*t234**(-1)*t1567**(-1) )
      a8a(1,1) = a8a(1,1) + izb(k2,k8)*t178**(-1)*s34**(-1)*s56**(-1)
     &  * ( zba2(k5,k2,k4,k3)*zba3(k2,k1,k7,k8,k6)*za(k1,k7)*za(k7,k8)*
     &    zb(k2,k4)*s78**(-1)*t234**(-1) )
      a8a(1,1) = a8a(1,1) + izb(k2,k8)*s34**(-1)*s56**(-1) * ( zba2(k5,
     &    k1,k6,k7)*za(k1,k6)*za(k4,k3)*za(k7,k8)*zb(k2,k4)**2*
     &    s78**(-1)*t156**(-1)*t234**(-1) )

      a8a(2,1)= + iza(k1,k7)*iza(k7,k8)*izb(k2,k8)*izb(k7,k8)*s34**(-1)
     & *s56**(-1) * ( zba3(k5,k1,k6,k7,k3)*za(k1,k6)*za(k1,k8)*zb(k2,k4
     &    )*zb(k2,k7)*t1567**(-1) )
      a8a(2,1) = a8a(2,1) + iza(k1,k7)*iza(k7,k8)*izb(k7,k8)*t178**(-1)
     & *s34**(-1)*s56**(-1) * ( zba2(k5,k2,k4,k3)*zba2(k7,k1,k8,k6)*za(
     &    k1,k8)**2*zb(k2,k4)*t234**(-1) )
      a8a(2,1) = a8a(2,1) + iza(k1,k7)*iza(k7,k8)*izb(k7,k8)*s34**(-1)*
     & s56**(-1) * ( zba2(k7,k2,k4,k3)*zba3(k5,k1,k6,k7,k8)*za(k1,k6)*
     &    za(k1,k8)*zb(k2,k4)*t234**(-1)*t1567**(-1) )
      a8a(2,1) = a8a(2,1) + iza(k7,k8)*izb(k2,k8)*izb(k7,k8)*t278**(-1)
     & *s34**(-1)*s56**(-1) * (  - zba2(k4,k2,k7,k8)*zba2(k5,k1,k6,k3)*
     &    za(k1,k6)*zb(k2,k7)**2*t156**(-1) )
      a8a(2,1) = a8a(2,1) + iza(k7,k8)*izb(k2,k8)*izb(k7,k8)*s34**(-1)*
     & s56**(-1) * ( zba2(k5,k1,k6,k8)*zba3(k7,k1,k5,k6,k3)*za(k1,k6)*
     &    zb(k2,k4)*zb(k2,k7)*t156**(-1)*t1567**(-1) )
      a8a(2,1) = a8a(2,1) + iza(k7,k8)*izb(k7,k8)*s34**(-1)*s56**(-1)
     &  * ( zba2(k5,k1,k6,k8)*zba2(k7,k2,k4,k3)*zba3(k7,k1,k5,k6,k8)*
     &    za(k1,k6)*zb(k2,k4)*t156**(-1)*t234**(-1)*t1567**(-1) )

      a8a(1,2)= + iza(k2,k8)*iza(k7,k8)*izb(k1,k7)*izb(k7,k8)*s34**(-1)
     & *s56**(-1) * ( zba2(k4,k2,k8,k7)*zba2(k8,k1,k7,k6)*zba3(k5,k1,k6
     &    ,k7,k3)*t1567**(-1) )
      a8a(1,2) = a8a(1,2) + iza(k2,k8)*iza(k7,k8)*izb(k7,k8)*t278**(-1)
     & *s34**(-1)*s56**(-1) * (  - zba2(k4,k2,k8,k7)*zba2(k5,k1,k6,k3)*
     &    za(k1,k6)*za(k2,k7)*zb(k2,k8)*t156**(-1) )
      a8a(1,2) = a8a(1,2) + iza(k2,k8)*iza(k7,k8)*izb(k7,k8)*s34**(-1)*
     & s56**(-1) * ( zba2(k4,k2,k8,k7)*zba2(k5,k1,k6,k7)*zba2(k8,k2,k4,
     &    k3)*za(k1,k6)*t156**(-1)*t1567**(-1) )
      a8a(1,2) = a8a(1,2) + iza(k7,k8)*izb(k1,k7)*izb(k7,k8)*t178**(-1)
     & *s34**(-1)*s56**(-1) * ( zba2(k5,k2,k4,k3)*zba2(k8,k1,k7,k6)*za(
     &    k1,k7)*zb(k1,k8)*zb(k2,k4)*t234**(-1) )
      a8a(1,2) = a8a(1,2) + iza(k7,k8)*izb(k1,k7)*izb(k7,k8)*s34**(-1)*
     & s56**(-1) * (  - zba2(k5,k1,k6,k7)*zba2(k8,k1,k7,k6)*zba2(k8,k2,
     &    k4,k3)*zb(k2,k4)*t234**(-1)*t1567**(-1) )
      a8a(1,2) = a8a(1,2) + iza(k7,k8)*izb(k7,k8)*s34**(-1)*s56**(-1)
     &  * ( zba2(k5,k1,k6,k7)*zba2(k8,k2,k4,k3)*zba3(k8,k1,k5,k6,k7)*
     &    za(k1,k6)*zb(k2,k4)*t156**(-1)*t234**(-1)*t1567**(-1) )

c     To fix with our standard notation removing the overall factor of i
      a8a(:,:)=-a8a(:,:)
      return
      end
