!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine ZZmbc(st,j1,j2,j3,j4,j5,j6,za,zb,bcoeff)
      implicit none
      include 'types.f'

c----Bubble coefficients extracted from BDK 11.5, 11.8
c----after performing a transformation
      integer:: j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'blabels.f'
      logical:: ggZZuse6d
      character(len=9):: st
      complex(dp):: bcoeff(8),zab2,funcb34pp,funcb134pp
      complex(dp):: funcb134mp,funcb34mp,
     & oldfuncb134mp,oldfuncb34mp
      real(dp):: IDelta
      common/ggZZuse6d/ggZZuse6d
!$omp threadprivate(/ggZZuse6d/)

c----statement functions
      zab2(j1,j2,j3,j4)=+za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      funcb34pp(j1,j2,j3,j4,j5,j6)=za(j3,j4)/(za(j1,j2)**2*za(j5,j6))
     & *((za(j1,j5)*zb(j1,j4)/(s(j1,j3)+s(j1,j4)))**2
     &  +(za(j2,j5)*zb(j2,j4)/(s(j2,j3)+s(j2,j4)))**2)
     & +two*za(j1,j5)*za(j2,j5)/(za(j1,j2)**3*za(j5,j6))
     & *(za(j3,j1)*zb(j4,j1)/(s(j1,j3)+s(j1,j4))
     &  +za(j2,j3)*zb(j4,j2)/(s(j2,j3)+s(j2,j4)))

      funcb134pp(j1,j2,j3,j4,j5,j6)=(two*za(j1,j3)*za(j2,j5)
     & *(za(j2,j3)*zb(j2,j6)*za(j5,j6)/(s(j2,j5)+s(j2,j6))
     &  -za(j1,j5)*zb(j1,j4)*za(j3,j4)/(s(j1,j3)+s(j1,j4)))
     & -za(j1,j2)
     & *((za(j3,j4)*za(j1,j5)*zb(j1,j4)/(s(j1,j3)+s(j1,j4)))**2
     &  +(za(j2,j3)*za(j5,j6)*zb(j2,j6)/(s(j2,j5)+s(j2,j6)))**2)
     & -(za(j1,j3)*za(j2,j5)+za(j2,j3)*za(j1,j5))*za(j3,j5))
     & /(za(j1,j2)**3*za(j3,j4)*za(j5,j6))
c---end statement functions


      IDelta=s(j1,j2)**2+s(j3,j4)**2+s(j5,j6)**2
     & -two*s(j1,j2)*s(j3,j4)
     & -two*s(j3,j4)*s(j5,j6)
     & -two*s(j1,j2)*s(j5,j6)
      IDelta=1._dp/IDelta

      bcoeff(:)=czip


      if     (st=='q+qb-g+g+') then

      bcoeff(b34)=funcb34pp(j1,j2,j3,j4,j5,j6)
      bcoeff(b56)=funcb34pp(j1,j2,j5,j6,j3,j4)
      bcoeff(b134)=funcb134pp(j1,j2,j3,j4,j5,j6)
c      bcoeff(b234)=funcb134pp(j2,j1,j3,j4,j5,j6)
      bcoeff(b234)=-bcoeff(b34)-bcoeff(b56)-bcoeff(b134)

      bcoeff(rat)=(
     &  (za(j1,j5)**2*zb(j1,j4)**2/(s(j1,j3)+s(j1,j4))
     &  +za(j2,j5)**2*zb(j2,j4)**2/(s(j2,j3)+s(j2,j4)))
     & /(za(j5,j6)*zb(j3,j4))
     & +(za(j3,j1)**2*zb(j1,j6)**2/(s(j1,j5)+s(j1,j6))
     &  +za(j3,j2)**2*zb(j2,j6)**2/(s(j2,j5)+s(j2,j6)))
     & /(zb(j5,j6)*za(j3,j4))
     & -za(j3,j5)**2/(za(j3,j4)*za(j5,j6))
     & +zb(j4,j6)**2/(zb(j3,j4)*zb(j5,j6)))/za(j1,j2)**2

      elseif (st=='q+qb-g-g+') then


c      bcoeff(b12zm)= + IDelta**2 * ( three*zba2(j4,j1,j2,j3)*zba2(j2,j3,
c     &    j4,j1)*zba2(j6,j1,j2,j5)/zba2(j1,j3,j4,j2)*s(j3,j4) + three*
c     &    zba2(j4,j1,j2,j3)*zba2(j2,j4,j3,j1)*zba2(j6,j1,j2,j5)/zba2(
c     &    j1,j3,j4,j2)*s(j1,j2) - three*zba2(j4,j1,j2,j3)*zba2(j2,j3,j4,
c     &    j1)*zba2(j6,j1,j2,j5)/zba2(j1,j3,j4,j2)*s(j5,j6) - three*
c     &    zba2(j4,j1,j2,j3)*zba2(j2,j5,j6,j1)*zba2(j6,j1,j2,j5)/zba2(
c     &    j1,j5,j6,j2)*s(j3,j4) + three*zba2(j4,j1,j2,j3)*zba2(j2,j5,j6,
c     &    j1)*zba2(j6,j1,j2,j5)/zba2(j1,j5,j6,j2)*s(j1,j2) + three*
c     &    zba2(j4,j1,j2,j3)*zba2(j2,j5,j6,j1)*zba2(j6,j1,j2,j5)/zba2(
c     &    j1,j5,j6,j2)*s(j5,j6) + three*zba2(j4,j1,j2,j3)*zba2(j2,j4,j3,
c     &    j1)*zba2(j6,j1,j2,j5)/zab2(j2,j3,j4,j1)*s(j3,j4) + three*
c     &    zba2(j4,j1,j2,j3)*zba2(j2,j4,j3,j1)*zba2(j6,j1,j2,j5)/zba2(
c     &    j1,j4,j3,j2)*s(j1,j2) - three*zba2(j4,j1,j2,j3)*zba2(j2,j4,j3,
c     &    j1)*zba2(j6,j1,j2,j5)/zab2(j2,j3,j4,j1)*s(j5,j6) - three*
c     &    zba2(j4,j1,j2,j3)*zba2(j2,j6,j5,j1)*zba2(j6,j1,j2,j5)/zba2(
c     &    j1,j6,j5,j2)*s(j3,j4) )
c     & + IDelta**2 * ( three*zba2(j4,j2,j1,
c     &    j3)*zba2(j2,j6,j5,j1)*zba2(j6,j1,j2,j5)/zba2(j1,j6,j5,j2)*s(
c     &    j1,j2) + three*zba2(j4,j1,j2,j3)*zba2(j2,j6,j5,j1)*zba2(j6,j2,
c     &    j1,j5)/zba2(j1,j6,j5,j2)*s(j5,j6) )
c     & + IDelta * (  - 2._dp*zb(j3,j4)*zb(
c     &    j4,j1)*zb(j6,j5)*za(j3,j4)*za(j3,j5)*za(j2,j5)/zba2(j1,j3,j4
c     &    ,j2)**3*t(j3,j4,j1) + 2._dp*zb(j3,j4)*zb(j4,j1)*zb(j6,j5)*za(
c     &    j3,j4)*za(j3,j5)*za(j2,j5)/zba2(j1,j3,j4,j2)**3*t(j3,j4,j2)
c     &     - 2._dp*zb(j3,j4)*zb(j4,j6)*zb(j1,j6)*za(j3,j4)*za(j3,j2)*za(
c     &    j6,j5)/zab2(j2,j3,j4,j1)**3*t(j4,j3,j1) + 2._dp*zb(j3,j4)*zb(
c     &    j4,j6)*zb(j1,j6)*za(j3,j4)*za(j3,j2)*za(j6,j5)/zba2(j1,j4,j3
c     &    ,j2)**3*t(j4,j3,j2) + zb(j3,j4)*zb(j4,j6)*za(j3,j4)*za(j3,j5)
c     &    /zba2(j1,j3,j4,j2)**2*s(j3,j4) - zb(j3,j4)*zb(j4,j6)*za(j3,
c     &    j4)*za(j3,j5)/zba2(j1,j3,j4,j2)**2*s(j1,j2) - zb(j3,j4)*zb(
c     &    j4,j6)*za(j3,j4)*za(j3,j5)/zba2(j1,j3,j4,j2)**2*s(j5,j6) +
c     &    zb(j3,j4)*zb(j4,j6)*za(j3,j4)*za(j3,j5)/zab2(j2,j3,j4,j1)**2
c     &    *s(j3,j4) - zb(j3,j4)*zb(j4,j6)*za(j3,j4)*za(j3,j5)/zba2(j1,
c     &    j4,j3,j2)**2*s(j1,j2) - zb(j3,j4)*zb(j4,j6)*za(j3,j4)*za(j3,
c     &    j5)/zab2(j2,j3,j4,j1)**2*s(j5,j6) - 2._dp*zb(j3,j4)*zb(j1,j6)
c     &    *zb(j6,j5)*za(j3,j2)*za(j3,j5)*za(j6,j5)/zba2(j1,j5,j6,j2)**
c     &    3*t(j5,j6,j1) )
c     & + IDelta * ( 2._dp*zb(j3,j4)*zb(j1,
c     &    j6)*zb(j6,j5)*za(j3,j2)*za(j3,j5)*za(j6,j5)/zba2(j1,j5,j6,j2
c     &    )**3*t(j5,j6,j2) - 2._dp*zb(j3,j4)*zb(j1,j6)*za(j3,j2)*za(j3,
c     &    j5)/zab2(j2,j3,j4,j1)**3*t(j4,j3,j1)*t(j4,j3,j2) + 2._dp*zb(
c     &    j3,j4)*zb(j1,j6)*za(j3,j2)*za(j3,j5)/zab2(j2,j3,j4,j1)**3*t(
c     &    j4,j3,j2)**2 + 2._dp*zb(j3,j4)*zb(j2,j6)*za(j3,j2)*za(j3,j5)*
c     &    izab2(j2,j3,j4,j1)**2*t(j4,j3,j1) - 2._dp*zb(j3,j4)*zb(j2,j6)*
c     &    za(j3,j2)*za(j3,j5)/zab2(j2,j3,j4,j1)**2*t(j4,j3,j2) - one/
c     &    2._dp*zb(j3,j4)*za(j3,j5)**2/za(j6,j5)*zba2(j2,j4,j3,j1)*
c     &    izba2(j1,j3,j4,j2) - one/2._dp*zb(j3,j4)*za(j3,j5)**2/za(j6,
c     &    j5)*zba2(j2,j4,j3,j1)/zab2(j2,j3,j4,j1) + zb(j3,j4)*za(j3,j5
c     &    )**2/za(j6,j5)/zab2(j2,j3,j4,j1)**2*s(j3,j4)*t(j4,j3,j2) -
c     &    zb(j3,j4)*za(j3,j5)**2/za(j6,j5)/zab2(j2,j3,j4,j1)**2*s(j1,
c     &    j2)*t(j4,j3,j2) - zb(j3,j4)*za(j3,j5)**2/za(j6,j5)/zba2(j1,
c     &    j4,j3,j2)**2*s(j5,j6)*t(j4,j3,j2) - 2._dp*zb(j4,j1)*zb(j4,j6)*
c     &    zb(j6,j5)*za(j3,j4)*za(j2,j5)*za(j6,j5)/zba2(j1,j6,j5,j2)**3
c     &    *t(j6,j5,j1) )
c     & + IDelta * ( 2._dp*zb(j4,j1)*zb(j4,
c     &    j6)*zb(j6,j5)*za(j3,j4)*za(j2,j5)*za(j6,j5)/zba2(j1,j6,j5,j2
c     &    )**3*t(j6,j5,j2) + 2._dp*zb(j4,j1)*zb(j4,j6)*za(j3,j4)*za(j1,
c     &    j5)/zba2(j1,j3,j4,j2)**2*t(j3,j4,j1) - 2._dp*zb(j4,j1)*zb(j4,
c     &    j6)*za(j3,j4)*za(j1,j5)/zba2(j1,j3,j4,j2)**2*t(j3,j4,j2) - 2.
c     &    D0*zb(j4,j1)*zb(j4,j6)*za(j3,j4)*za(j2,j5)/zba2(j1,j3,j4,j2)
c     &    **3*t(j3,j4,j1)**2 + 2._dp*zb(j4,j1)*zb(j4,j6)*za(j3,j4)*za(j2
c     &    ,j5)/zba2(j1,j3,j4,j2)**3*t(j3,j4,j1)*t(j3,j4,j2) - 2._dp*zb(
c     &    j4,j1)*zb(j6,j5)*za(j3,j5)*za(j2,j5)/zba2(j1,j6,j5,j2)**3*t(
c     &    j6,j5,j1)*t(j6,j5,j2) + 2._dp*zb(j4,j1)*zb(j6,j5)*za(j3,j5)*
c     &    za(j2,j5)/zba2(j1,j6,j5,j2)**3*t(j6,j5,j2)**2 + 2._dp*zb(j4,
c     &    j2)*zb(j6,j5)*za(j3,j5)*za(j2,j5)/zba2(j1,j6,j5,j2)**2*t(j6,
c     &    j5,j1) - 2._dp*zb(j4,j2)*zb(j6,j5)*za(j3,j5)*za(j2,j5)/zba2(
c     &    j1,j6,j5,j2)**2*t(j6,j5,j2) - one/2._dp*zb(j4,j6)**2*za(j3,j4
c     &    )/zb(j6,j5)*zba2(j2,j4,j3,j1)/zba2(j1,j3,j4,j2) - one/2._dp
c     &    *zb(j4,j6)**2*za(j3,j4)/zb(j6,j5)*zba2(j2,j4,j3,j1)/zba2(j1
c     &    ,j4,j3,j2) )
c     & + IDelta * ( zb(j4,j6)**2*za(j3,j4)
c     &    /zb(j6,j5)/zba2(j1,j3,j4,j2)**2*s(j3,j4)*t(j3,j4,j1) - zb(
c     &    j4,j6)**2*za(j3,j4)/zb(j6,j5)/zba2(j1,j3,j4,j2)**2*s(j1,j2)
c     &    *t(j3,j4,j1) - zb(j4,j6)**2*za(j3,j4)/zb(j6,j5)/zba2(j1,j3,
c     &    j4,j2)**2*s(j5,j6)*t(j3,j4,j1) - one/2._dp*zb(j4,j6)**2*za(j6
c     &    ,j5)/zb(j3,j4)*zba2(j2,j6,j5,j1)/zba2(j1,j6,j5,j2) - one/2.
c     &    D0*zb(j4,j6)**2*za(j6,j5)/zb(j3,j4)*zba2(j2,j5,j6,j1)/zba2(
c     &    j1,j5,j6,j2) - zb(j4,j6)**2*za(j6,j5)/zb(j3,j4)/zba2(j1,j5,
c     &    j6,j2)**2*s(j3,j4)*t(j5,j6,j1) - zb(j4,j6)**2*za(j6,j5)/zb(
c     &    j3,j4)/zba2(j1,j5,j6,j2)**2*s(j1,j2)*t(j5,j6,j1) + zb(j4,j6)
c     &    **2*za(j6,j5)/zb(j3,j4)/zba2(j1,j5,j6,j2)**2*s(j5,j6)*t(j5,
c     &    j6,j1) + 2._dp*zb(j4,j6)*zb(j1,j6)*za(j3,j1)*za(j6,j5)/zba2(
c     &    j1,j5,j6,j2)**2*t(j5,j6,j1) - 2._dp*zb(j4,j6)*zb(j1,j6)*za(j3,
c     &    j1)*za(j6,j5)/zba2(j1,j5,j6,j2)**2*t(j5,j6,j2) - 2._dp*zb(j4,
c     &    j6)*zb(j1,j6)*za(j3,j2)*za(j6,j5)/zba2(j1,j5,j6,j2)**3*t(j5,
c     &    j6,j1)**2 )
c     & + IDelta * ( 2._dp*zb(j4,j6)*zb(j1,
c     &    j6)*za(j3,j2)*za(j6,j5)/zba2(j1,j5,j6,j2)**3*t(j5,j6,j1)*t(
c     &    j5,j6,j2) - zb(j4,j6)*zb(j6,j5)*za(j3,j5)*za(j6,j5)/zba2(j1,
c     &    j6,j5,j2)**2*s(j3,j4) - zb(j4,j6)*zb(j6,j5)*za(j3,j5)*za(j6,
c     &    j5)/zba2(j1,j6,j5,j2)**2*s(j1,j2) + zb(j4,j6)*zb(j6,j5)*za(
c     &    j3,j5)*za(j6,j5)/zba2(j1,j6,j5,j2)**2*s(j5,j6) - zb(j4,j6)*
c     &    zb(j6,j5)*za(j3,j5)*za(j6,j5)/zba2(j1,j5,j6,j2)**2*s(j3,j4)
c     &     - zb(j4,j6)*zb(j6,j5)*za(j3,j5)*za(j6,j5)/zba2(j1,j5,j6,j2)
c     &    **2*s(j1,j2) + zb(j4,j6)*zb(j6,j5)*za(j3,j5)*za(j6,j5)/zba2(
c     &    j1,j5,j6,j2)**2*s(j5,j6) + zb(j4,j6)*za(j3,j5)*zba2(j2,j3,j4,
c     &    j1)/zba2(j1,j3,j4,j2) + zb(j4,j6)*za(j3,j5)*zba2(j2,j4,j3,j1
c     &    )/zab2(j2,j3,j4,j1) + zb(j4,j6)*za(j3,j5)*zba2(j2,j6,j5,j1)*
c     &    izba2(j1,j6,j5,j2) + zb(j4,j6)*za(j3,j5)*zba2(j2,j5,j6,j1)*
c     &    izba2(j1,j5,j6,j2) - one/2._dp*zb(j6,j5)*za(j3,j5)**2/za(j3,
c     &    j4)*zba2(j2,j6,j5,j1)/zba2(j1,j6,j5,j2) - one/2._dp*zb(j6,j5
c     &    )*za(j3,j5)**2/za(j3,j4)*zba2(j2,j5,j6,j1)/zba2(j1,j5,j6,j2
c     &    ) )
c     & + IDelta * (  - zb(j6,j5)*za(j3,j5)
c     &    **2/za(j3,j4)/zba2(j1,j6,j5,j2)**2*s(j3,j4)*t(j6,j5,j2) -
c     &    zb(j6,j5)*za(j3,j5)**2/za(j3,j4)/zba2(j1,j6,j5,j2)**2*s(j1,
c     &    j2)*t(j6,j5,j2) + zb(j6,j5)*za(j3,j5)**2/za(j3,j4)/zba2(j1,
c     &    j6,j5,j2)**2*s(j5,j6)*t(j6,j5,j2) )
c     & - 2._dp*zb(j4,j1)*zb(j4,j2)*za(j1,j5
c     & )*za(j2,j5)/zb(j3,j4)/za(j6,j5)/zba2(j1,j3,j4,j2)**2 + 2._dp*
c     &    zb(j4,j1)*za(j2,j5)/zb(j3,j4)/za(j6,j5)*zba2(j4,j3,j2,j5)*
c     &    izba2(j1,j3,j4,j2)**3*t(j3,j4,j1) - 2._dp*zb(j1,j6)*zb(j2,j6)*
c     &    za(j3,j1)*za(j3,j2)/zb(j6,j5)/za(j3,j4)/zab2(j2,j3,j4,j1)
c     &    **2 + 2._dp*zb(j1,j6)*za(j3,j2)/zb(j6,j5)/za(j3,j4)*zba2(j6,
c     &    j4,j1,j3)/zab2(j2,j3,j4,j1)**3*t(j4,j3,j2) - izb(j3,j4)/za(
c     &    j6,j5)*zba2(j4,j3,j2,j5)**2/zba2(j1,j3,j4,j2)**2 - izb(j6,j5
c     &    )/za(j3,j4)*zba2(j6,j4,j1,j3)**2/zab2(j2,j3,j4,j1)**2


c      write(6,*) 'old:b12zm',bcoeff(b12zm)

c      bcoeff(b12zm)= + IDelta**2
c     & * ( three*zba2(j4,j1,j2,j3)*zab2(j1,j3,j4,j2)*zba2(j6,j1,j2,j5)
c     & /zba2(j1,j3,j4,j2)*s(j3,j4)
c     &  + three*zba2(j4,j1,j2,j3)*zba2(j2,j4,j3,j1)*zba2(j6,j1,j2,j5)
c     & /zba2(j1,j3,j4,j2)*s(j1,j2)
c     &  -three*zba2(j4,j1,j2,j3)*zab2(j1,j3,j4,j2)*zba2(j6,j1,j2,j5)
c     & /zba2(j1,j3,j4,j2)*s(j5,j6)
c     & - three*zba2(j4,j1,j2,j3)*zba2(j2,j5,j6,j1)*zba2(j6,j1,j2,j5)
c     & /zba2(j1,j5,j6,j2)*s(j3,j4)
c     &  + three*zba2(j4,j1,j2,j3)*zba2(j2,j5,j6,j1)*zba2(j6,j1,j2,j5)
c     & /zba2(j1,j5,j6,j2)*s(j1,j2)
c     &  + three*zba2(j4,j1,j2,j3)*zba2(j2,j5,j6,j1)*zba2(j6,j1,j2,j5)
c     & /zba2(j1,j5,j6,j2)*s(j5,j6)
c     &  + three*zba2(j4,j1,j2,j3)*zba2(j2,j4,j3,j1)*zba2(j6,j1,j2,j5)
c     & /zab2(j2,j3,j4,j1)*s(j3,j4)
c     & + three*zba2(j4,j1,j2,j3)*zba2(j2,j4,j3,j1)*zba2(j6,j1,j2,j5)
c     & /zab2(j2,j3,j4,j1)*s(j1,j2)
c     & - three*zba2(j4,j1,j2,j3)*zba2(j2,j4,j3,j1)*zba2(j6,j1,j2,j5)
c     & /zab2(j2,j3,j4,j1)*s(j5,j6)
c     &  - three*zba2(j4,j1,j2,j3)*zba2(j2,j6,j5,j1)*zba2(j6,j1,j2,j5)
c     & /zba2(j1,j6,j5,j2)*s(j3,j4))

c     & + IDelta**2
c     & * ( three*zba2(j4,j2,j1,j3)*zba2(j2,j6,j5,j1)*zba2(j6,j1,j2,j5)
c     & /zba2(j1,j6,j5,j2)*s(j1,j2)
c     &  + three*zba2(j4,j1,j2,j3)*zba2(j2,j6,j5,j1)*zba2(j6,j2,j1,j5)
c     & /zba2(j1,j6,j5,j2)*s(j5,j6))

c     & + IDelta
c     &  * (-2._dp*zb(j3,j4)*zb(j4,j1)*zb(j6,j5)*za(j3,j4)*za(j3,j5)
c     & *za(j2,j5)/zba2(j1,j3,j4,j2)**3*t(j3,j4,j1)
c     &  + 2._dp*zb(j3,j4)*zb(j4,j1)*zb(j6,j5)*za(j3,j4)*za(j3,j5)
c     & *za(j2,j5)/zba2(j1,j3,j4,j2)**3*t(j3,j4,j2)
c     & -2._dp*zb(j3,j4)*zb(j4,j6)*zb(j1,j6)*za(j3,j4)*za(j3,j2)
c     & *za(j6,j5)/zab2(j2,j3,j4,j1)**3*t(j4,j3,j1)
c     &  +2._dp*zb(j3,j4)*zb(j4,j6)*zb(j1,j6)*za(j3,j4)*za(j3,j2)
c     &  *za(j6,j5)/zab2(j2,j3,j4,j1)**3*t(j4,j3,j2)
c     &  + zb(j3,j4)*zb(j4,j6)*za(j3,j4)*za(j3,j5)
c     & /zba2(j1,j3,j4,j2)**2*s(j3,j4)
c     & - zb(j3,j4)*zb(j4,j6)*za(j3,j4)*za(j3,j5)
c     & /zba2(j1,j3,j4,j2)**2*s(j1,j2)
c     &  - zb(j3,j4)*zb(j4,j6)*za(j3,j4)*za(j3,j5)
c     & /zba2(j1,j3,j4,j2)**2*s(j5,j6)
c     &  + zb(j3,j4)*zb(j4,j6)*za(j3,j4)*za(j3,j5)
c     & /zab2(j2,j3,j4,j1)**2*s(j3,j4)
c     &  -zb(j3,j4)*zb(j4,j6)*za(j3,j4)*za(j3,j5)
c     & /zab2(j2,j3,j4,j1)**2*s(j1,j2)
c     &  -zb(j3,j4)*zb(j4,j6)*za(j3,j4)*za(j3,j5)
c     & /zab2(j2,j3,j4,j1)**2*s(j5,j6)
c     &  -2._dp*zb(j3,j4)*zb(j1,j6)*zb(j6,j5)*za(j3,j2)*za(j3,j5)
c     & *za(j6,j5)/zba2(j1,j5,j6,j2)**3*t(j5,j6,j1))

c     & + IDelta * (
c     &  2._dp*zb(j3,j4)*zb(j1,j6)*zb(j6,j5)*za(j3,j2)*za(j3,j5)
c     &  *za(j6,j5)/zba2(j1,j5,j6,j2)**3*t(j5,j6,j2)
c     &  -2._dp*zb(j3,j4)*zb(j1,j6)*za(j3,j2)*za(j3,j5)
c     & /zab2(j2,j3,j4,j1)**3*t(j4,j3,j1)*t(j4,j3,j2)
c     & + 2._dp*zb(j3,j4)*zb(j1,j6)*za(j3,j2)*za(j3,j5)
c     & /zab2(j2,j3,j4,j1)**3*t(j4,j3,j2)**2
c     &  + 2._dp*zb(j3,j4)*zb(j2,j6)*za(j3,j2)*za(j3,j5)*
c     &    izab2(j2,j3,j4,j1)**2*t(j4,j3,j1)
c     &  -2._dp*zb(j3,j4)*zb(j2,j6)*za(j3,j2)*za(j3,j5)
c     & /zab2(j2,j3,j4,j1)**2*t(j4,j3,j2)
c     &  - one/2._dp*zb(j3,j4)*za(j3,j5)**2/za(j6,j5)*zba2(j2,j4,j3,j1)*
c     &    izba2(j1,j3,j4,j2)
c     & - one/2._dp*zb(j3,j4)*za(j3,j5)**2/za(j6,j5)*zba2(j2,j4,j3,j1)
c     &  /zab2(j2,j3,j4,j1)
c     & + zb(j3,j4)*za(j3,j5)**2/za(j6,j5)
c     & /zab2(j2,j3,j4,j1)**2*s(j3,j4)*t(j4,j3,j2)
c     & -zb(j3,j4)*za(j3,j5)**2
c     & /za(j6,j5)/zab2(j2,j3,j4,j1)**2*s(j1,j2)*t(j4,j3,j2)
c     &  -zb(j3,j4)*za(j3,j5)**2/za(j6,j5)/zab2(j2,j3,j4,j1)**2
c     & *s(j5,j6)*t(j4,j3,j2)
c     &  -2._dp*zb(j4,j1)*zb(j4,j6)*zb(j6,j5)*za(j3,j4)*za(j2,j5)
c     & *za(j6,j5)/zba2(j1,j6,j5,j2)**3*t(j6,j5,j1))

c     & + IDelta
c     & * ( 2._dp*zb(j4,j1)*zb(j4,j6)*zb(j6,j5)*za(j3,j4)*za(j2,j5)
c     &   *za(j6,j5)/zba2(j1,j6,j5,j2)**3*t(j6,j5,j2)
c     &  + 2._dp*zb(j4,j1)*zb(j4,j6)*za(j3,j4)*za(j1,j5)
c     &   /zba2(j1,j3,j4,j2)**2*t(j3,j4,j1)
c     &  - 2._dp*zb(j4,j1)*zb(j4,j6)*za(j3,j4)*za(j1,j5)
c     &   /zba2(j1,j3,j4,j2)**2*t(j3,j4,j2)
c     &  - 2._dp*zb(j4,j1)*zb(j4,j6)*za(j3,j4)*za(j2,j5)
c     &  /zba2(j1,j3,j4,j2)**3*t(j3,j4,j1)**2
c     &  +2._dp*zb(j4,j1)*zb(j4,j6)*za(j3,j4)*za(j2,j5)
c     & /zba2(j1,j3,j4,j2)**3*t(j3,j4,j1)*t(j3,j4,j2)
c     &  - 2._dp*zb(j4,j1)*zb(j6,j5)*za(j3,j5)*za(j2,j5)
c     & /zba2(j1,j6,j5,j2)**3*t(j6,j5,j1)*t(j6,j5,j2)
c     &  + 2._dp*zb(j4,j1)*zb(j6,j5)*za(j3,j5)*za(j2,j5)
c     & /zba2(j1,j6,j5,j2)**3*t(j6,j5,j2)**2
c     &  + 2._dp*zb(j4,j2)*zb(j6,j5)*za(j3,j5)*za(j2,j5)
c     & /zba2(j1,j6,j5,j2)**2*t(j6,j5,j1)
c     &  - 2._dp*zb(j4,j2)*zb(j6,j5)*za(j3,j5)*za(j2,j5)
c     & /zba2(j1,j6,j5,j2)**2*t(j6,j5,j2)
c     &  - one/2._dp*zb(j4,j6)**2*za(j3,j4)/zb(j6,j5)*zba2(j2,j4,j3,j1)
c     &  /zba2(j1,j3,j4,j2)
c     &  - one/2._dp*zb(j4,j6)**2*za(j3,j4)/zb(j6,j5)*zba2(j2,j4,j3,j1)
c     & /zab2(j2,j3,j4,j1))
c     & + IDelta
c     & * ( zb(j4,j6)**2*za(j3,j4)/zb(j6,j5)/zba2(j1,j3,j4,j2)**2
c     & *s(j3,j4)*t(j3,j4,j1)
c     &  -zb(j4,j6)**2*za(j3,j4)/zb(j6,j5)
c     & /zba2(j1,j3,j4,j2)**2*s(j1,j2)*t(j3,j4,j1)
c     &  - zb(j4,j6)**2*za(j3,j4)/zb(j6,j5)
c     & /zba2(j1,j3,j4,j2)**2*s(j5,j6)*t(j3,j4,j1)
c     &  -one/2._dp*zb(j4,j6)**2*za(j6,j5)/zb(j3,j4)*zba2(j2,j6,j5,j1)
c     & /zba2(j1,j6,j5,j2)
c     &  - one/2._dp*zb(j4,j6)**2*za(j6,j5)/zb(j3,j4)*zba2(j2,j5,j6,j1)
c     & /zba2(j1,j5,j6,j2)
c     &  - zb(j4,j6)**2*za(j6,j5)/zb(j3,j4)
c     & /zba2(j1,j5,j6,j2)**2*s(j3,j4)*t(j5,j6,j1)
c     &  -zb(j4,j6)**2*za(j6,j5)/zb(j3,j4)
c     & /zba2(j1,j5,j6,j2)**2*s(j1,j2)*t(j5,j6,j1)
c     & + zb(j4,j6)**2*za(j6,j5)/zb(j3,j4)
c     & /zba2(j1,j5,j6,j2)**2*s(j5,j6)*t(j5,j6,j1)
c     & + 2._dp*zb(j4,j6)*zb(j1,j6)*za(j3,j1)*za(j6,j5)
c     &  /zba2(j1,j5,j6,j2)**2*t(j5,j6,j1)
c     &  - 2._dp*zb(j4,j6)*zb(j1,j6)*za(j3,j1)*za(j6,j5)
c     &  /zba2(j1,j5,j6,j2)**2*t(j5,j6,j2)
c     &  - 2._dp*zb(j4,j6)*zb(j1,j6)*za(j3,j2)*za(j6,j5)
c     & /zba2(j1,j5,j6,j2)**3*t(j5,j6,j1)**2)
c     & + IDelta
c     &  * ( 2._dp*zb(j4,j6)*zb(j1,j6)*za(j3,j2)*za(j6,j5)
c     &  /zba2(j1,j5,j6,j2)**3*t(j5,j6,j1)*t(j5,j6,j2)
c     &  -zb(j4,j6)*zb(j6,j5)*za(j3,j5)*za(j6,j5)
c     &  /zba2(j1,j6,j5,j2)**2*s(j3,j4)
c     &  - zb(j4,j6)*zb(j6,j5)*za(j3,j5)*za(j6,j5)
c     &  /zba2(j1,j6,j5,j2)**2*s(j1,j2)
c     &  + zb(j4,j6)*zb(j6,j5)*za(j3,j5)*za(j6,j5)
c     &  /zba2(j1,j6,j5,j2)**2*s(j5,j6)
c     & - zb(j4,j6)*zb(j6,j5)*za(j3,j5)*za(j6,j5)
c     & /zba2(j1,j5,j6,j2)**2*s(j3,j4)
c     &  -zb(j4,j6)*zb(j6,j5)*za(j3,j5)*za(j6,j5)
c     &  /zba2(j1,j5,j6,j2)**2*s(j1,j2)
c     &  +zb(j4,j6)*zb(j6,j5)*za(j3,j5)*za(j6,j5)
c     & /zba2(j1,j5,j6,j2)**2*s(j5,j6)
c     &  +zb(j4,j6)*za(j3,j5)*zab2(j1,j3,j4,j2)/zba2(j1,j3,j4,j2)
c     &  +zb(j4,j6)*za(j3,j5)*zba2(j2,j4,j3,j1)/zab2(j2,j3,j4,j1)
c     &  +zb(j4,j6)*za(j3,j5)*zba2(j2,j6,j5,j1)*izba2(j1,j6,j5,j2)
c     &  +zb(j4,j6)*za(j3,j5)*zba2(j2,j5,j6,j1)*izba2(j1,j5,j6,j2)
c     &  - one/2._dp*zb(j6,j5)*za(j3,j5)**2/za(j3,j4)*zba2(j2,j6,j5,j1)
c     & /zba2(j1,j6,j5,j2)
c     &  - one/2._dp*zb(j6,j5)*za(j3,j5)**2/za(j3,j4)*zba2(j2,j5,j6,j1)
c     & /zba2(j1,j5,j6,j2))

c     & + IDelta
c     & * (-zb(j6,j5)*za(j3,j5)**2
c     & /za(j3,j4)/zba2(j1,j6,j5,j2)**2*s(j3,j4)*t(j6,j5,j2)
c     & -zb(j6,j5)*za(j3,j5)**2
c     & /za(j3,j4)/zba2(j1,j6,j5,j2)**2*s(j1,j2)*t(j6,j5,j2)
c     &  +zb(j6,j5)*za(j3,j5)**2
c     & /za(j3,j4)/zba2(j1,j6,j5,j2)**2*s(j5,j6)*t(j6,j5,j2))

c     & - 2._dp*zb(j4,j1)*zb(j4,j2)*za(j1,j5)*za(j2,j5)
c     &  /zb(j3,j4)/za(j6,j5)/zba2(j1,j3,j4,j2)**2
c     &  + 2._dp*zb(j4,j1)*za(j2,j5)
c     & /zb(j3,j4)/za(j6,j5)*zba2(j4,j3,j2,j5)
c     & /zba2(j1,j3,j4,j2)**3*t(j3,j4,j1)
c     &  -2._dp*zb(j1,j6)*zb(j2,j6)*za(j3,j1)*za(j3,j2)
c     & /zb(j6,j5)/za(j3,j4)/zab2(j2,j3,j4,j1)**2
c     & + 2._dp*zb(j1,j6)*za(j3,j2)/zb(j6,j5)/za(j3,j4)
c     & *zba2(j6,j4,j1,j3)/zab2(j2,j3,j4,j1)**3*t(j4,j3,j2)
c     &  -izb(j3,j4)/za(j6,j5)*zba2(j4,j3,j2,j5)**2/zba2(j1,j3,j4,j2)**2
c     &   -izb(j6,j5)/za(j3,j4)*zba2(j6,j4,j1,j3)**2/zab2(j2,j3,j4,j1)**2

c      write(6,*) 'new:b12zm',bcoeff(b12zm)
c      pause

c      write(6,*) 'old b134',oldfuncb134mp(j1,j2,j3,j4,j5,j6,za,zb)
c      write(6,*) 'NEW b134',funcb134mp(j1,j2,j3,j4,j5,j6,za,zb)
c      write(6,*) 'ratio',funcb134mp(j1,j2,j3,j4,j5,j6,za,zb)
c     &                  /oldfuncb134mp(j1,j2,j3,j4,j5,j6,za,zb)
c      write(6,*)
c      write(6,*) 'old b234',oldfuncb134mp(j2,j1,j4,j3,j6,j5,zb,za)
c      write(6,*) 'NEW b234',funcb134mp(j2,j1,j4,j3,j6,j5,zb,za)
c      write(6,*) 'ratio',funcb134mp(j2,j1,j4,j3,j6,j5,zb,za)
c     &                  /oldfuncb134mp(j2,j1,j4,j3,j6,j5,zb,za)
c      pause

c      write(6,*) 'old b34',oldfuncb34mp(j1,j2,j3,j4,j5,j6,za,zb)
c      write(6,*) 'NEW b34',funcb34mp(j1,j2,j3,j4,j5,j6,za,zb)
c      write(6,*) 'ratio',funcb34mp(j1,j2,j3,j4,j5,j6,za,zb)
c     &                  /oldfuncb34mp(j1,j2,j3,j4,j5,j6,za,zb)
c      write(6,*)
c      write(6,*) 'old b56',oldfuncb34mp(j2,j1,j4,j3,j6,j5,zb,za)
c      write(6,*) 'NEW b56',funcb34mp(j2,j1,j4,j3,j6,j5,zb,za)
c      write(6,*) 'ratio',funcb34mp(j2,j1,j4,j3,j6,j5,zb,za)
c     &                  /oldfuncb34mp(j2,j1,j4,j3,j6,j5,zb,za)
c      pause

      if (ggZZuse6d) then
      bcoeff(b34)=funcb34mp(j1,j2,j3,j4,j5,j6,za,zb)
      bcoeff(b56)=funcb34mp(j1,j2,j5,j6,j3,j4,za,zb)
      bcoeff(b134)=funcb134mp(j1,j2,j3,j4,j5,j6,za,zb)
      bcoeff(b234)=funcb134mp(j2,j1,j4,j3,j6,j5,zb,za)
      else
      bcoeff(b34)=oldfuncb34mp(j1,j2,j3,j4,j5,j6,za,zb)
      bcoeff(b56)=oldfuncb34mp(j1,j2,j5,j6,j3,j4,za,zb)
      bcoeff(b134)=oldfuncb134mp(j1,j2,j3,j4,j5,j6,za,zb)
      bcoeff(b234)=oldfuncb134mp(j2,j1,j4,j3,j6,j5,zb,za)
      endif
      bcoeff(b12)=-bcoeff(b34)-bcoeff(b56)
     & -bcoeff(b134)-bcoeff(b234)

      bcoeff(rat)=
     & ((zb(j1,j4)*za(j1,j5))**2*za(j4,j3)
     & /(za(j6,j5)*(s(j1,j4)+s(j1,j3)))
     & +(zb(j2,j6)*za(j2,j3))**2*zb(j4,j3)
     & /(zb(j6,j5)*(s(j2,j4)+s(j2,j3)))
     & +(zb(j1,j6)*za(j1,j3))**2*za(j6,j5)
     & /(za(j4,j3)*(s(j1,j6)+s(j1,j5)))
     & +(zb(j2,j4)*za(j2,j5))**2*zb(j6,j5)
     & /(zb(j4,j3)*(s(j2,j6)+s(j2,j5)))
     & -zab2(j5,j2,j3,j4)**2/zb(j4,j3)/za(j6,j5)
     & -zab2(j3,j1,j4,j6)**2/za(j4,j3)/zb(j6,j5)
     & -two*zb(j4,j6)*za(j3,j5))/zab2(j2,j3,j4,j1)**2
     & -IDelta*zab2(j1,j3,j4,j2)/zab2(j2,j3,j4,j1)*(
     & +four*zb(j4,j6)*za(j3,j5)
     & +(s(j1,j2)-s(j3,j4)-s(j5,j6))
     & *(zb(j4,j6)**2/zb(j4,j3)/zb(j6,j5)
     &  +za(j3,j5)**2/za(j4,j3)/za(j6,j5)))


      endif

      return
      end



c--- This form has been massaged to remove the explicit <2|3+4|1]^3 poles
      function funcb34mp(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: funcb34mp
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t,sbar134,sbar234,del,s134,s234
      complex(dp):: zab2
      real(dp):: IDelta
c---statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      t(j1,j2,j3)=s(j1,j2)+s(j2,j3)+s(j1,j3)
      del(j1,j2,j3,j4,j5,j6)=s(j1,j2)-s(j3,j4)-s(j5,j6)
c---end statement functions

      IDelta=s(j1,j2)**2+s(j3,j4)**2+s(j5,j6)**2
     & -two*s(j1,j2)*s(j3,j4)
     & -two*s(j3,j4)*s(j5,j6)
     & -two*s(j1,j2)*s(j5,j6)

      IDelta=one/IDelta
      s134=t(j1,j3,j4)
      s234=t(j2,j3,j4)
      sbar134=s134-s(j3,j4)
      sbar234=s234-s(j3,j4)


      funcb34mp=
     & -za(j4,j3)*za(j1,j5)**2*zb(j1,j4)**2*s(j3,j4)
     & /sbar134**2/za(j6,j5)/zab2(j2,j3,j4,j1)**2

     & -zb(j4,j3)*zb(j2,j6)**2*za(j2,j3)**2*s(j3,j4)
     & /sbar234**2/zb(j6,j5)/zab2(j2,j3,j4,j1)**2

     & -two*za(j3,j2)*zb(j1,j4)*s(j3,j4)
     & /s(j5,j6)/zab2(j2,j3,j4,j1)**3/sbar134/sbar234
     & *(s(j3,j4)*(za(j2,j5)*zb(j1,j6)*zab2(j1,j3,j4,j2)
     &            -za(j1,j5)*zb(j2,j6)*zab2(j2,j3,j4,j1))
     &  -za(j1,j5)*zb(j5,j6)*zab2(j5,j3,j4,j2)*zab2(j2,j3,j4,j1)
     &  +za(j5,j6)*zb(j2,j6)*zab2(j1,j3,j4,j6)*zab2(j2,j3,j4,j1))

     & +6._dp*IDelta**2
     & *zab2(j3,j1,j2,j4)*zab2(j5,j1,j2,j6)*del(j5,j6,j1,j2,j3,j4)
     & *zab2(j1,j3,j4,j2)/zab2(j2,j3,j4,j1)

     & + IDelta *(
     & +two*(s134-s234)/s(j5,j6)*zab2(j1,j3,j4,j2)/zab2(j2,j3,j4,j1)*(
     &  +s234*zb(j1,j4)*zb(j1,j6)*za(j2,j3)*za(j2,j5)
     &       /zab2(j2,j3,j4,j1)**2
     &  +s234*zb(j1,j4)*za(j1,j5)*zab2(j3,j2,j5,j6)
     &       /zab2(j1,j3,j4,j2)/zab2(j2,j3,j4,j1)
     &  +zb(j4,j6)*za(j2,j5)*zab2(j3,j2,j4,j1)/zab2(j2,j3,j4,j1)
     &  -s134*za(j2,j3)*za(j2,j5)*zb(j1,j4)*zb(j1,j6)
     &       /zab2(j2,j3,j4,j1)**2
     &  -s134*za(j2,j3)*zb(j2,j6)*zab2(j5,j1,j6,j4)
     &       /zab2(j1,j3,j4,j2)/zab2(j2,j3,j4,j1)
     &  -za(j3,j5)*zb(j1,j6)*zab2(j2,j1,j3,j4)/zab2(j2,j3,j4,j1) )

     & -eight*za(j1,j5)*za(j2,j5)*zb(j1,j4)**2*za(j3,j4)
     &  /za(j5,j6)/zab2(j2,j3,j4,j1)**2*zab2(j1,j3,j4,j2)

     & -eight*za(j2,j3)**2*zb(j3,j4)*zb(j2,j6)*zb(j1,j6)
     &  /zb(j5,j6)/zab2(j2,j3,j4,j1)**2*zab2(j1,j3,j4,j2)

     & +two*za(j3,j4)*za(j1,j5)*zb(j4,j1)*zb(j4,j6)*(s234-s134)
     & /zab2(j2,j3,j4,j1)**2
     & +two*za(j3,j2)*za(j3,j5)*zb(j3,j4)*zb(j2,j6)*(s234-s134)
     & /zab2(j2,j3,j4,j1)**2
     & +two*za(j3,j5)*zb(j4,j6)*s(j3,j4)*del(j3,j4,j1,j2,j5,j6)
     & /zab2(j2,j3,j4,j1)**2
     & -za(j3,j5)**2*zb(j3,j4)*del(j3,j4,j1,j2,j5,j6)*s234
     & /za(j6,j5)/zab2(j2,j3,j4,j1)**2
     & -za(j3,j4)*zb(j4,j6)**2*del(j3,j4,j1,j2,j5,j6)*s134
     & /zb(j6,j5)/zab2(j2,j3,j4,j1)**2

     & +za(j3,j5)**2*zb(j3,j4)*zab2(j1,j3,j4,j2)
     & /za(j6,j5)/zab2(j2,j3,j4,j1)

     & +za(j3,j4)*zb(j4,j6)**2*zab2(j1,j3,j4,j2)
     & /zb(j6,j5)/zab2(j2,j3,j4,j1)

     & -two*za(j3,j5)*zb(j4,j6)*zb(j6,j5)*zab2(j1,j3,j4,j2)
     & /zb(j6,j5)/zab2(j2,j3,j4,j1)
     & )

      return
      end


c--- This form is extracted from BDK
      function oldfuncb34mp(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: oldfuncb34mp
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t,sbar134,sbar234,del,s134,s234
      complex(dp):: zab2
      real(dp):: IDelta
c---statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      t(j1,j2,j3)=s(j1,j2)+s(j2,j3)+s(j1,j3)
      del(j1,j2,j3,j4,j5,j6)=s(j1,j2)-s(j3,j4)-s(j5,j6)
c---end statement functions

      IDelta=s(j1,j2)**2+s(j3,j4)**2+s(j5,j6)**2
     & -two*s(j1,j2)*s(j3,j4)
     & -two*s(j3,j4)*s(j5,j6)
     & -two*s(j1,j2)*s(j5,j6)

      IDelta=one/IDelta
      s134=t(j1,j3,j4)
      s234=t(j2,j3,j4)
      sbar134=s134-s(j3,j4)
      sbar234=s234-s(j3,j4)


      oldfuncb34mp=
     & -za(j4,j3)*za(j1,j5)**2*zb(j1,j4)**2*s(j3,j4)
     & /sbar134**2/za(j6,j5)/zab2(j2,j3,j4,j1)**2

     & -zb(j4,j3)*zb(j2,j6)**2*za(j2,j3)**2*s(j3,j4)
     & /sbar234**2/zb(j6,j5)/zab2(j2,j3,j4,j1)**2

     & -two*za(j1,j5)*za(j3,j2)*zb(j4,j1)*zab2(j5,j3,j4,j1)*s(j3,j4)
     & /sbar134/za(j6,j5)/zab2(j2,j3,j4,j1)**3

     & -two*za(j3,j2)*zb(j2,j6)*zb(j1,j4)*zab2(j2,j3,j4,j6)*s(j3,j4)
     &  /sbar234/zb(j6,j5)/zab2(j2,j3,j4,j1)**3

     & -two*za(j1,j5)*za(j2,j5)*zb(j1,j4)**2*za(j4,j3)
     &  /za(j6,j5)/zab2(j2,j3,j4,j1)**3

     & -two*za(j2,j3)**2*zb(j4,j3)*zb(j2,j6)*zb(j1,j6)
     &  /zb(j6,j5)/zab2(j2,j3,j4,j1)**3

     & +6._dp*IDelta**2
     & *zab2(j3,j1,j2,j4)*zab2(j5,j1,j2,j6)*del(j5,j6,j1,j2,j3,j4)
     & *zab2(j1,j3,j4,j2)/zab2(j2,j3,j4,j1)

     & + IDelta *(
     & +two*za(j3,j5)*za(j2,j5)*zb(j4,j1)*zb(j6,j5)*s(j3,j4)*(s234-s134)
     & /zab2(j2,j3,j4,j1)**3
     & -two*za(j3,j4)*za(j2,j5)*zb(j4,j1)*zb(j4,j6)*s134*(s234-s134)
     & /zab2(j2,j3,j4,j1)**3
     & +two*za(j3,j2)*za(j6,j5)*zb(j4,j6)*zb(j1,j6)*s(j3,j4)*(s234-s134)
     & /zab2(j2,j3,j4,j1)**3
     & -two*za(j3,j2)*za(j3,j5)*zb(j3,j4)*zb(j1,j6)*(s234-s134)*s234
     & /zab2(j2,j3,j4,j1)**3

     & +two*za(j3,j4)*za(j1,j5)*zb(j4,j1)*zb(j4,j6)*(s234-s134)
     & /zab2(j2,j3,j4,j1)**2
     & +two*za(j3,j2)*za(j3,j5)*zb(j3,j4)*zb(j2,j6)*(s234-s134)
     & /zab2(j2,j3,j4,j1)**2
     & +two*za(j3,j5)*zb(j4,j6)*s(j3,j4)*del(j3,j4,j1,j2,j5,j6)
     & /zab2(j2,j3,j4,j1)**2
     & -za(j3,j5)**2*zb(j3,j4)*del(j3,j4,j1,j2,j5,j6)*s234
     & /za(j6,j5)/zab2(j2,j3,j4,j1)**2
     & -za(j3,j4)*zb(j4,j6)**2*del(j3,j4,j1,j2,j5,j6)*s134
     & /zb(j6,j5)/zab2(j2,j3,j4,j1)**2

     & +za(j3,j5)**2*zb(j3,j4)*zab2(j1,j3,j4,j2)
     & /za(j6,j5)/zab2(j2,j3,j4,j1)

     & +za(j3,j4)*zb(j4,j6)**2*zab2(j1,j3,j4,j2)
     & /zb(j6,j5)/zab2(j2,j3,j4,j1)

     & -two*za(j3,j5)*zb(j4,j6)*zb(j6,j5)*zab2(j1,j3,j4,j2)
     & /zb(j6,j5)/zab2(j2,j3,j4,j1)
     & )

      return
      end


c--- This form has been massaged to remove the explicit <2|3+4|1]^3 poles
      function funcb134mp(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: funcb134mp
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: j1,j2,j4,j3,j6,j5
      complex(dp):: zab2,zbb22
c---statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zbb22(j1,j2,j3,j4,j5,j6)
     & =zb(j1,j2)*zab2(j2,j4,j5,j6)+zb(j1,j3)*zab2(j3,j4,j5,j6)
c---end statement functions

      funcb134mp=

     &+ zab2(j1,j3,j4,j2)/zab2(j2,j3,j4,j1)**3 * (
     & + two*za(j2,j5)**2*zb(j1,j4)/za(j1,j2)/za(j5,j6)/zb(j1,j2)
     &   /zb(j3,j4)*zbb22(j4,j5,j6,j3,j4,j1))

     &+ 1._dp/zab2(j2,j3,j4,j1)**2 * (
     & + s(j3,j4)*za(j1,j5)**2*za(j3,j4)*zb(j1,j4)**2
     &   /za(j5,j6)/zab2(j1,j3,j4,j1)**2
     & + s(j5,j6)*za(j2,j5)**2*zb(j2,j4)**2*zb(j5,j6)
     &   /zb(j3,j4)/zab2(j2,j5,j6,j2)**2
     & - two*s(j3,j4)*za(j1,j5)**2*za(j2,j3)*zb(j1,j4)
     &   /za(j1,j2)/za(j5,j6)/zab2(j1,j3,j4,j1)
     & - two*s(j5,j6)*za(j2,j5)*zb(j1,j6)*zb(j2,j4)**2
     &   /zb(j1,j2)/zb(j3,j4)/zab2(j2,j5,j6,j2)
     & + two*za(j1,j5)*za(j2,j5)*zb(j1,j4)*zbb22(j2,j3,j4,j5,j6,j4)
     &   /za(j1,j2)/zb(j1,j2)/zb(j3,j4)/za(j5,j6)
     & - 1._dp/za(j5,j6)/zb(j3,j4)*zab2(j5,j2,j3,j4)**2)


      return
      end


c--- This form is extracted from BDK
      function oldfuncb134mp(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: oldfuncb134mp
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: j1,j2,j4,j3,j6,j5
      real(dp):: t
      complex(dp):: zab2
c---statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      t(j1,j2,j4)=s(j1,j2)+s(j2,j4)+s(j1,j4)
c---end statement functions

      oldfuncb134mp=
     &  +za(j3,j4)*za(j1,j5)**2*zb(j1,j4)**2*s(j3,j4)
     & /((s(j1,j4)+s(j1,j3))**2*za(j5,j6)*zab2(j2,j3,j4,j1)**2)

     & +za(j2,j5)**2*zb(j2,j4)**2*zb(j5,j6)*s(j5,j6)
     & /((s(j2,j6)+s(j2,j5))**2*zb(j3,j4)*zab2(j2,j5,j6,j1)**2)

     & -two*za(j2,j3)*zb(j1,j4)*za(j1,j5)*zab2(j5,j3,j4,j1)*s(j3,j4)
     & /((s(j1,j4)+s(j1,j3))*za(j5,j6)*zab2(j2,j3,j4,j1)**3)

     & +two*za(j2,j5)*zb(j2,j4)*zb(j1,j6)*zab2(j2,j5,j6,j4)*s(j5,j6)
     & /((s(j2,j6)+s(j2,j5))*zb(j3,j4)*zab2(j2,j5,j6,j1)**3)

     & +(two*zb(j1,j4)*za(j2,j5)*(
     & +zb(j1,j4)*za(j1,j5)*t(j2,j4,j3)
     & -zb(j2,j4)*za(j2,j5)*s(j5,j6)
     & +zb(j2,j1)*za(j1,j5)*zb(j4,j3)*za(j3,j2)
     & +t(j1,j4,j3)*zab2(j5,j2,j3,j4))
     & +zab2(j5,j2,j3,j4)**2*zab2(j2,j3,j4,j1))
     & /zb(j3,j4)/za(j6,j5)/zab2(j2,j3,j4,j1)**3

      return
      end

