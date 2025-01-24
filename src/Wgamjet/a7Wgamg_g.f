!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a7Wgamg_g(p1,p2,p3,p4,p5,p6,p7,za,zb,ab,ba)
      implicit none
c     Colour ordered amplitude with a factor of
c     i*e*gw^2*gs^2/(2*rt2) T^CA T^CB removed
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer::p1,p2,p3,p4,p5,p6,p7,h6,h7
      complex(dp):: ab(2,2,2),ba(2,2,2),tmp(2,2,2)
c     amplitude ab(h5,h6,h7),ba(h5,h6,h7)
      call a7Wgamg_gunsym(p1,p2,p3,p4,p5,p6,p7,za,zb,ab)
      call a7Wgamg_gunsym(p1,p2,p3,p4,p5,p7,p6,za,zb,tmp)
      do h7=1,2
      do h6=1,2
      ba(:,h6,h7)=tmp(:,h7,h6)
      enddo
      enddo
      return
      end

      subroutine a7Wgamg_gunsym(p1,p2,p3,p4,p5,p6,p7,za,zb,a)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      integer::p1,p2,p3,p4,p5,p6,p7
      real(dp),parameter::Qsum=1._dp/3._dp
      real(dp):: s3,s134,s234,s157,s167,s256,s267,s34,s345
      complex(dp):: iza,izb,zab2,zab3,a(2,2,2),prop34,prop345
c      amplitude a(h5,h6,h7)

      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p1,p3)
      zab2(p1,p2,p3,p4)=+za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zab3(p1,p2,p3,p4,p5)=
     & +za(p1,p2)*zb(p2,p5)+za(p1,p3)*zb(p3,p5)+za(p1,p4)*zb(p4,p5)
      iza(p1,p2)=cone/za(p1,p2)
      izb(p1,p2)=cone/zb(p1,p2)

      s134=s3(p1,p3,p4)
      s157=s3(p1,p5,p7)
      s167=s3(p1,p6,p7)
      s234=s3(p2,p3,p4)
      s256=s3(p2,p5,p6)
      s267=s3(p2,p6,p7)
      s34=s(p3,p4)
      s345=s3(p3,p4,p5)
      prop34=cmplx(s34-wmass**2,wmass*wwidth,kind=dp)
      prop345=cmplx(s345-wmass**2,wmass*wwidth,kind=dp)

      a(1,1,1)= + prop34**(-1)*prop345**(-1)*s267**(-1) * (  - za(p2,p6
     &    )*za(p3,p5)*zb(p1,p4)**2*iza(p6,p7)*izb(p1,p5)*izb(p1,p6)*
     &    izb(p1,p7)*izb(p6,p7)*zab2(p7,p6,p7,p1)*zab3(p4,p2,p6,p7,p1)
     &     + za(p2,p6)*za(p3,p5)*zb(p1,p4)*iza(p6,p7)*izb(p1,p6)*izb(p1
     &    ,p7)*izb(p6,p7)*zab2(p7,p6,p7,p1)*zab3(p5,p2,p6,p7,p1) + za(
     &    p2,p6)*za(p4,p5)*zb(p1,p4)**2*iza(p6,p7)*izb(p1,p5)*izb(p1,p6
     &    )*izb(p1,p7)*izb(p6,p7)*zab2(p7,p6,p7,p1)*zab3(p3,p2,p6,p7,p1
     &    ) + za(p2,p6)*zb(p1,p4)*iza(p6,p7)*izb(p1,p5)*izb(p1,p6)*izb(
     &    p1,p7)*izb(p6,p7)*zab2(p5,p3,p4,p1)*zab2(p7,p6,p7,p1)*zab3(p3
     &    ,p2,p6,p7,p1) - za(p2,p7)*za(p3,p5)*zb(p1,p4)**2*izb(p1,p5)*
     &    izb(p1,p6)*izb(p6,p7)*zab3(p4,p2,p6,p7,p1) + za(p2,p7)*za(p3,
     &    p5)*zb(p1,p4)*izb(p1,p6)*izb(p6,p7)*zab3(p5,p2,p6,p7,p1) +
     &    za(p2,p7)*za(p4,p5)*zb(p1,p4)**2*izb(p1,p5)*izb(p1,p6)*izb(p6
     &    ,p7)*zab3(p3,p2,p6,p7,p1) + za(p2,p7)*zb(p1,p4)*izb(p1,p5)*
     &    izb(p1,p6)*izb(p6,p7)*zab2(p5,p3,p4,p1)*zab3(p3,p2,p6,p7,p1)
     &     - za(p3,p5)*zb(p1,p4)**2*izb(p1,p5)*izb(p1,p6)*izb(p1,p7)*
     &    izb(p2,p6)*zab2(p7,p2,p6,p1)*zab3(p4,p2,p6,p7,p1) )
      a(1,1,1) = a(1,1,1) + prop34**(-1)*prop345**(-1)*s267**(-1) * (
     &    za(p3,p5)*zb(p1,p4)*izb(p1,p6)*izb(p1,p7)*izb(p2,p6)*zab2(p7,
     &    p2,p6,p1)*zab3(p5,p2,p6,p7,p1) + za(p4,p5)*zb(p1,p4)**2*izb(
     &    p1,p5)*izb(p1,p6)*izb(p1,p7)*izb(p2,p6)*zab2(p7,p2,p6,p1)*
     &    zab3(p3,p2,p6,p7,p1) + zb(p1,p4)*izb(p1,p5)*izb(p1,p6)*izb(p1
     &    ,p7)*izb(p2,p6)*zab2(p5,p3,p4,p1)*zab2(p7,p2,p6,p1)*zab3(p3,
     &    p2,p6,p7,p1) )
      a(1,1,1) = a(1,1,1) + prop34**(-1)*s267**(-1) * ( za(p2,p6)*zb(p1
     &    ,p4)*iza(p6,p7)*izb(p1,p5)*izb(p1,p6)*izb(p1,p7)*izb(p6,p7)*
     &    zab2(p3,p1,p4,p1)*zab2(p7,p6,p7,p1)*zab3(p5,p2,p6,p7,p1)*
     &    s134**(-1) + za(p2,p7)*zb(p1,p4)*izb(p1,p5)*izb(p1,p6)*izb(p6
     &    ,p7)*zab2(p3,p1,p4,p1)*zab3(p5,p2,p6,p7,p1)*s134**(-1) + zb(
     &    p1,p4)*izb(p1,p5)*izb(p1,p6)*izb(p1,p7)*izb(p2,p6)*zab2(p3,p1
     &    ,p4,p1)*zab2(p7,p2,p6,p1)*zab3(p5,p2,p6,p7,p1)*s134**(-1) )
      a(1,1,1) = a(1,1,1) + prop34**(-1)*s267**(-1)*Qsum * ( za(p2,p6)*
     &    zb(p1,p4)*iza(p6,p7)*izb(p1,p5)*izb(p1,p6)*izb(p1,p7)*izb(p6,
     &    p7)*zab2(p3,p1,p4,p1)*zab2(p7,p6,p7,p1)*zab3(p5,p2,p6,p7,p1)*
     &    s134**(-1) + za(p2,p7)*zb(p1,p4)*izb(p1,p5)*izb(p1,p6)*izb(p6
     &    ,p7)*zab2(p3,p1,p4,p1)*zab3(p5,p2,p6,p7,p1)*s134**(-1) + zb(
     &    p1,p4)*izb(p1,p5)*izb(p1,p6)*izb(p1,p7)*izb(p2,p6)*zab2(p3,p1
     &    ,p4,p1)*zab2(p7,p2,p6,p1)*zab3(p5,p2,p6,p7,p1)*s134**(-1) )
      a(1,1,1) = a(1,1,1) + prop34**(-1) * (  - zb(p1,p4)*iza(p6,p7)*
     &    izb(p1,p5)*izb(p1,p6)*izb(p1,p7)*izb(p2,p5)*izb(p6,p7)*zab2(
     &    p3,p1,p4,p1)*zab2(p6,p2,p5,p1)*zab2(p7,p6,p7,p1)*s134**(-1)
     &     + zb(p1,p4)*izb(p1,p5)*izb(p1,p6)*izb(p1,p7)*izb(p2,p5)*
     &    zab2(p3,p1,p4,p1)*zab2(p6,p2,p5,p1)*zab3(p7,p2,p5,p6,p1)*
     &    s134**(-1)*s256**(-1) + zb(p1,p4)*izb(p1,p5)*izb(p1,p6)*izb(
     &    p1,p7)*izb(p2,p6)*zab2(p3,p1,p4,p1)*zab2(p5,p2,p6,p1)*zab3(p7
     &    ,p2,p5,p6,p1)*s134**(-1)*s256**(-1) - zb(p1,p4)*izb(p1,p5)*
     &    izb(p1,p6)*izb(p2,p5)*izb(p6,p7)*zab2(p3,p1,p4,p1)*zab2(p7,p2
     &    ,p5,p1)*s134**(-1) )
      a(1,1,1) = a(1,1,1) + prop34**(-1)*Qsum * (  - zb(p1,p4)*iza(p6,
     &    p7)*izb(p1,p5)*izb(p1,p6)*izb(p1,p7)*izb(p2,p5)*izb(p6,p7)*
     &    zab2(p3,p1,p4,p1)*zab2(p6,p2,p5,p1)*zab2(p7,p6,p7,p1)*
     &    s134**(-1) + zb(p1,p4)*izb(p1,p5)*izb(p1,p6)*izb(p1,p7)*izb(
     &    p2,p5)*zab2(p3,p1,p4,p1)*zab2(p6,p2,p5,p1)*zab3(p7,p2,p5,p6,
     &    p1)*s134**(-1)*s256**(-1) + zb(p1,p4)*izb(p1,p5)*izb(p1,p6)*
     &    izb(p1,p7)*izb(p2,p6)*zab2(p3,p1,p4,p1)*zab2(p5,p2,p6,p1)*
     &    zab3(p7,p2,p5,p6,p1)*s134**(-1)*s256**(-1) - zb(p1,p4)*izb(p1
     &    ,p5)*izb(p1,p6)*izb(p2,p5)*izb(p6,p7)*zab2(p3,p1,p4,p1)*zab2(
     &    p7,p2,p5,p1)*s134**(-1) )
      a(1,1,1) = a(1,1,1) + prop345**(-1)*s267**(-1) * (  - 2*za(p2,p6)
     &    *zb(p1,p3)*zb(p1,p4)*iza(p6,p7)*izb(p1,p5)*izb(p1,p6)*izb(p1,
     &    p7)*izb(p3,p5)*izb(p6,p7)*zab2(p7,p6,p7,p1)*zab3(p3,p2,p6,p7,
     &    p1) - 2*za(p2,p6)*zb(p1,p4)*iza(p6,p7)*izb(p1,p6)*izb(p1,p7)*
     &    izb(p3,p5)*izb(p6,p7)*zab2(p7,p6,p7,p1)*zab3(p5,p2,p6,p7,p1)
     &     - 2*za(p2,p7)*zb(p1,p3)*zb(p1,p4)*izb(p1,p5)*izb(p1,p6)*izb(
     &    p3,p5)*izb(p6,p7)*zab3(p3,p2,p6,p7,p1) - 2*za(p2,p7)*zb(p1,p4
     &    )*izb(p1,p6)*izb(p3,p5)*izb(p6,p7)*zab3(p5,p2,p6,p7,p1) - 2*
     &    zb(p1,p3)*zb(p1,p4)*izb(p1,p5)*izb(p1,p6)*izb(p1,p7)*izb(p2,
     &    p6)*izb(p3,p5)*zab2(p7,p2,p6,p1)*zab3(p3,p2,p6,p7,p1) - 2*zb(
     &    p1,p4)*izb(p1,p6)*izb(p1,p7)*izb(p2,p6)*izb(p3,p5)*zab2(p7,p2
     &    ,p6,p1)*zab3(p5,p2,p6,p7,p1) )
      a(1,1,2)= + prop34**(-1)*prop345**(-1)*s167**(-1) * (  - za(p1,p6
     &    )*za(p2,p3)*za(p3,p5)*zb(p1,p7)*zb(p3,p7)*iza(p1,p7)*iza(p6,
     &    p7)*izb(p5,p7)*izb(p6,p7)*zab2(p6,p1,p7,p4) + za(p1,p6)*za(p2
     &    ,p3)*za(p3,p5)*zb(p1,p7)*zb(p4,p7)*iza(p1,p7)*iza(p6,p7)*izb(
     &    p5,p7)*izb(p6,p7)*zab2(p6,p1,p7,p3) - za(p1,p6)*za(p2,p3)*za(
     &    p4,p5)*zb(p1,p7)*zb(p4,p7)*iza(p1,p7)*iza(p6,p7)*izb(p5,p7)*
     &    izb(p6,p7)*zab2(p6,p1,p7,p4) + za(p1,p6)*za(p2,p3)*zb(p1,p7)*
     &    iza(p1,p7)*iza(p6,p7)*izb(p5,p7)*izb(p6,p7)*zab2(p5,p3,p4,p7)
     &    *zab2(p6,p1,p7,p4) + za(p1,p6)*za(p2,p4)*za(p3,p5)*zb(p1,p7)*
     &    zb(p4,p7)*iza(p1,p7)*iza(p6,p7)*izb(p5,p7)*izb(p6,p7)*zab2(p6
     &    ,p1,p7,p4) + 2*za(p1,p6)*za(p2,p5)*za(p3,p5)*zb(p1,p7)*zb(p4,
     &    p5)*iza(p1,p7)*iza(p6,p7)*izb(p5,p7)*izb(p6,p7)*zab2(p6,p1,p7
     &    ,p7) - za(p1,p6)*za(p2,p5)*za(p3,p5)*zb(p1,p7)*zb(p4,p7)*iza(
     &    p1,p7)*iza(p6,p7)*izb(p5,p7)*izb(p6,p7)*zab2(p6,p1,p7,p5) )
      a(1,1,2) = a(1,1,2) + prop34**(-1)*prop345**(-1)*s267**(-1) * (
     &     - za(p3,p5)*zb(p1,p3)*zb(p4,p7)*iza(p6,p7)*izb(p2,p6)*izb(p5
     &    ,p7)*izb(p6,p7)*zab2(p3,p2,p6,p7)*zab2(p6,p2,p6,p7) + za(p3,
     &    p5)*zb(p1,p4)*zb(p3,p7)*iza(p6,p7)*izb(p2,p6)*izb(p5,p7)*izb(
     &    p6,p7)*zab2(p3,p2,p6,p7)*zab2(p6,p2,p6,p7) - za(p3,p5)*zb(p1,
     &    p4)*zb(p4,p7)*iza(p6,p7)*izb(p2,p6)*izb(p5,p7)*izb(p6,p7)*
     &    zab2(p4,p2,p6,p7)*zab2(p6,p2,p6,p7) + za(p3,p5)*zb(p1,p5)*zb(
     &    p4,p7)*iza(p6,p7)*izb(p2,p6)*izb(p5,p7)*izb(p6,p7)*zab2(p5,p2
     &    ,p6,p7)*zab2(p6,p2,p6,p7) - 2*za(p3,p5)*zb(p1,p7)*zb(p4,p5)*
     &    iza(p6,p7)*izb(p2,p6)*izb(p5,p7)*izb(p6,p7)*zab2(p5,p2,p6,p7)
     &    *zab2(p6,p2,p6,p7) + za(p4,p5)*zb(p1,p4)*zb(p4,p7)*iza(p6,p7)
     &    *izb(p2,p6)*izb(p5,p7)*izb(p6,p7)*zab2(p3,p2,p6,p7)*zab2(p6,
     &    p2,p6,p7) - zb(p1,p4)*iza(p6,p7)*izb(p2,p6)*izb(p5,p7)*izb(p6
     &    ,p7)*zab2(p3,p2,p6,p7)*zab2(p5,p3,p4,p7)*zab2(p6,p2,p6,p7) )
      a(1,1,2) = a(1,1,2) + prop34**(-1)*prop345**(-1) * ( 2*za(p1,p6)*
     &    za(p3,p5)*zb(p1,p7)*zb(p4,p5)*iza(p1,p7)*iza(p6,p7)*izb(p2,p6
     &    )*izb(p5,p7)*izb(p6,p7)*zab2(p5,p2,p6,p7) + za(p3,p5)*zb(p3,
     &    p7)*iza(p1,p7)*iza(p6,p7)*izb(p2,p6)*izb(p5,p7)*izb(p6,p7)*
     &    zab2(p3,p2,p6,p7)*zab2(p6,p1,p7,p4) - za(p3,p5)*zb(p4,p7)*
     &    iza(p1,p7)*iza(p6,p7)*izb(p2,p6)*izb(p5,p7)*izb(p6,p7)*zab2(
     &    p3,p2,p6,p7)*zab2(p6,p1,p7,p3) - za(p3,p5)*zb(p4,p7)*iza(p1,
     &    p7)*iza(p6,p7)*izb(p2,p6)*izb(p5,p7)*izb(p6,p7)*zab2(p4,p2,p6
     &    ,p7)*zab2(p6,p1,p7,p4) + za(p3,p5)*zb(p4,p7)*iza(p1,p7)*iza(
     &    p6,p7)*izb(p2,p6)*izb(p5,p7)*izb(p6,p7)*zab2(p5,p2,p6,p7)*
     &    zab2(p6,p1,p7,p5) + za(p4,p5)*zb(p4,p7)*iza(p1,p7)*iza(p6,p7)
     &    *izb(p2,p6)*izb(p5,p7)*izb(p6,p7)*zab2(p3,p2,p6,p7)*zab2(p6,
     &    p1,p7,p4) - iza(p1,p7)*iza(p6,p7)*izb(p2,p6)*izb(p5,p7)*izb(
     &    p6,p7)*zab2(p3,p2,p6,p7)*zab2(p5,p3,p4,p7)*zab2(p6,p1,p7,p4)
     &     )
      a(1,1,2) = a(1,1,2) + prop34**(-1)*s167**(-1) * ( za(p1,p6)*za(p2
     &    ,p3)*zb(p1,p7)*iza(p1,p7)*iza(p6,p7)*izb(p5,p7)*izb(p6,p7)*
     &    zab2(p5,p2,p3,p4)*zab2(p6,p1,p7,p7)*s234**(-1) + za(p1,p6)*
     &    zb(p1,p7)*iza(p1,p7)*iza(p6,p7)*izb(p2,p5)*izb(p5,p7)*izb(p6,
     &    p7)*zab2(p3,p2,p5,p7)*zab2(p6,p1,p7,p4) )
      a(1,1,2) = a(1,1,2) + prop34**(-1)*s167**(-1)*Qsum * (  - za(p1,
     &    p6)*za(p2,p3)*zb(p1,p7)*iza(p1,p7)*iza(p6,p7)*izb(p5,p7)*izb(
     &    p6,p7)*zab2(p5,p2,p3,p4)*zab2(p6,p1,p7,p7)*s234**(-1) + za(p1
     &    ,p6)*zb(p1,p7)*iza(p1,p7)*iza(p6,p7)*izb(p2,p5)*izb(p5,p7)*
     &    izb(p6,p7)*zab2(p3,p2,p5,p7)*zab2(p6,p1,p7,p4) )
      a(1,1,2) = a(1,1,2) + prop34**(-1)*s267**(-1) * (  - zb(p1,p4)*
     &    iza(p6,p7)*izb(p2,p6)*izb(p5,p7)*izb(p6,p7)*zab2(p3,p1,p4,p7)
     &    *zab2(p5,p2,p6,p7)*zab2(p6,p2,p6,p7)*s134**(-1) - zb(p1,p7)*
     &    iza(p1,p5)*iza(p6,p7)*izb(p1,p5)*izb(p2,p6)*izb(p5,p7)*izb(p6
     &    ,p7)*zab2(p3,p2,p6,p7)*zab2(p5,p1,p5,p4)*zab2(p6,p2,p6,p7) )
      a(1,1,2) = a(1,1,2) + prop34**(-1)*s267**(-1)*Qsum * (  - zb(p1,
     &    p4)*iza(p6,p7)*izb(p2,p6)*izb(p5,p7)*izb(p6,p7)*zab2(p3,p1,p4
     &    ,p7)*zab2(p5,p2,p6,p7)*zab2(p6,p2,p6,p7)*s134**(-1) + zb(p1,
     &    p7)*iza(p1,p5)*iza(p6,p7)*izb(p1,p5)*izb(p2,p6)*izb(p5,p7)*
     &    izb(p6,p7)*zab2(p3,p2,p6,p7)*zab2(p5,p1,p5,p4)*zab2(p6,p2,p6,
     &    p7) )
      a(1,1,2) = a(1,1,2) + prop34**(-1) * ( za(p1,p6)*za(p2,p3)*zb(p1,
     &    p7)*iza(p1,p7)*iza(p6,p7)*izb(p5,p7)*izb(p6,p7)*zab2(p5,p1,p7
     &    ,p7)*zab2(p6,p2,p3,p4)*s234**(-1)*s157**(-1) - za(p1,p6)*zb(
     &    p1,p7)*iza(p1,p7)*iza(p6,p7)*izb(p2,p6)*izb(p5,p7)*izb(p6,p7)
     &    *zab2(p3,p2,p6,p7)*zab2(p5,p1,p7,p4)*s157**(-1) - za(p2,p3)*
     &    zb(p1,p7)*iza(p1,p5)*iza(p6,p7)*izb(p1,p5)*izb(p5,p7)*izb(p6,
     &    p7)*zab2(p5,p1,p5,p7)*zab2(p6,p1,p5,p7)*zab2(p6,p2,p3,p4)*
     &    s234**(-1)*s157**(-1) - zb(p1,p4)*iza(p6,p7)*izb(p2,p5)*izb(
     &    p5,p7)*izb(p6,p7)*zab2(p3,p1,p4,p7)*zab2(p6,p2,p5,p7)**2*
     &    s134**(-1)*s256**(-1) - zb(p1,p4)*iza(p6,p7)*izb(p2,p6)*izb(
     &    p5,p7)*izb(p6,p7)*zab2(p3,p1,p4,p7)*zab2(p5,p2,p6,p7)*zab2(p6
     &    ,p2,p5,p7)*s134**(-1)*s256**(-1) + zb(p1,p7)*iza(p1,p5)*iza(
     &    p6,p7)*izb(p1,p5)*izb(p2,p6)*izb(p5,p7)*izb(p6,p7)*zab2(p3,p2
     &    ,p6,p7)*zab2(p5,p1,p5,p7)*zab3(p6,p1,p5,p7,p4)*s157**(-1) +
     &    iza(p1,p7)*iza(p6,p7)*izb(p2,p5)*izb(p5,p7)*izb(p6,p7)*zab2(
     &    p6,p1,p7,p4)*zab2(p6,p2,p5,p7)*zab3(p3,p2,p5,p6,p7)*
     &    s256**(-1) )
      a(1,1,2) = a(1,1,2) + prop34**(-1) * ( iza(p1,p7)*iza(p6,p7)*izb(
     &    p2,p6)*izb(p5,p7)*izb(p6,p7)*zab2(p5,p2,p6,p7)*zab2(p6,p1,p7,
     &    p4)*zab3(p3,p2,p5,p6,p7)*s256**(-1) )
      a(1,1,2) = a(1,1,2) + prop34**(-1)*Qsum * (  - za(p1,p6)*za(p2,p3
     &    )*zb(p1,p7)*iza(p1,p7)*iza(p6,p7)*izb(p5,p7)*izb(p6,p7)*zab2(
     &    p5,p1,p7,p7)*zab2(p6,p2,p3,p4)*s234**(-1)*s157**(-1) + za(p1,
     &    p6)*zb(p1,p7)*iza(p1,p7)*iza(p6,p7)*izb(p2,p6)*izb(p5,p7)*
     &    izb(p6,p7)*zab2(p3,p2,p6,p7)*zab2(p5,p1,p7,p4)*s157**(-1) +
     &    za(p2,p3)*zb(p1,p7)*iza(p1,p5)*iza(p6,p7)*izb(p1,p5)*izb(p5,
     &    p7)*izb(p6,p7)*zab2(p5,p1,p5,p7)*zab2(p6,p1,p5,p7)*zab2(p6,p2
     &    ,p3,p4)*s234**(-1)*s157**(-1) - zb(p1,p4)*iza(p6,p7)*izb(p2,
     &    p5)*izb(p5,p7)*izb(p6,p7)*zab2(p3,p1,p4,p7)*zab2(p6,p2,p5,p7)
     &    **2*s134**(-1)*s256**(-1) - zb(p1,p4)*iza(p6,p7)*izb(p2,p6)*
     &    izb(p5,p7)*izb(p6,p7)*zab2(p3,p1,p4,p7)*zab2(p5,p2,p6,p7)*
     &    zab2(p6,p2,p5,p7)*s134**(-1)*s256**(-1) - zb(p1,p7)*iza(p1,p5
     &    )*iza(p6,p7)*izb(p1,p5)*izb(p2,p6)*izb(p5,p7)*izb(p6,p7)*
     &    zab2(p3,p2,p6,p7)*zab2(p5,p1,p5,p7)*zab3(p6,p1,p5,p7,p4)*
     &    s157**(-1) + iza(p1,p7)*iza(p6,p7)*izb(p2,p5)*izb(p5,p7)*izb(
     &    p6,p7)*zab2(p6,p1,p7,p4)*zab2(p6,p2,p5,p7)*zab3(p3,p2,p5,p6,
     &    p7)*s256**(-1) )
      a(1,1,2) = a(1,1,2) + prop34**(-1)*Qsum * ( iza(p1,p7)*iza(p6,p7)
     &    *izb(p2,p6)*izb(p5,p7)*izb(p6,p7)*zab2(p5,p2,p6,p7)*zab2(p6,
     &    p1,p7,p4)*zab3(p3,p2,p5,p6,p7)*s256**(-1) )
      a(1,1,2) = a(1,1,2) + prop345**(-1)*s167**(-1) * ( 2*za(p1,p6)*
     &    za(p2,p3)*zb(p1,p7)*zb(p3,p7)*iza(p1,p7)*iza(p6,p7)*izb(p3,p5
     &    )*izb(p5,p7)*izb(p6,p7)*zab2(p6,p1,p7,p4) + 2*za(p1,p6)*za(p2
     &    ,p5)*zb(p1,p7)*iza(p1,p7)*iza(p6,p7)*izb(p3,p5)*izb(p6,p7)*
     &    zab2(p6,p1,p7,p4) )
      a(1,1,2) = a(1,1,2) + prop345**(-1)*s267**(-1) * (  - 2*zb(p1,p4)
     &    *zb(p3,p7)*iza(p6,p7)*izb(p2,p6)*izb(p3,p5)*izb(p5,p7)*izb(p6
     &    ,p7)*zab2(p3,p2,p6,p7)*zab2(p6,p2,p6,p7) - 2*zb(p1,p4)*iza(p6
     &    ,p7)*izb(p2,p6)*izb(p3,p5)*izb(p6,p7)*zab2(p5,p2,p6,p7)*zab2(
     &    p6,p2,p6,p7) )
      a(1,1,2) = a(1,1,2) + prop345**(-1) * (  - 2*zb(p3,p7)*iza(p1,p7)
     &    *iza(p6,p7)*izb(p2,p6)*izb(p3,p5)*izb(p5,p7)*izb(p6,p7)*zab2(
     &    p3,p2,p6,p7)*zab2(p6,p1,p7,p4) - 2*iza(p1,p7)*iza(p6,p7)*izb(
     &    p2,p6)*izb(p3,p5)*izb(p6,p7)*zab2(p5,p2,p6,p7)*zab2(p6,p1,p7,
     &    p4) )
      a(1,2,1)= + prop34**(-1)*prop345**(-1)*s167**(-1) * ( za(p2,p3)*
     &    za(p3,p5)*zb(p1,p6)*zb(p3,p7)*iza(p1,p7)*iza(p6,p7)*izb(p1,p7
     &    )*izb(p5,p7)*izb(p6,p7)*zab2(p7,p1,p6,p4)*zab2(p7,p1,p7,p6)
     &     - za(p2,p3)*za(p3,p5)*zb(p1,p6)*zb(p3,p7)*iza(p6,p7)**2*izb(
     &    p5,p7)*izb(p6,p7)**2*zab2(p7,p1,p6,p4)*zab2(p7,p6,p7,p6) -
     &    za(p2,p3)*za(p3,p5)*zb(p1,p6)*zb(p4,p7)*iza(p1,p7)*iza(p6,p7)
     &    *izb(p1,p7)*izb(p5,p7)*izb(p6,p7)*zab2(p7,p1,p6,p3)*zab2(p7,
     &    p1,p7,p6) + za(p2,p3)*za(p3,p5)*zb(p1,p6)*zb(p4,p7)*iza(p6,p7
     &    )**2*izb(p5,p7)*izb(p6,p7)**2*zab2(p7,p1,p6,p3)*zab2(p7,p6,p7
     &    ,p6) + za(p2,p3)*za(p4,p5)*zb(p1,p6)*zb(p4,p7)*iza(p1,p7)*
     &    iza(p6,p7)*izb(p1,p7)*izb(p5,p7)*izb(p6,p7)*zab2(p7,p1,p6,p4)
     &    *zab2(p7,p1,p7,p6) - za(p2,p3)*za(p4,p5)*zb(p1,p6)*zb(p4,p7)*
     &    iza(p6,p7)**2*izb(p5,p7)*izb(p6,p7)**2*zab2(p7,p1,p6,p4)*
     &    zab2(p7,p6,p7,p6) - za(p2,p3)*zb(p1,p6)*iza(p1,p7)*iza(p6,p7)
     &    *izb(p1,p7)*izb(p5,p7)*izb(p6,p7)*zab2(p5,p3,p4,p7)*zab2(p7,
     &    p1,p6,p4)*zab2(p7,p1,p7,p6) )
      a(1,2,1) = a(1,2,1) + prop34**(-1)*prop345**(-1)*s167**(-1) * (
     &    za(p2,p3)*zb(p1,p6)*iza(p6,p7)**2*izb(p5,p7)*izb(p6,p7)**2*
     &    zab2(p5,p3,p4,p7)*zab2(p7,p1,p6,p4)*zab2(p7,p6,p7,p6) - za(p2
     &    ,p4)*za(p3,p5)*zb(p1,p6)*zb(p4,p7)*iza(p1,p7)*iza(p6,p7)*izb(
     &    p1,p7)*izb(p5,p7)*izb(p6,p7)*zab2(p7,p1,p6,p4)*zab2(p7,p1,p7,
     &    p6) + za(p2,p4)*za(p3,p5)*zb(p1,p6)*zb(p4,p7)*iza(p6,p7)**2*
     &    izb(p5,p7)*izb(p6,p7)**2*zab2(p7,p1,p6,p4)*zab2(p7,p6,p7,p6)
     &     - 2*za(p2,p5)*za(p3,p5)*zb(p1,p6)*zb(p4,p5)*iza(p1,p7)*iza(
     &    p6,p7)*izb(p1,p7)*izb(p5,p7)*izb(p6,p7)*zab2(p7,p1,p6,p7)*
     &    zab2(p7,p1,p7,p6) + 2*za(p2,p5)*za(p3,p5)*zb(p1,p6)*zb(p4,p5)
     &    *iza(p6,p7)**2*izb(p5,p7)*izb(p6,p7)**2*zab2(p7,p1,p6,p7)*
     &    zab2(p7,p6,p7,p6) + za(p2,p5)*za(p3,p5)*zb(p1,p6)*zb(p4,p7)*
     &    iza(p1,p7)*iza(p6,p7)*izb(p1,p7)*izb(p5,p7)*izb(p6,p7)*zab2(
     &    p7,p1,p6,p5)*zab2(p7,p1,p7,p6) - za(p2,p5)*za(p3,p5)*zb(p1,p6
     &    )*zb(p4,p7)*iza(p6,p7)**2*izb(p5,p7)*izb(p6,p7)**2*zab2(p7,p1
     &    ,p6,p5)*zab2(p7,p6,p7,p6) )
      a(1,2,1) = a(1,2,1) + prop34**(-1)*prop345**(-1)*s267**(-1) * (
     &    za(p2,p7)**2*za(p3,p5)*zb(p1,p3)*zb(p4,p7)*iza(p2,p6)*iza(p6,
     &    p7)*izb(p5,p7)*izb(p6,p7)*zab2(p3,p2,p7,p6) - za(p2,p7)**2*
     &    za(p3,p5)*zb(p1,p4)*zb(p3,p7)*iza(p2,p6)*iza(p6,p7)*izb(p5,p7
     &    )*izb(p6,p7)*zab2(p3,p2,p7,p6) + za(p2,p7)**2*za(p3,p5)*zb(p1
     &    ,p4)*zb(p4,p7)*iza(p2,p6)*iza(p6,p7)*izb(p5,p7)*izb(p6,p7)*
     &    zab2(p4,p2,p7,p6) - za(p2,p7)**2*za(p3,p5)*zb(p1,p5)*zb(p4,p7
     &    )*iza(p2,p6)*iza(p6,p7)*izb(p5,p7)*izb(p6,p7)*zab2(p5,p2,p7,
     &    p6) + 2*za(p2,p7)**2*za(p3,p5)*zb(p1,p7)*zb(p4,p5)*iza(p2,p6)
     &    *iza(p6,p7)*izb(p5,p7)*izb(p6,p7)*zab2(p5,p2,p7,p6) - za(p2,
     &    p7)**2*za(p4,p5)*zb(p1,p4)*zb(p4,p7)*iza(p2,p6)*iza(p6,p7)*
     &    izb(p5,p7)*izb(p6,p7)*zab2(p3,p2,p7,p6) + za(p2,p7)**2*zb(p1,
     &    p4)*iza(p2,p6)*iza(p6,p7)*izb(p5,p7)*izb(p6,p7)*zab2(p3,p2,p7
     &    ,p6)*zab2(p5,p3,p4,p7) - za(p2,p7)*za(p3,p5)*zb(p1,p3)*zb(p4,
     &    p7)*iza(p6,p7)**2*izb(p5,p7)*izb(p6,p7)**2*zab2(p3,p2,p7,p6)*
     &    zab2(p7,p6,p7,p6) )
      a(1,2,1) = a(1,2,1) + prop34**(-1)*prop345**(-1)*s267**(-1) * (
     &    za(p2,p7)*za(p3,p5)*zb(p1,p4)*zb(p3,p7)*iza(p6,p7)**2*izb(p5,
     &    p7)*izb(p6,p7)**2*zab2(p3,p2,p7,p6)*zab2(p7,p6,p7,p6) - za(p2
     &    ,p7)*za(p3,p5)*zb(p1,p4)*zb(p4,p7)*iza(p6,p7)**2*izb(p5,p7)*
     &    izb(p6,p7)**2*zab2(p4,p2,p7,p6)*zab2(p7,p6,p7,p6) + za(p2,p7)
     &    *za(p3,p5)*zb(p1,p5)*zb(p4,p7)*iza(p6,p7)**2*izb(p5,p7)*izb(
     &    p6,p7)**2*zab2(p5,p2,p7,p6)*zab2(p7,p6,p7,p6) - 2*za(p2,p7)*
     &    za(p3,p5)*zb(p1,p7)*zb(p4,p5)*iza(p6,p7)**2*izb(p5,p7)*izb(p6
     &    ,p7)**2*zab2(p5,p2,p7,p6)*zab2(p7,p6,p7,p6) + za(p2,p7)*za(p4
     &    ,p5)*zb(p1,p4)*zb(p4,p7)*iza(p6,p7)**2*izb(p5,p7)*izb(p6,p7)
     &    **2*zab2(p3,p2,p7,p6)*zab2(p7,p6,p7,p6) - za(p2,p7)*zb(p1,p4)
     &    *iza(p6,p7)**2*izb(p5,p7)*izb(p6,p7)**2*zab2(p3,p2,p7,p6)*
     &    zab2(p5,p3,p4,p7)*zab2(p7,p6,p7,p6) )
      a(1,2,1) = a(1,2,1) + prop34**(-1)*prop345**(-1) * (  - za(p2,p3)
     &    *za(p2,p7)*za(p3,p5)*zb(p1,p6)*zb(p3,p7)*iza(p1,p7)*iza(p2,p6
     &    )*iza(p6,p7)*izb(p1,p7)*izb(p5,p7)*izb(p6,p7)*zab2(p7,p1,p7,
     &    p4) - za(p2,p3)*za(p2,p7)*za(p4,p5)*zb(p1,p6)*zb(p4,p7)*iza(
     &    p1,p7)*iza(p2,p6)*iza(p6,p7)*izb(p1,p7)*izb(p5,p7)*izb(p6,p7)
     &    *zab2(p7,p1,p7,p4) + za(p2,p3)*za(p2,p7)*zb(p1,p6)*iza(p1,p7)
     &    *iza(p2,p6)*iza(p6,p7)*izb(p1,p7)*izb(p5,p7)*izb(p6,p7)*zab2(
     &    p5,p3,p4,p7)*zab2(p7,p1,p7,p4) - 2*za(p2,p5)*za(p2,p7)*za(p3,
     &    p5)*zb(p1,p6)*zb(p4,p5)*iza(p2,p6)*iza(p6,p7)*izb(p5,p7)*izb(
     &    p6,p7) - za(p2,p7)*za(p3,p2)*za(p3,p5)*zb(p1,p6)*zb(p2,p6)*
     &    zb(p4,p7)*iza(p1,p7)*iza(p2,p6)*iza(p6,p7)*izb(p1,p7)*izb(p2,
     &    p6)*izb(p5,p7)*izb(p6,p7)*zab2(p7,p1,p7,p3) - za(p2,p7)*za(p3
     &    ,p5)*za(p4,p2)*zb(p1,p6)*zb(p2,p6)*zb(p4,p7)*iza(p1,p7)*iza(
     &    p2,p6)*iza(p6,p7)*izb(p1,p7)*izb(p2,p6)*izb(p5,p7)*izb(p6,p7)
     &    *zab2(p7,p1,p7,p4) + za(p2,p7)*za(p3,p5)*za(p5,p2)*zb(p1,p6)*
     &    zb(p2,p6)*zb(p4,p7)*iza(p1,p7)*iza(p2,p6)*iza(p6,p7)*izb(p1,
     &    p7)*izb(p2,p6)*izb(p5,p7)*izb(p6,p7)*zab2(p7,p1,p7,p5) )
      a(1,2,1) = a(1,2,1) + prop34**(-1)*s167**(-1) * (  - za(p2,p3)*
     &    zb(p1,p6)*iza(p1,p7)*iza(p6,p7)*izb(p1,p7)*izb(p5,p7)*izb(p6,
     &    p7)*zab2(p5,p2,p3,p4)*zab2(p7,p1,p6,p7)*zab2(p7,p1,p7,p6)*
     &    s234**(-1) + za(p2,p3)*zb(p1,p6)*iza(p6,p7)**2*izb(p5,p7)*
     &    izb(p6,p7)**2*zab2(p5,p2,p3,p4)*zab2(p7,p1,p6,p7)*zab2(p7,p6,
     &    p7,p6)*s234**(-1) - zb(p1,p6)*iza(p1,p7)*iza(p6,p7)*izb(p1,p7
     &    )*izb(p2,p5)*izb(p5,p7)*izb(p6,p7)*zab2(p3,p2,p5,p7)*zab2(p7,
     &    p1,p6,p4)*zab2(p7,p1,p7,p6) + zb(p1,p6)*iza(p6,p7)**2*izb(p2,
     &    p5)*izb(p5,p7)*izb(p6,p7)**2*zab2(p3,p2,p5,p7)*zab2(p7,p1,p6,
     &    p4)*zab2(p7,p6,p7,p6) )
      a(1,2,1) = a(1,2,1) + prop34**(-1)*s167**(-1)*Qsum * ( za(p2,p3)*
     &    zb(p1,p6)*iza(p1,p7)*iza(p6,p7)*izb(p1,p7)*izb(p5,p7)*izb(p6,
     &    p7)*zab2(p5,p2,p3,p4)*zab2(p7,p1,p6,p7)*zab2(p7,p1,p7,p6)*
     &    s234**(-1) - za(p2,p3)*zb(p1,p6)*iza(p6,p7)**2*izb(p5,p7)*
     &    izb(p6,p7)**2*zab2(p5,p2,p3,p4)*zab2(p7,p1,p6,p7)*zab2(p7,p6,
     &    p7,p6)*s234**(-1) - zb(p1,p6)*iza(p1,p7)*iza(p6,p7)*izb(p1,p7
     &    )*izb(p2,p5)*izb(p5,p7)*izb(p6,p7)*zab2(p3,p2,p5,p7)*zab2(p7,
     &    p1,p6,p4)*zab2(p7,p1,p7,p6) + zb(p1,p6)*iza(p6,p7)**2*izb(p2,
     &    p5)*izb(p5,p7)*izb(p6,p7)**2*zab2(p3,p2,p5,p7)*zab2(p7,p1,p6,
     &    p4)*zab2(p7,p6,p7,p6) )
      a(1,2,1) = a(1,2,1) + prop34**(-1)*s267**(-1) * ( za(p2,p7)**2*
     &    zb(p1,p4)*iza(p2,p6)*iza(p6,p7)*izb(p5,p7)*izb(p6,p7)*zab2(p3
     &    ,p1,p4,p7)*zab2(p5,p2,p7,p6)*s134**(-1) + za(p2,p7)**2*zb(p1,
     &    p7)*iza(p1,p5)*iza(p2,p6)*iza(p6,p7)*izb(p1,p5)*izb(p5,p7)*
     &    izb(p6,p7)*zab2(p3,p2,p7,p6)*zab2(p5,p1,p5,p4) - za(p2,p7)*
     &    zb(p1,p4)*iza(p6,p7)**2*izb(p5,p7)*izb(p6,p7)**2*zab2(p3,p1,
     &    p4,p7)*zab2(p5,p2,p7,p6)*zab2(p7,p6,p7,p6)*s134**(-1) - za(p2
     &    ,p7)*zb(p1,p7)*iza(p1,p5)*iza(p6,p7)**2*izb(p1,p5)*izb(p5,p7)
     &    *izb(p6,p7)**2*zab2(p3,p2,p7,p6)*zab2(p5,p1,p5,p4)*zab2(p7,p6
     &    ,p7,p6) )
      a(1,2,1) = a(1,2,1) + prop34**(-1)*s267**(-1)*Qsum * ( za(p2,p7)
     &    **2*zb(p1,p4)*iza(p2,p6)*iza(p6,p7)*izb(p5,p7)*izb(p6,p7)*
     &    zab2(p3,p1,p4,p7)*zab2(p5,p2,p7,p6)*s134**(-1) - za(p2,p7)**2
     &    *zb(p1,p7)*iza(p1,p5)*iza(p2,p6)*iza(p6,p7)*izb(p1,p5)*izb(p5
     &    ,p7)*izb(p6,p7)*zab2(p3,p2,p7,p6)*zab2(p5,p1,p5,p4) - za(p2,
     &    p7)*zb(p1,p4)*iza(p6,p7)**2*izb(p5,p7)*izb(p6,p7)**2*zab2(p3,
     &    p1,p4,p7)*zab2(p5,p2,p7,p6)*zab2(p7,p6,p7,p6)*s134**(-1) +
     &    za(p2,p7)*zb(p1,p7)*iza(p1,p5)*iza(p6,p7)**2*izb(p1,p5)*izb(
     &    p5,p7)*izb(p6,p7)**2*zab2(p3,p2,p7,p6)*zab2(p5,p1,p5,p4)*
     &    zab2(p7,p6,p7,p6) )
      a(1,2,1) = a(1,2,1) + prop34**(-1) * ( za(p2,p3)*za(p2,p7)*zb(p1,
     &    p6)*iza(p2,p6)*iza(p6,p7)*izb(p5,p7)*izb(p6,p7)*zab2(p5,p1,p7
     &    ,p4)*s157**(-1) - za(p2,p3)*za(p2,p7)*zb(p1,p7)*iza(p1,p5)*
     &    iza(p2,p6)*iza(p6,p7)*izb(p1,p5)*izb(p5,p7)*izb(p6,p7)*zab2(
     &    p5,p1,p5,p6)*zab2(p7,p1,p5,p4)*s157**(-1) + za(p2,p3)*zb(p1,
     &    p6)*iza(p6,p7)*izb(p5,p7)*izb(p6,p7)*zab2(p5,p1,p7,p6)*zab2(
     &    p7,p2,p3,p4)*s234**(-1)*s157**(-1) - za(p2,p3)*zb(p1,p7)*iza(
     &    p1,p5)*iza(p6,p7)**2*izb(p1,p5)*izb(p5,p7)*izb(p6,p7)**2*
     &    zab2(p5,p1,p5,p6)*zab2(p7,p2,p3,p4)*zab2(p7,p6,p7,p6)*
     &    s234**(-1) - za(p2,p3)*zb(p1,p7)*iza(p1,p5)*iza(p6,p7)*izb(p1
     &    ,p5)*izb(p5,p7)*izb(p6,p7)*zab2(p5,p1,p5,p6)*zab2(p7,p1,p5,p6
     &    )*zab2(p7,p2,p3,p4)*s234**(-1)*s157**(-1) + za(p2,p5)*za(p2,
     &    p7)*zb(p1,p4)*iza(p2,p6)*iza(p6,p7)*izb(p5,p7)*izb(p6,p7)*
     &    zab2(p3,p1,p4,p6)*zab3(p7,p2,p5,p6,p7)*s134**(-1)*s256**(-1)
     &     - za(p2,p5)*za(p2,p7)*zb(p1,p6)*iza(p1,p7)*iza(p2,p6)*iza(p6
     &    ,p7)*izb(p1,p7)*izb(p5,p7)*izb(p6,p7)*zab2(p7,p1,p7,p4)*zab3(
     &    p3,p2,p5,p6,p7)*s256**(-1) )
      a(1,2,1) = a(1,2,1) + prop34**(-1) * ( zb(p1,p4)*iza(p6,p7)**2*
     &    izb(p2,p5)*izb(p5,p7)*izb(p6,p7)**2*zab2(p3,p1,p4,p6)*zab2(p7
     &    ,p2,p5,p7)*zab2(p7,p6,p7,p6)*s134**(-1) - zb(p1,p4)*iza(p6,p7
     &    )*izb(p2,p5)*izb(p5,p7)*izb(p6,p7)*zab2(p3,p1,p4,p6)*zab2(p7,
     &    p2,p5,p6)*zab2(p7,p2,p5,p7)*s134**(-1)*s256**(-1) + zb(p1,p6)
     &    *iza(p1,p7)*iza(p6,p7)*izb(p1,p7)*izb(p2,p5)*izb(p5,p7)*izb(
     &    p6,p7)*zab2(p3,p2,p5,p6)*zab2(p7,p1,p7,p4)*zab2(p7,p2,p5,p7)*
     &    s256**(-1) )
      a(1,2,1) = a(1,2,1) + prop34**(-1)*Qsum * (  - za(p2,p3)*za(p2,p7
     &    )*zb(p1,p6)*iza(p2,p6)*iza(p6,p7)*izb(p5,p7)*izb(p6,p7)*zab2(
     &    p5,p1,p7,p4)*s157**(-1) + za(p2,p3)*za(p2,p7)*zb(p1,p7)*iza(
     &    p1,p5)*iza(p2,p6)*iza(p6,p7)*izb(p1,p5)*izb(p5,p7)*izb(p6,p7)
     &    *zab2(p5,p1,p5,p6)*zab2(p7,p1,p5,p4)*s157**(-1) - za(p2,p3)*
     &    zb(p1,p6)*iza(p6,p7)*izb(p5,p7)*izb(p6,p7)*zab2(p5,p1,p7,p6)*
     &    zab2(p7,p2,p3,p4)*s234**(-1)*s157**(-1) + za(p2,p3)*zb(p1,p7)
     &    *iza(p1,p5)*iza(p6,p7)**2*izb(p1,p5)*izb(p5,p7)*izb(p6,p7)**2
     &    *zab2(p5,p1,p5,p6)*zab2(p7,p2,p3,p4)*zab2(p7,p6,p7,p6)*
     &    s234**(-1) + za(p2,p3)*zb(p1,p7)*iza(p1,p5)*iza(p6,p7)*izb(p1
     &    ,p5)*izb(p5,p7)*izb(p6,p7)*zab2(p5,p1,p5,p6)*zab2(p7,p1,p5,p6
     &    )*zab2(p7,p2,p3,p4)*s234**(-1)*s157**(-1) + za(p2,p5)*za(p2,
     &    p7)*zb(p1,p4)*iza(p2,p6)*iza(p6,p7)*izb(p5,p7)*izb(p6,p7)*
     &    zab2(p3,p1,p4,p6)*zab3(p7,p2,p5,p6,p7)*s134**(-1)*s256**(-1)
     &     - za(p2,p5)*za(p2,p7)*zb(p1,p6)*iza(p1,p7)*iza(p2,p6)*iza(p6
     &    ,p7)*izb(p1,p7)*izb(p5,p7)*izb(p6,p7)*zab2(p7,p1,p7,p4)*zab3(
     &    p3,p2,p5,p6,p7)*s256**(-1) )
      a(1,2,1) = a(1,2,1) + prop34**(-1)*Qsum * ( zb(p1,p4)*iza(p6,p7)
     &    **2*izb(p2,p5)*izb(p5,p7)*izb(p6,p7)**2*zab2(p3,p1,p4,p6)*
     &    zab2(p7,p2,p5,p7)*zab2(p7,p6,p7,p6)*s134**(-1) - zb(p1,p4)*
     &    iza(p6,p7)*izb(p2,p5)*izb(p5,p7)*izb(p6,p7)*zab2(p3,p1,p4,p6)
     &    *zab2(p7,p2,p5,p6)*zab2(p7,p2,p5,p7)*s134**(-1)*s256**(-1) +
     &    zb(p1,p6)*iza(p1,p7)*iza(p6,p7)*izb(p1,p7)*izb(p2,p5)*izb(p5,
     &    p7)*izb(p6,p7)*zab2(p3,p2,p5,p6)*zab2(p7,p1,p7,p4)*zab2(p7,p2
     &    ,p5,p7)*s256**(-1) )
      a(1,2,1) = a(1,2,1) + prop345**(-1)*s167**(-1) * (  - 2*za(p2,p3)
     &    *zb(p1,p6)*zb(p3,p7)*iza(p1,p7)*iza(p6,p7)*izb(p1,p7)*izb(p3,
     &    p5)*izb(p5,p7)*izb(p6,p7)*zab2(p7,p1,p6,p4)*zab2(p7,p1,p7,p6)
     &     + 2*za(p2,p3)*zb(p1,p6)*zb(p3,p7)*iza(p6,p7)**2*izb(p3,p5)*
     &    izb(p5,p7)*izb(p6,p7)**2*zab2(p7,p1,p6,p4)*zab2(p7,p6,p7,p6)
     &     - 2*za(p2,p5)*zb(p1,p6)*iza(p1,p7)*iza(p6,p7)*izb(p1,p7)*
     &    izb(p3,p5)*izb(p6,p7)*zab2(p7,p1,p6,p4)*zab2(p7,p1,p7,p6) + 2
     &    *za(p2,p5)*zb(p1,p6)*iza(p6,p7)**2*izb(p3,p5)*izb(p6,p7)**2*
     &    zab2(p7,p1,p6,p4)*zab2(p7,p6,p7,p6) )
      a(1,2,1) = a(1,2,1) + prop345**(-1)*s267**(-1) * ( 2*za(p2,p7)**2
     &    *zb(p1,p4)*zb(p3,p7)*iza(p2,p6)*iza(p6,p7)*izb(p3,p5)*izb(p5,
     &    p7)*izb(p6,p7)*zab2(p3,p2,p7,p6) + 2*za(p2,p7)**2*zb(p1,p4)*
     &    iza(p2,p6)*iza(p6,p7)*izb(p3,p5)*izb(p6,p7)*zab2(p5,p2,p7,p6)
     &     - 2*za(p2,p7)*zb(p1,p4)*zb(p3,p7)*iza(p6,p7)**2*izb(p3,p5)*
     &    izb(p5,p7)*izb(p6,p7)**2*zab2(p3,p2,p7,p6)*zab2(p7,p6,p7,p6)
     &     - 2*za(p2,p7)*zb(p1,p4)*iza(p6,p7)**2*izb(p3,p5)*izb(p6,p7)
     &    **2*zab2(p5,p2,p7,p6)*zab2(p7,p6,p7,p6) )
      a(1,2,1) = a(1,2,1) + prop345**(-1) * (  - 2*za(p2,p7)*za(p3,p2)*
     &    zb(p1,p6)*zb(p2,p6)*zb(p3,p7)*iza(p1,p7)*iza(p2,p6)*iza(p6,p7
     &    )*izb(p1,p7)*izb(p2,p6)*izb(p3,p5)*izb(p5,p7)*izb(p6,p7)*
     &    zab2(p7,p1,p7,p4) - 2*za(p2,p7)*za(p5,p2)*zb(p1,p6)*zb(p2,p6)
     &    *iza(p1,p7)*iza(p2,p6)*iza(p6,p7)*izb(p1,p7)*izb(p2,p6)*izb(
     &    p3,p5)*izb(p6,p7)*zab2(p7,p1,p7,p4) )
      a(1,2,2)= + prop34**(-1)*prop345**(-1)*s167**(-1) * (  - za(p2,p3
     &    )*za(p3,p5)*zb(p1,p3)*zb(p1,p6)*iza(p1,p7)*iza(p6,p7)*izb(p1,
     &    p5)*zab2(p1,p6,p7,p4) - za(p2,p3)*za(p3,p5)*zb(p1,p3)*zb(p1,
     &    p7)*iza(p1,p6)*iza(p6,p7)*izb(p1,p5)*zab2(p1,p6,p7,p4) - za(
     &    p2,p3)*za(p3,p5)*zb(p1,p3)*zb(p6,p7)*iza(p1,p6)*iza(p1,p7)*
     &    izb(p1,p5)*zab2(p1,p6,p7,p4) + za(p2,p3)*za(p3,p5)*zb(p1,p4)*
     &    zb(p1,p6)*iza(p1,p7)*iza(p6,p7)*izb(p1,p5)*zab2(p1,p6,p7,p3)
     &     + za(p2,p3)*za(p3,p5)*zb(p1,p4)*zb(p1,p7)*iza(p1,p6)*iza(p6,
     &    p7)*izb(p1,p5)*zab2(p1,p6,p7,p3) + za(p2,p3)*za(p3,p5)*zb(p1,
     &    p4)*zb(p6,p7)*iza(p1,p6)*iza(p1,p7)*izb(p1,p5)*zab2(p1,p6,p7,
     &    p3) - za(p2,p3)*za(p4,p5)*zb(p1,p4)*zb(p1,p6)*iza(p1,p7)*iza(
     &    p6,p7)*izb(p1,p5)*zab2(p1,p6,p7,p4) - za(p2,p3)*za(p4,p5)*zb(
     &    p1,p4)*zb(p1,p7)*iza(p1,p6)*iza(p6,p7)*izb(p1,p5)*zab2(p1,p6,
     &    p7,p4) - za(p2,p3)*za(p4,p5)*zb(p1,p4)*zb(p6,p7)*iza(p1,p6)*
     &    iza(p1,p7)*izb(p1,p5)*zab2(p1,p6,p7,p4) - za(p2,p3)*zb(p1,p6)
     &    *iza(p1,p7)*iza(p6,p7)*izb(p1,p5)*zab2(p1,p6,p7,p4)*zab2(p5,
     &    p3,p4,p1) )
      a(1,2,2) = a(1,2,2) + prop34**(-1)*prop345**(-1)*s167**(-1) * (
     &     - za(p2,p3)*zb(p1,p7)*iza(p1,p6)*iza(p6,p7)*izb(p1,p5)*zab2(
     &    p1,p6,p7,p4)*zab2(p5,p3,p4,p1) - za(p2,p3)*zb(p6,p7)*iza(p1,
     &    p6)*iza(p1,p7)*izb(p1,p5)*zab2(p1,p6,p7,p4)*zab2(p5,p3,p4,p1)
     &     + za(p2,p4)*za(p3,p5)*zb(p1,p4)*zb(p1,p6)*iza(p1,p7)*iza(p6,
     &    p7)*izb(p1,p5)*zab2(p1,p6,p7,p4) + za(p2,p4)*za(p3,p5)*zb(p1,
     &    p4)*zb(p1,p7)*iza(p1,p6)*iza(p6,p7)*izb(p1,p5)*zab2(p1,p6,p7,
     &    p4) + za(p2,p4)*za(p3,p5)*zb(p1,p4)*zb(p6,p7)*iza(p1,p6)*iza(
     &    p1,p7)*izb(p1,p5)*zab2(p1,p6,p7,p4) - za(p2,p5)*za(p3,p5)*zb(
     &    p1,p4)*zb(p1,p6)*iza(p1,p7)*iza(p6,p7)*izb(p1,p5)*zab2(p1,p6,
     &    p7,p5) - za(p2,p5)*za(p3,p5)*zb(p1,p4)*zb(p1,p7)*iza(p1,p6)*
     &    iza(p6,p7)*izb(p1,p5)*zab2(p1,p6,p7,p5) - za(p2,p5)*za(p3,p5)
     &    *zb(p1,p4)*zb(p6,p7)*iza(p1,p6)*iza(p1,p7)*izb(p1,p5)*zab2(p1
     &    ,p6,p7,p5) - 2*za(p2,p5)*za(p3,p5)*zb(p1,p6)*zb(p4,p5)*iza(p1
     &    ,p7)*iza(p6,p7)*izb(p1,p5)*zab2(p1,p6,p7,p1) - 2*za(p2,p5)*
     &    za(p3,p5)*zb(p1,p7)*zb(p4,p5)*iza(p1,p6)*iza(p6,p7)*izb(p1,p5
     &    )*zab2(p1,p6,p7,p1) )
      a(1,2,2) = a(1,2,2) + prop34**(-1)*prop345**(-1)*s167**(-1) * (
     &     - 2*za(p2,p5)*za(p3,p5)*zb(p4,p5)*zb(p6,p7)*iza(p1,p6)*iza(
     &    p1,p7)*izb(p1,p5)*zab2(p1,p6,p7,p1) )
      a(1,2,2) = a(1,2,2) + prop34**(-1)*prop345**(-1)*s267**(-1) * (
     &    za(p1,p2)**2*za(p3,p5)*zb(p1,p4)**2*iza(p1,p6)*iza(p1,p7)*
     &    iza(p2,p6)*izb(p1,p5)*zab2(p4,p2,p6,p7) - za(p1,p2)**2*za(p3,
     &    p5)*zb(p1,p4)*iza(p1,p6)*iza(p1,p7)*iza(p2,p6)*zab2(p5,p2,p6,
     &    p7) - za(p1,p2)**2*za(p4,p5)*zb(p1,p4)**2*iza(p1,p6)*iza(p1,
     &    p7)*iza(p2,p6)*izb(p1,p5)*zab2(p3,p2,p6,p7) - za(p1,p2)**2*
     &    zb(p1,p4)*iza(p1,p6)*iza(p1,p7)*iza(p2,p6)*izb(p1,p5)*zab2(p3
     &    ,p2,p6,p7)*zab2(p5,p3,p4,p1) + za(p1,p2)*za(p3,p5)*zb(p1,p4)
     &    **2*iza(p1,p6)*iza(p6,p7)*izb(p1,p5)*zab2(p4,p2,p6,p7) + za(
     &    p1,p2)*za(p3,p5)*zb(p1,p4)**2*iza(p1,p7)*iza(p6,p7)*izb(p1,p5
     &    )*zab2(p4,p2,p7,p6) - za(p1,p2)*za(p3,p5)*zb(p1,p4)*iza(p1,p6
     &    )*iza(p6,p7)*zab2(p5,p2,p6,p7) - za(p1,p2)*za(p3,p5)*zb(p1,p4
     &    )*iza(p1,p7)*iza(p6,p7)*zab2(p5,p2,p7,p6) - za(p1,p2)*za(p4,
     &    p5)*zb(p1,p4)**2*iza(p1,p6)*iza(p6,p7)*izb(p1,p5)*zab2(p3,p2,
     &    p6,p7) - za(p1,p2)*za(p4,p5)*zb(p1,p4)**2*iza(p1,p7)*iza(p6,
     &    p7)*izb(p1,p5)*zab2(p3,p2,p7,p6) )
      a(1,2,2) = a(1,2,2) + prop34**(-1)*prop345**(-1)*s267**(-1) * (
     &     - za(p1,p2)*zb(p1,p4)*iza(p1,p6)*iza(p6,p7)*izb(p1,p5)*zab2(
     &    p3,p2,p6,p7)*zab2(p5,p3,p4,p1) - za(p1,p2)*zb(p1,p4)*iza(p1,
     &    p7)*iza(p6,p7)*izb(p1,p5)*zab2(p3,p2,p7,p6)*zab2(p5,p3,p4,p1)
     &     )
      a(1,2,2) = a(1,2,2) + prop34**(-1)*prop345**(-1) * ( za(p1,p2)*
     &    za(p1,p7)*za(p3,p2)*za(p3,p5)*zb(p1,p4)*zb(p2,p6)*zb(p7,p3)*
     &    iza(p1,p6)*iza(p1,p7)**2*iza(p2,p6)*izb(p1,p5)*izb(p2,p6) +
     &    za(p1,p2)*za(p1,p7)*za(p3,p5)*za(p4,p2)*zb(p1,p4)*zb(p2,p6)*
     &    zb(p7,p4)*iza(p1,p6)*iza(p1,p7)**2*iza(p2,p6)*izb(p1,p5)*izb(
     &    p2,p6) - za(p1,p2)*za(p1,p7)*za(p3,p5)*za(p5,p2)*zb(p1,p4)*
     &    zb(p2,p6)*zb(p7,p5)*iza(p1,p6)*iza(p1,p7)**2*iza(p2,p6)*izb(
     &    p1,p5)*izb(p2,p6) - za(p1,p2)*za(p2,p3)*za(p3,p5)*zb(p1,p3)*
     &    zb(p4,p7)*iza(p1,p6)*iza(p1,p7)*iza(p2,p6)*izb(p1,p5) - za(p1
     &    ,p2)*za(p2,p3)*za(p4,p5)*zb(p1,p4)*zb(p4,p7)*iza(p1,p6)*iza(
     &    p1,p7)*iza(p2,p6)*izb(p1,p5) - za(p1,p2)*za(p2,p3)*zb(p4,p7)*
     &    iza(p1,p6)*iza(p1,p7)*iza(p2,p6)*izb(p1,p5)*zab2(p5,p3,p4,p1)
     &     - 2*za(p1,p2)*za(p2,p5)*za(p3,p5)*zb(p1,p7)*zb(p4,p5)*iza(p1
     &    ,p6)*iza(p1,p7)*iza(p2,p6)*izb(p1,p5) )
      a(1,2,2) = a(1,2,2) + prop34**(-1)*s167**(-1) * (  - za(p2,p3)*
     &    zb(p1,p6)*iza(p1,p7)*iza(p6,p7)*izb(p1,p5)*zab2(p1,p6,p7,p1)*
     &    zab2(p5,p2,p3,p4)*s234**(-1) - za(p2,p3)*zb(p1,p7)*iza(p1,p6)
     &    *iza(p6,p7)*izb(p1,p5)*zab2(p1,p6,p7,p1)*zab2(p5,p2,p3,p4)*
     &    s234**(-1) - za(p2,p3)*zb(p6,p7)*iza(p1,p6)*iza(p1,p7)*izb(p1
     &    ,p5)*zab2(p1,p6,p7,p1)*zab2(p5,p2,p3,p4)*s234**(-1) - zb(p1,
     &    p6)*iza(p1,p7)*iza(p6,p7)*izb(p1,p5)*izb(p2,p5)*zab2(p1,p6,p7
     &    ,p4)*zab2(p3,p2,p5,p1) - zb(p1,p7)*iza(p1,p6)*iza(p6,p7)*izb(
     &    p1,p5)*izb(p2,p5)*zab2(p1,p6,p7,p4)*zab2(p3,p2,p5,p1) - zb(p6
     &    ,p7)*iza(p1,p6)*iza(p1,p7)*izb(p1,p5)*izb(p2,p5)*zab2(p1,p6,
     &    p7,p4)*zab2(p3,p2,p5,p1) )
      a(1,2,2) = a(1,2,2) + prop34**(-1)*s167**(-1)*Qsum * ( za(p2,p3)*
     &    zb(p1,p6)*iza(p1,p7)*iza(p6,p7)*izb(p1,p5)*zab2(p1,p6,p7,p1)*
     &    zab2(p5,p2,p3,p4)*s234**(-1) + za(p2,p3)*zb(p1,p7)*iza(p1,p6)
     &    *iza(p6,p7)*izb(p1,p5)*zab2(p1,p6,p7,p1)*zab2(p5,p2,p3,p4)*
     &    s234**(-1) + za(p2,p3)*zb(p6,p7)*iza(p1,p6)*iza(p1,p7)*izb(p1
     &    ,p5)*zab2(p1,p6,p7,p1)*zab2(p5,p2,p3,p4)*s234**(-1) - zb(p1,
     &    p6)*iza(p1,p7)*iza(p6,p7)*izb(p1,p5)*izb(p2,p5)*zab2(p1,p6,p7
     &    ,p4)*zab2(p3,p2,p5,p1) - zb(p1,p7)*iza(p1,p6)*iza(p6,p7)*izb(
     &    p1,p5)*izb(p2,p5)*zab2(p1,p6,p7,p4)*zab2(p3,p2,p5,p1) - zb(p6
     &    ,p7)*iza(p1,p6)*iza(p1,p7)*izb(p1,p5)*izb(p2,p5)*zab2(p1,p6,
     &    p7,p4)*zab2(p3,p2,p5,p1) )
      a(1,2,2) = a(1,2,2) + prop34**(-1)*s267**(-1) * (  - za(p1,p2)**2
     &    *zb(p1,p4)*iza(p1,p6)*iza(p1,p7)*iza(p2,p6)*izb(p1,p5)*zab2(
     &    p3,p1,p4,p1)*zab2(p5,p2,p6,p7)*s134**(-1) - za(p1,p2)*zb(p1,
     &    p4)*iza(p1,p6)*iza(p6,p7)*izb(p1,p5)*zab2(p3,p1,p4,p1)*zab2(
     &    p5,p2,p6,p7)*s134**(-1) - za(p1,p2)*zb(p1,p4)*iza(p1,p7)*iza(
     &    p6,p7)*izb(p1,p5)*zab2(p3,p1,p4,p1)*zab2(p5,p2,p7,p6)*
     &    s134**(-1) )
      a(1,2,2) = a(1,2,2) + prop34**(-1)*s267**(-1)*Qsum * (  - za(p1,
     &    p2)**2*zb(p1,p4)*iza(p1,p6)*iza(p1,p7)*iza(p2,p6)*izb(p1,p5)*
     &    zab2(p3,p1,p4,p1)*zab2(p5,p2,p6,p7)*s134**(-1) - za(p1,p2)*
     &    zb(p1,p4)*iza(p1,p6)*iza(p6,p7)*izb(p1,p5)*zab2(p3,p1,p4,p1)*
     &    zab2(p5,p2,p6,p7)*s134**(-1) - za(p1,p2)*zb(p1,p4)*iza(p1,p7)
     &    *iza(p6,p7)*izb(p1,p5)*zab2(p3,p1,p4,p1)*zab2(p5,p2,p7,p6)*
     &    s134**(-1) )
      a(1,2,2) = a(1,2,2) + prop34**(-1) * ( za(p1,p2)*za(p2,p3)*zb(p1,
     &    p7)*iza(p1,p6)*iza(p1,p7)*iza(p2,p6)*izb(p1,p5)*zab2(p5,p1,p7
     &    ,p4)*s157**(-1) + za(p1,p2)*za(p2,p5)*zb(p1,p4)*iza(p1,p6)*
     &    iza(p1,p7)*iza(p2,p6)*izb(p1,p5)*zab2(p3,p1,p4,p7)*zab3(p1,p2
     &    ,p5,p6,p1)*s134**(-1)*s256**(-1) + za(p1,p2)*za(p2,p5)*zb(p4,
     &    p7)*iza(p1,p6)*iza(p1,p7)*iza(p2,p6)*izb(p1,p5)*zab3(p3,p2,p5
     &    ,p6,p1)*s256**(-1) - za(p2,p3)*zb(p1,p7)*iza(p1,p6)*iza(p1,p7
     &    )*izb(p1,p5)*zab2(p1,p2,p3,p4)*zab2(p5,p1,p7,p6)*s234**(-1)*
     &    s157**(-1) + zb(p1,p4)*iza(p1,p6)*iza(p1,p7)*izb(p1,p5)*izb(
     &    p2,p5)*zab2(p1,p2,p5,p1)*zab2(p1,p2,p5,p6)*zab2(p3,p1,p4,p7)*
     &    s134**(-1)*s256**(-1) - zb(p1,p4)*iza(p1,p6)*iza(p6,p7)*izb(
     &    p1,p5)*izb(p2,p5)*zab2(p1,p2,p5,p1)*zab2(p3,p1,p4,p7)*
     &    s134**(-1) - zb(p1,p4)*iza(p1,p7)*iza(p6,p7)*izb(p1,p5)*izb(
     &    p2,p5)*zab2(p1,p2,p5,p1)*zab2(p3,p1,p4,p6)*s134**(-1) + zb(p4
     &    ,p7)*iza(p1,p6)*iza(p1,p7)*izb(p1,p5)*izb(p2,p5)*zab2(p1,p2,
     &    p5,p1)*zab2(p3,p2,p5,p6)*s256**(-1) )
      a(1,2,2) = a(1,2,2) + prop34**(-1)*Qsum * (  - za(p1,p2)*za(p2,p3
     &    )*zb(p1,p7)*iza(p1,p6)*iza(p1,p7)*iza(p2,p6)*izb(p1,p5)*zab2(
     &    p5,p1,p7,p4)*s157**(-1) + za(p1,p2)*za(p2,p5)*zb(p1,p4)*iza(
     &    p1,p6)*iza(p1,p7)*iza(p2,p6)*izb(p1,p5)*zab2(p3,p1,p4,p7)*
     &    zab3(p1,p2,p5,p6,p1)*s134**(-1)*s256**(-1) + za(p1,p2)*za(p2,
     &    p5)*zb(p4,p7)*iza(p1,p6)*iza(p1,p7)*iza(p2,p6)*izb(p1,p5)*
     &    zab3(p3,p2,p5,p6,p1)*s256**(-1) + za(p2,p3)*zb(p1,p7)*iza(p1,
     &    p6)*iza(p1,p7)*izb(p1,p5)*zab2(p1,p2,p3,p4)*zab2(p5,p1,p7,p6)
     &    *s234**(-1)*s157**(-1) + zb(p1,p4)*iza(p1,p6)*iza(p1,p7)*izb(
     &    p1,p5)*izb(p2,p5)*zab2(p1,p2,p5,p1)*zab2(p1,p2,p5,p6)*zab2(p3
     &    ,p1,p4,p7)*s134**(-1)*s256**(-1) - zb(p1,p4)*iza(p1,p6)*iza(
     &    p6,p7)*izb(p1,p5)*izb(p2,p5)*zab2(p1,p2,p5,p1)*zab2(p3,p1,p4,
     &    p7)*s134**(-1) - zb(p1,p4)*iza(p1,p7)*iza(p6,p7)*izb(p1,p5)*
     &    izb(p2,p5)*zab2(p1,p2,p5,p1)*zab2(p3,p1,p4,p6)*s134**(-1) +
     &    zb(p4,p7)*iza(p1,p6)*iza(p1,p7)*izb(p1,p5)*izb(p2,p5)*zab2(p1
     &    ,p2,p5,p1)*zab2(p3,p2,p5,p6)*s256**(-1) )
      a(1,2,2) = a(1,2,2) + prop345**(-1)*s167**(-1) * ( 2*za(p2,p3)*
     &    zb(p1,p3)*zb(p1,p6)*iza(p1,p7)*iza(p6,p7)*izb(p1,p5)*izb(p3,
     &    p5)*zab2(p1,p6,p7,p4) + 2*za(p2,p3)*zb(p1,p3)*zb(p1,p7)*iza(
     &    p1,p6)*iza(p6,p7)*izb(p1,p5)*izb(p3,p5)*zab2(p1,p6,p7,p4) + 2
     &    *za(p2,p3)*zb(p1,p3)*zb(p6,p7)*iza(p1,p6)*iza(p1,p7)*izb(p1,
     &    p5)*izb(p3,p5)*zab2(p1,p6,p7,p4) + 2*za(p2,p5)*zb(p1,p6)*iza(
     &    p1,p7)*iza(p6,p7)*izb(p3,p5)*zab2(p1,p6,p7,p4) + 2*za(p2,p5)*
     &    zb(p1,p7)*iza(p1,p6)*iza(p6,p7)*izb(p3,p5)*zab2(p1,p6,p7,p4)
     &     + 2*za(p2,p5)*zb(p6,p7)*iza(p1,p6)*iza(p1,p7)*izb(p3,p5)*
     &    zab2(p1,p6,p7,p4) )
      a(1,2,2) = a(1,2,2) + prop345**(-1)*s267**(-1) * ( 2*za(p1,p2)**2
     &    *zb(p1,p3)*zb(p1,p4)*iza(p1,p6)*iza(p1,p7)*iza(p2,p6)*izb(p1,
     &    p5)*izb(p3,p5)*zab2(p3,p2,p6,p7) + 2*za(p1,p2)**2*zb(p1,p4)*
     &    iza(p1,p6)*iza(p1,p7)*iza(p2,p6)*izb(p3,p5)*zab2(p5,p2,p6,p7)
     &     + 2*za(p1,p2)*zb(p1,p3)*zb(p1,p4)*iza(p1,p6)*iza(p6,p7)*izb(
     &    p1,p5)*izb(p3,p5)*zab2(p3,p2,p6,p7) + 2*za(p1,p2)*zb(p1,p3)*
     &    zb(p1,p4)*iza(p1,p7)*iza(p6,p7)*izb(p1,p5)*izb(p3,p5)*zab2(p3
     &    ,p2,p7,p6) + 2*za(p1,p2)*zb(p1,p4)*iza(p1,p6)*iza(p6,p7)*izb(
     &    p3,p5)*zab2(p5,p2,p6,p7) + 2*za(p1,p2)*zb(p1,p4)*iza(p1,p7)*
     &    iza(p6,p7)*izb(p3,p5)*zab2(p5,p2,p7,p6) )
      a(1,2,2) = a(1,2,2) + prop345**(-1) * ( 2*za(p1,p2)*za(p1,p7)*za(
     &    p3,p2)*zb(p1,p3)*zb(p2,p6)*zb(p7,p4)*iza(p1,p6)*iza(p1,p7)**2
     &    *iza(p2,p6)*izb(p1,p5)*izb(p2,p6)*izb(p3,p5) + 2*za(p1,p2)*
     &    za(p1,p7)*za(p5,p2)*zb(p2,p6)*zb(p7,p4)*iza(p1,p6)*iza(p1,p7)
     &    **2*iza(p2,p6)*izb(p2,p6)*izb(p3,p5) )

      a(2,1,1)= + prop34**(-1)*prop345**(-1)*s167**(-1) * (  - za(p2,p3
     &    )**2*zb(p1,p2)*zb(p3,p5)*iza(p1,p7)*iza(p2,p5)*izb(p1,p7)*
     &    izb(p2,p6)*izb(p2,p7)*zab2(p6,p1,p7,p4)*zab2(p7,p1,p7,p2) +
     &    za(p2,p3)**2*zb(p1,p2)*zb(p3,p5)*iza(p2,p5)*iza(p6,p7)*izb(p2
     &    ,p6)*izb(p2,p7)*izb(p6,p7)*zab2(p6,p1,p7,p4)*zab2(p7,p6,p7,p2
     &    ) + za(p2,p3)**2*zb(p1,p2)*zb(p3,p5)*iza(p2,p5)*izb(p2,p6)*
     &    izb(p6,p7)*zab2(p7,p1,p6,p4) + za(p2,p3)**2*zb(p1,p2)*zb(p4,
     &    p5)*iza(p1,p7)*iza(p2,p5)*izb(p1,p7)*izb(p2,p6)*izb(p2,p7)*
     &    zab2(p6,p1,p7,p3)*zab2(p7,p1,p7,p2) - za(p2,p3)**2*zb(p1,p2)*
     &    zb(p4,p5)*iza(p2,p5)*iza(p6,p7)*izb(p2,p6)*izb(p2,p7)*izb(p6,
     &    p7)*zab2(p6,p1,p7,p3)*zab2(p7,p6,p7,p2) - za(p2,p3)**2*zb(p1,
     &    p2)*zb(p4,p5)*iza(p2,p5)*izb(p2,p6)*izb(p6,p7)*zab2(p7,p1,p6,
     &    p3) - za(p2,p3)*zb(p1,p2)*zb(p4,p5)*iza(p1,p7)*izb(p1,p7)*
     &    izb(p2,p6)*izb(p2,p7)*zab2(p6,p1,p7,p5)*zab2(p7,p1,p7,p2) +
     &    za(p2,p3)*zb(p1,p2)*zb(p4,p5)*iza(p6,p7)*izb(p2,p6)*izb(p2,p7
     &    )*izb(p6,p7)*zab2(p6,p1,p7,p5)*zab2(p7,p6,p7,p2) )
      a(2,1,1) = a(2,1,1) + prop34**(-1)*prop345**(-1)*s167**(-1) * (
     &    za(p2,p3)*zb(p1,p2)*zb(p4,p5)*izb(p2,p6)*izb(p6,p7)*zab2(p7,
     &    p1,p6,p5) - za(p2,p3)*zb(p1,p2)*iza(p1,p7)*iza(p2,p5)*izb(p1,
     &    p7)*izb(p2,p6)*izb(p2,p7)*zab2(p2,p3,p4,p5)*zab2(p6,p1,p7,p4)
     &    *zab2(p7,p1,p7,p2) + za(p2,p3)*zb(p1,p2)*iza(p2,p5)*iza(p6,p7
     &    )*izb(p2,p6)*izb(p2,p7)*izb(p6,p7)*zab2(p2,p3,p4,p5)*zab2(p6,
     &    p1,p7,p4)*zab2(p7,p6,p7,p2) + za(p2,p3)*zb(p1,p2)*iza(p2,p5)*
     &    izb(p2,p6)*izb(p6,p7)*zab2(p2,p3,p4,p5)*zab2(p7,p1,p6,p4) )
      a(2,1,1) = a(2,1,1) + prop34**(-1)*prop345**(-1)*s267**(-1) * (
     &    za(p2,p3)*za(p2,p6)*zb(p1,p3)*zb(p4,p5)*iza(p2,p5)*iza(p6,p7)
     &    *izb(p2,p6)*izb(p2,p7)*izb(p6,p7)*zab2(p3,p6,p7,p2)*zab2(p7,
     &    p6,p7,p2) - za(p2,p3)*za(p2,p6)*zb(p1,p4)*zb(p3,p5)*iza(p2,p5
     &    )*iza(p6,p7)*izb(p2,p6)*izb(p2,p7)*izb(p6,p7)*zab2(p3,p6,p7,
     &    p2)*zab2(p7,p6,p7,p2) + za(p2,p3)*za(p2,p6)*zb(p1,p4)*zb(p4,
     &    p5)*iza(p2,p5)*iza(p6,p7)*izb(p2,p6)*izb(p2,p7)*izb(p6,p7)*
     &    zab2(p4,p6,p7,p2)*zab2(p7,p6,p7,p2) - za(p2,p3)*za(p2,p6)*zb(
     &    p1,p5)*zb(p4,p5)*iza(p2,p5)*iza(p6,p7)*izb(p2,p6)*izb(p2,p7)*
     &    izb(p6,p7)*zab2(p5,p6,p7,p2)*zab2(p7,p6,p7,p2) + za(p2,p3)*
     &    za(p2,p7)*zb(p1,p3)*zb(p4,p5)*iza(p2,p5)*izb(p2,p6)*izb(p6,p7
     &    )*zab2(p3,p6,p7,p2) - za(p2,p3)*za(p2,p7)*zb(p1,p4)*zb(p3,p5)
     &    *iza(p2,p5)*izb(p2,p6)*izb(p6,p7)*zab2(p3,p6,p7,p2) + za(p2,
     &    p3)*za(p2,p7)*zb(p1,p4)*zb(p4,p5)*iza(p2,p5)*izb(p2,p6)*izb(
     &    p6,p7)*zab2(p4,p6,p7,p2) - za(p2,p3)*za(p2,p7)*zb(p1,p5)*zb(
     &    p4,p5)*iza(p2,p5)*izb(p2,p6)*izb(p6,p7)*zab2(p5,p6,p7,p2) )
      a(2,1,1) = a(2,1,1) + prop34**(-1)*prop345**(-1)*s267**(-1) * (
     &    za(p2,p3)*zb(p1,p3)*zb(p4,p5)*iza(p2,p5)*izb(p2,p6)**2*izb(p2
     &    ,p7)*zab2(p3,p6,p7,p2)*zab2(p7,p2,p6,p2) - za(p2,p3)*zb(p1,p4
     &    )*zb(p3,p5)*iza(p2,p5)*izb(p2,p6)**2*izb(p2,p7)*zab2(p3,p6,p7
     &    ,p2)*zab2(p7,p2,p6,p2) + za(p2,p3)*zb(p1,p4)*zb(p4,p5)*iza(p2
     &    ,p5)*izb(p2,p6)**2*izb(p2,p7)*zab2(p4,p6,p7,p2)*zab2(p7,p2,p6
     &    ,p2) - za(p2,p3)*zb(p1,p5)*zb(p4,p5)*iza(p2,p5)*izb(p2,p6)**2
     &    *izb(p2,p7)*zab2(p5,p6,p7,p2)*zab2(p7,p2,p6,p2) - za(p2,p4)*
     &    za(p2,p6)*zb(p1,p4)*zb(p4,p5)*iza(p2,p5)*iza(p6,p7)*izb(p2,p6
     &    )*izb(p2,p7)*izb(p6,p7)*zab2(p3,p6,p7,p2)*zab2(p7,p6,p7,p2)
     &     - za(p2,p4)*za(p2,p7)*zb(p1,p4)*zb(p4,p5)*iza(p2,p5)*izb(p2,
     &    p6)*izb(p6,p7)*zab2(p3,p6,p7,p2) - za(p2,p4)*zb(p1,p4)*zb(p4,
     &    p5)*iza(p2,p5)*izb(p2,p6)**2*izb(p2,p7)*zab2(p3,p6,p7,p2)*
     &    zab2(p7,p2,p6,p2) - 2*za(p2,p6)*za(p3,p5)*zb(p1,p5)*zb(p4,p5)
     &    *iza(p2,p5)*iza(p6,p7)*izb(p2,p6)*izb(p2,p7)*izb(p6,p7)*zab2(
     &    p2,p6,p7,p2)*zab2(p7,p6,p7,p2) )
      a(2,1,1) = a(2,1,1) + prop34**(-1)*prop345**(-1)*s267**(-1) * (
     &     - za(p2,p6)*zb(p1,p4)*iza(p2,p5)*iza(p6,p7)*izb(p2,p6)*izb(
     &    p2,p7)*izb(p6,p7)*zab2(p2,p3,p4,p5)*zab2(p3,p6,p7,p2)*zab2(p7
     &    ,p6,p7,p2) - 2*za(p2,p7)*za(p3,p5)*zb(p1,p5)*zb(p4,p5)*iza(p2
     &    ,p5)*izb(p2,p6)*izb(p6,p7)*zab2(p2,p6,p7,p2) - za(p2,p7)*zb(
     &    p1,p4)*iza(p2,p5)*izb(p2,p6)*izb(p6,p7)*zab2(p2,p3,p4,p5)*
     &    zab2(p3,p6,p7,p2) - 2*za(p3,p5)*zb(p1,p5)*zb(p4,p5)*iza(p2,p5
     &    )*izb(p2,p6)**2*izb(p2,p7)*zab2(p2,p6,p7,p2)*zab2(p7,p2,p6,p2
     &    ) - zb(p1,p4)*iza(p2,p5)*izb(p2,p6)**2*izb(p2,p7)*zab2(p2,p3,
     &    p4,p5)*zab2(p3,p6,p7,p2)*zab2(p7,p2,p6,p2) )
      a(2,1,1) = a(2,1,1) + prop34**(-1)*prop345**(-1) * (  - za(p2,p3)
     &    *zb(p1,p2)*zb(p3,p5)*iza(p1,p7)*iza(p2,p5)*izb(p1,p7)*izb(p2,
     &    p6)**2*izb(p2,p7)*zab2(p3,p2,p6,p2)*zab2(p7,p1,p7,p4) + za(p2
     &    ,p3)*zb(p1,p2)*zb(p4,p5)*iza(p1,p7)*iza(p2,p5)*izb(p1,p7)*
     &    izb(p2,p6)**2*izb(p2,p7)*zab2(p3,p2,p6,p2)*zab2(p7,p1,p7,p3)
     &     + za(p2,p3)*zb(p1,p2)*zb(p4,p5)*iza(p1,p7)*iza(p2,p5)*izb(p1
     &    ,p7)*izb(p2,p6)**2*izb(p2,p7)*zab2(p4,p2,p6,p2)*zab2(p7,p1,p7
     &    ,p4) - za(p2,p3)*zb(p1,p2)*zb(p4,p5)*iza(p1,p7)*iza(p2,p5)*
     &    izb(p1,p7)*izb(p2,p6)**2*izb(p2,p7)*zab2(p5,p2,p6,p2)*zab2(p7
     &    ,p1,p7,p5) - za(p2,p4)*zb(p1,p2)*zb(p4,p5)*iza(p1,p7)*iza(p2,
     &    p5)*izb(p1,p7)*izb(p2,p6)**2*izb(p2,p7)*zab2(p3,p2,p6,p2)*
     &    zab2(p7,p1,p7,p4) + 2*za(p2,p6)*za(p3,p5)*zb(p1,p2)*zb(p4,p5)
     &    *iza(p1,p7)*iza(p2,p5)*izb(p1,p7)*izb(p2,p6)*izb(p2,p7)*zab2(
     &    p7,p1,p7,p5) - zb(p1,p2)*iza(p1,p7)*iza(p2,p5)*izb(p1,p7)*
     &    izb(p2,p6)**2*izb(p2,p7)*zab2(p2,p3,p4,p5)*zab2(p3,p2,p6,p2)*
     &    zab2(p7,p1,p7,p4) )
      a(2,1,1) = a(2,1,1) + prop34**(-1)*s167**(-1) * (  - za(p2,p3)*
     &    zb(p1,p2)*iza(p1,p7)*iza(p2,p5)*izb(p1,p7)*izb(p2,p6)*izb(p2,
     &    p7)*zab2(p2,p3,p4,p4)*zab2(p6,p1,p7,p5)*zab2(p7,p1,p7,p2)*
     &    s234**(-1) + za(p2,p3)*zb(p1,p2)*iza(p2,p5)*iza(p6,p7)*izb(p2
     &    ,p6)*izb(p2,p7)*izb(p6,p7)*zab2(p2,p3,p4,p4)*zab2(p6,p1,p7,p5
     &    )*zab2(p7,p6,p7,p2)*s234**(-1) + za(p2,p3)*zb(p1,p2)*iza(p2,
     &    p5)*izb(p2,p6)*izb(p6,p7)*zab2(p2,p3,p4,p4)*zab2(p7,p1,p6,p5)
     &    *s234**(-1) )
      a(2,1,1) = a(2,1,1) + prop34**(-1)*s167**(-1)*Qsum * ( za(p2,p3)*
     &    zb(p1,p2)*iza(p1,p7)*iza(p2,p5)*izb(p1,p7)*izb(p2,p6)*izb(p2,
     &    p7)*zab2(p2,p3,p4,p4)*zab2(p6,p1,p7,p5)*zab2(p7,p1,p7,p2)*
     &    s234**(-1) - za(p2,p3)*zb(p1,p2)*iza(p2,p5)*iza(p6,p7)*izb(p2
     &    ,p6)*izb(p2,p7)*izb(p6,p7)*zab2(p2,p3,p4,p4)*zab2(p6,p1,p7,p5
     &    )*zab2(p7,p6,p7,p2)*s234**(-1) - za(p2,p3)*zb(p1,p2)*iza(p2,
     &    p5)*izb(p2,p6)*izb(p6,p7)*zab2(p2,p3,p4,p4)*zab2(p7,p1,p6,p5)
     &    *s234**(-1) )
      a(2,1,1) = a(2,1,1) + prop34**(-1)*s267**(-1) * (  - za(p2,p6)*
     &    zb(p1,p4)*iza(p2,p5)*iza(p6,p7)*izb(p2,p6)*izb(p2,p7)*izb(p6,
     &    p7)*zab2(p2,p6,p7,p2)*zab2(p3,p1,p4,p5)*zab2(p7,p6,p7,p2)*
     &    s134**(-1) - za(p2,p6)*iza(p1,p5)*iza(p2,p5)*iza(p6,p7)*izb(
     &    p2,p6)*izb(p2,p7)*izb(p6,p7)*zab2(p2,p1,p5,p4)*zab2(p3,p6,p7,
     &    p2)*zab2(p7,p6,p7,p2) - za(p2,p7)*zb(p1,p4)*iza(p2,p5)*izb(p2
     &    ,p6)*izb(p6,p7)*zab2(p2,p6,p7,p2)*zab2(p3,p1,p4,p5)*
     &    s134**(-1) - za(p2,p7)*iza(p1,p5)*iza(p2,p5)*izb(p2,p6)*izb(
     &    p6,p7)*zab2(p2,p1,p5,p4)*zab2(p3,p6,p7,p2) - zb(p1,p4)*iza(p2
     &    ,p5)*izb(p2,p6)**2*izb(p2,p7)*zab2(p2,p6,p7,p2)*zab2(p3,p1,p4
     &    ,p5)*zab2(p7,p2,p6,p2)*s134**(-1) - iza(p1,p5)*iza(p2,p5)*
     &    izb(p2,p6)**2*izb(p2,p7)*zab2(p2,p1,p5,p4)*zab2(p3,p6,p7,p2)*
     &    zab2(p7,p2,p6,p2) )
      a(2,1,1) = a(2,1,1) + prop34**(-1)*s267**(-1)*Qsum * (  - za(p2,
     &    p6)*zb(p1,p4)*iza(p2,p5)*iza(p6,p7)*izb(p2,p6)*izb(p2,p7)*
     &    izb(p6,p7)*zab2(p2,p6,p7,p2)*zab2(p3,p1,p4,p5)*zab2(p7,p6,p7,
     &    p2)*s134**(-1) + za(p2,p6)*iza(p1,p5)*iza(p2,p5)*iza(p6,p7)*
     &    izb(p2,p6)*izb(p2,p7)*izb(p6,p7)*zab2(p2,p1,p5,p4)*zab2(p3,p6
     &    ,p7,p2)*zab2(p7,p6,p7,p2) - za(p2,p7)*zb(p1,p4)*iza(p2,p5)*
     &    izb(p2,p6)*izb(p6,p7)*zab2(p2,p6,p7,p2)*zab2(p3,p1,p4,p5)*
     &    s134**(-1) + za(p2,p7)*iza(p1,p5)*iza(p2,p5)*izb(p2,p6)*izb(
     &    p6,p7)*zab2(p2,p1,p5,p4)*zab2(p3,p6,p7,p2) - zb(p1,p4)*iza(p2
     &    ,p5)*izb(p2,p6)**2*izb(p2,p7)*zab2(p2,p6,p7,p2)*zab2(p3,p1,p4
     &    ,p5)*zab2(p7,p2,p6,p2)*s134**(-1) + iza(p1,p5)*iza(p2,p5)*
     &    izb(p2,p6)**2*izb(p2,p7)*zab2(p2,p1,p5,p4)*zab2(p3,p6,p7,p2)*
     &    zab2(p7,p2,p6,p2) )
      a(2,1,1) = a(2,1,1) + prop34**(-1) * (  - za(p2,p3)*zb(p1,p2)*
     &    iza(p1,p7)*iza(p2,p5)*izb(p1,p7)*izb(p2,p6)*izb(p2,p7)*zab2(
     &    p6,p2,p3,p4)*zab2(p7,p1,p7,p5)*zab3(p2,p1,p5,p7,p2)*
     &    s234**(-1)*s157**(-1) - za(p2,p3)*iza(p1,p5)*iza(p2,p5)*iza(
     &    p6,p7)*izb(p2,p6)*izb(p2,p7)*izb(p6,p7)*zab2(p2,p1,p5,p2)*
     &    zab2(p6,p2,p3,p4)*zab2(p7,p6,p7,p2)*s234**(-1) - za(p2,p3)*
     &    iza(p1,p5)*iza(p2,p5)*izb(p2,p6)*izb(p2,p7)*zab2(p2,p1,p5,p2)
     &    *zab2(p6,p2,p3,p4)*zab2(p7,p1,p5,p2)*s234**(-1)*s157**(-1) -
     &    za(p2,p3)*iza(p1,p5)*iza(p2,p5)*izb(p2,p6)*izb(p6,p7)*zab2(p2
     &    ,p1,p5,p2)*zab2(p7,p2,p3,p4)*s234**(-1) - za(p2,p6)*zb(p1,p2)
     &    *iza(p1,p7)*iza(p2,p5)*izb(p1,p7)*izb(p2,p6)*izb(p2,p7)*zab2(
     &    p3,p2,p6,p5)*zab2(p7,p1,p7,p4)*s256**(-1) + za(p2,p6)*zb(p1,
     &    p4)*iza(p2,p5)*izb(p2,p6)*izb(p2,p7)*zab2(p3,p1,p4,p2)*zab2(
     &    p7,p2,p6,p5)*s134**(-1)*s256**(-1) + zb(p1,p2)*iza(p1,p7)*
     &    iza(p2,p5)*izb(p1,p7)*izb(p2,p6)**2*izb(p2,p7)*zab2(p3,p2,p6,
     &    p2)*zab2(p7,p1,p7,p5)*zab3(p2,p1,p5,p7,p4)*s157**(-1) )
      a(2,1,1) = a(2,1,1) + prop34**(-1) * ( iza(p1,p5)*iza(p2,p5)*izb(
     &    p2,p6)**2*izb(p2,p7)*zab2(p2,p1,p5,p2)*zab2(p3,p2,p6,p2)*
     &    zab2(p7,p1,p5,p4)*s157**(-1) )
      a(2,1,1) = a(2,1,1) + prop34**(-1)*Qsum * ( za(p2,p3)*zb(p1,p2)*
     &    iza(p1,p7)*iza(p2,p5)*izb(p1,p7)*izb(p2,p6)*izb(p2,p7)*zab2(
     &    p6,p2,p3,p4)*zab2(p7,p1,p7,p5)*zab3(p2,p1,p5,p7,p2)*
     &    s234**(-1)*s157**(-1) + za(p2,p3)*iza(p1,p5)*iza(p2,p5)*iza(
     &    p6,p7)*izb(p2,p6)*izb(p2,p7)*izb(p6,p7)*zab2(p2,p1,p5,p2)*
     &    zab2(p6,p2,p3,p4)*zab2(p7,p6,p7,p2)*s234**(-1) + za(p2,p3)*
     &    iza(p1,p5)*iza(p2,p5)*izb(p2,p6)*izb(p2,p7)*zab2(p2,p1,p5,p2)
     &    *zab2(p6,p2,p3,p4)*zab2(p7,p1,p5,p2)*s234**(-1)*s157**(-1) +
     &    za(p2,p3)*iza(p1,p5)*iza(p2,p5)*izb(p2,p6)*izb(p6,p7)*zab2(p2
     &    ,p1,p5,p2)*zab2(p7,p2,p3,p4)*s234**(-1) - za(p2,p6)*zb(p1,p2)
     &    *iza(p1,p7)*iza(p2,p5)*izb(p1,p7)*izb(p2,p6)*izb(p2,p7)*zab2(
     &    p3,p2,p6,p5)*zab2(p7,p1,p7,p4)*s256**(-1) + za(p2,p6)*zb(p1,
     &    p4)*iza(p2,p5)*izb(p2,p6)*izb(p2,p7)*zab2(p3,p1,p4,p2)*zab2(
     &    p7,p2,p6,p5)*s134**(-1)*s256**(-1) - zb(p1,p2)*iza(p1,p7)*
     &    iza(p2,p5)*izb(p1,p7)*izb(p2,p6)**2*izb(p2,p7)*zab2(p3,p2,p6,
     &    p2)*zab2(p7,p1,p7,p5)*zab3(p2,p1,p5,p7,p4)*s157**(-1) )
      a(2,1,1) = a(2,1,1) + prop34**(-1)*Qsum * (  - iza(p1,p5)*iza(p2,
     &    p5)*izb(p2,p6)**2*izb(p2,p7)*zab2(p2,p1,p5,p2)*zab2(p3,p2,p6,
     &    p2)*zab2(p7,p1,p5,p4)*s157**(-1) )
      a(2,1,1) = a(2,1,1) + prop345**(-1)*s167**(-1) * ( 2*za(p2,p3)**2
     &    *zb(p1,p2)*iza(p1,p7)*iza(p2,p5)*iza(p3,p5)*izb(p1,p7)*izb(p2
     &    ,p6)*izb(p2,p7)*zab2(p6,p1,p7,p4)*zab2(p7,p1,p7,p2) - 2*za(p2
     &    ,p3)**2*zb(p1,p2)*iza(p2,p5)*iza(p3,p5)*iza(p6,p7)*izb(p2,p6)
     &    *izb(p2,p7)*izb(p6,p7)*zab2(p6,p1,p7,p4)*zab2(p7,p6,p7,p2) -
     &    2*za(p2,p3)**2*zb(p1,p2)*iza(p2,p5)*iza(p3,p5)*izb(p2,p6)*
     &    izb(p6,p7)*zab2(p7,p1,p6,p4) )
      a(2,1,1) = a(2,1,1) + prop345**(-1)*s267**(-1) * ( 2*za(p2,p3)*
     &    za(p2,p6)*zb(p1,p4)*iza(p2,p5)*iza(p3,p5)*iza(p6,p7)*izb(p2,
     &    p6)*izb(p2,p7)*izb(p6,p7)*zab2(p3,p6,p7,p2)*zab2(p7,p6,p7,p2)
     &     + 2*za(p2,p3)*za(p2,p7)*zb(p1,p4)*iza(p2,p5)*iza(p3,p5)*izb(
     &    p2,p6)*izb(p6,p7)*zab2(p3,p6,p7,p2) + 2*za(p2,p3)*zb(p1,p4)*
     &    iza(p2,p5)*iza(p3,p5)*izb(p2,p6)**2*izb(p2,p7)*zab2(p3,p6,p7,
     &    p2)*zab2(p7,p2,p6,p2) )
      a(2,1,1) = a(2,1,1) + prop345**(-1) * ( 2*za(p2,p3)*zb(p1,p2)*
     &    iza(p1,p7)*iza(p2,p5)*iza(p3,p5)*izb(p1,p7)*izb(p2,p6)**2*
     &    izb(p2,p7)*zab2(p3,p2,p6,p2)*zab2(p7,p1,p7,p4) )
      a(2,1,2)= + prop34**(-1)*prop345**(-1)*s167**(-1) * ( za(p1,p6)*
     &    za(p2,p3)*za(p3,p7)*zb(p1,p7)*zb(p3,p5)*iza(p1,p7)*iza(p5,p7)
     &    *iza(p6,p7)*izb(p6,p7)*zab2(p6,p1,p7,p4) - za(p1,p6)*za(p2,p3
     &    )*za(p3,p7)*zb(p1,p7)*zb(p4,p5)*iza(p1,p7)*iza(p5,p7)*iza(p6,
     &    p7)*izb(p6,p7)*zab2(p6,p1,p7,p3) + za(p1,p6)*za(p2,p3)*za(p4,
     &    p7)*zb(p1,p7)*zb(p4,p5)*iza(p1,p7)*iza(p5,p7)*iza(p6,p7)*izb(
     &    p6,p7)*zab2(p6,p1,p7,p4) - za(p1,p6)*za(p2,p3)*zb(p1,p7)*iza(
     &    p1,p7)*iza(p5,p7)*iza(p6,p7)*izb(p6,p7)*zab2(p6,p1,p7,p4)*
     &    zab2(p7,p3,p4,p5) - za(p1,p6)*za(p2,p4)*za(p3,p7)*zb(p1,p7)*
     &    zb(p4,p5)*iza(p1,p7)*iza(p5,p7)*iza(p6,p7)*izb(p6,p7)*zab2(p6
     &    ,p1,p7,p4) + za(p1,p6)*za(p2,p5)*za(p3,p7)*zb(p1,p7)*zb(p4,p5
     &    )*iza(p1,p7)*iza(p5,p7)*iza(p6,p7)*izb(p6,p7)*zab2(p6,p1,p7,
     &    p5) - 2*za(p1,p6)*za(p2,p7)*za(p3,p5)*zb(p1,p7)*zb(p4,p5)*
     &    iza(p1,p7)*iza(p5,p7)*iza(p6,p7)*izb(p6,p7)*zab2(p6,p1,p7,p5)
     &     )
      a(2,1,2) = a(2,1,2) + prop34**(-1)*prop345**(-1)*s267**(-1) * ( 2
     &    *za(p3,p5)*zb(p1,p5)*zb(p4,p5)*iza(p5,p7)*iza(p6,p7)*izb(p2,
     &    p6)*izb(p6,p7)*zab2(p6,p2,p6,p7)*zab2(p7,p2,p6,p7) + za(p3,p7
     &    )*zb(p1,p3)*zb(p4,p5)*iza(p5,p7)*iza(p6,p7)*izb(p2,p6)*izb(p6
     &    ,p7)*zab2(p3,p2,p6,p7)*zab2(p6,p2,p6,p7) - za(p3,p7)*zb(p1,p4
     &    )*zb(p3,p5)*iza(p5,p7)*iza(p6,p7)*izb(p2,p6)*izb(p6,p7)*zab2(
     &    p3,p2,p6,p7)*zab2(p6,p2,p6,p7) + za(p3,p7)*zb(p1,p4)*zb(p4,p5
     &    )*iza(p5,p7)*iza(p6,p7)*izb(p2,p6)*izb(p6,p7)*zab2(p4,p2,p6,
     &    p7)*zab2(p6,p2,p6,p7) - za(p3,p7)*zb(p1,p5)*zb(p4,p5)*iza(p5,
     &    p7)*iza(p6,p7)*izb(p2,p6)*izb(p6,p7)*zab2(p5,p2,p6,p7)*zab2(
     &    p6,p2,p6,p7) - za(p4,p7)*zb(p1,p4)*zb(p4,p5)*iza(p5,p7)*iza(
     &    p6,p7)*izb(p2,p6)*izb(p6,p7)*zab2(p3,p2,p6,p7)*zab2(p6,p2,p6,
     &    p7) + zb(p1,p4)*iza(p5,p7)*iza(p6,p7)*izb(p2,p6)*izb(p6,p7)*
     &    zab2(p3,p2,p6,p7)*zab2(p6,p2,p6,p7)*zab2(p7,p3,p4,p5) )
      a(2,1,2) = a(2,1,2) + prop34**(-1)*prop345**(-1) * ( 2*za(p3,p5)*
     &    zb(p4,p5)*iza(p1,p7)*iza(p5,p7)*iza(p6,p7)*izb(p2,p6)*izb(p6,
     &    p7)*zab2(p6,p1,p7,p5)*zab2(p7,p2,p6,p7) - za(p3,p7)*zb(p3,p5)
     &    *iza(p1,p7)*iza(p5,p7)*iza(p6,p7)*izb(p2,p6)*izb(p6,p7)*zab2(
     &    p3,p2,p6,p7)*zab2(p6,p1,p7,p4) + za(p3,p7)*zb(p4,p5)*iza(p1,
     &    p7)*iza(p5,p7)*iza(p6,p7)*izb(p2,p6)*izb(p6,p7)*zab2(p3,p2,p6
     &    ,p7)*zab2(p6,p1,p7,p3) + za(p3,p7)*zb(p4,p5)*iza(p1,p7)*iza(
     &    p5,p7)*iza(p6,p7)*izb(p2,p6)*izb(p6,p7)*zab2(p4,p2,p6,p7)*
     &    zab2(p6,p1,p7,p4) - za(p3,p7)*zb(p4,p5)*iza(p1,p7)*iza(p5,p7)
     &    *iza(p6,p7)*izb(p2,p6)*izb(p6,p7)*zab2(p5,p2,p6,p7)*zab2(p6,
     &    p1,p7,p5) - za(p4,p7)*zb(p4,p5)*iza(p1,p7)*iza(p5,p7)*iza(p6,
     &    p7)*izb(p2,p6)*izb(p6,p7)*zab2(p3,p2,p6,p7)*zab2(p6,p1,p7,p4)
     &     + iza(p1,p7)*iza(p5,p7)*iza(p6,p7)*izb(p2,p6)*izb(p6,p7)*
     &    zab2(p3,p2,p6,p7)*zab2(p6,p1,p7,p4)*zab2(p7,p3,p4,p5) )
      a(2,1,2) = a(2,1,2) + prop34**(-1)*s167**(-1) * ( za(p1,p6)*za(p2
     &    ,p3)*za(p2,p7)*zb(p1,p7)*iza(p1,p7)*iza(p2,p5)*iza(p5,p7)*
     &    iza(p6,p7)*izb(p6,p7)*zab2(p6,p1,p7,p4) - za(p1,p6)*za(p2,p3)
     &    *zb(p1,p7)*iza(p1,p7)*iza(p5,p7)*iza(p6,p7)*izb(p6,p7)*zab2(
     &    p6,p1,p7,p5)*zab2(p7,p2,p3,p4)*s234**(-1) )
      a(2,1,2) = a(2,1,2) + prop34**(-1)*s167**(-1)*Qsum * ( za(p1,p6)*
     &    za(p2,p3)*za(p2,p7)*zb(p1,p7)*iza(p1,p7)*iza(p2,p5)*iza(p5,p7
     &    )*iza(p6,p7)*izb(p6,p7)*zab2(p6,p1,p7,p4) + za(p1,p6)*za(p2,
     &    p3)*zb(p1,p7)*iza(p1,p7)*iza(p5,p7)*iza(p6,p7)*izb(p6,p7)*
     &    zab2(p6,p1,p7,p5)*zab2(p7,p2,p3,p4)*s234**(-1) )
      a(2,1,2) = a(2,1,2) + prop34**(-1)*s267**(-1) * ( zb(p1,p4)*iza(
     &    p5,p7)*iza(p6,p7)*izb(p2,p6)*izb(p6,p7)*zab2(p3,p1,p4,p5)*
     &    zab2(p6,p2,p6,p7)*zab2(p7,p2,p6,p7)*s134**(-1) + iza(p1,p5)*
     &    iza(p5,p7)*iza(p6,p7)*izb(p2,p6)*izb(p6,p7)*zab2(p3,p2,p6,p7)
     &    *zab2(p6,p2,p6,p7)*zab2(p7,p1,p5,p4) )
      a(2,1,2) = a(2,1,2) + prop34**(-1)*s267**(-1)*Qsum * ( zb(p1,p4)*
     &    iza(p5,p7)*iza(p6,p7)*izb(p2,p6)*izb(p6,p7)*zab2(p3,p1,p4,p5)
     &    *zab2(p6,p2,p6,p7)*zab2(p7,p2,p6,p7)*s134**(-1) - iza(p1,p5)*
     &    iza(p5,p7)*iza(p6,p7)*izb(p2,p6)*izb(p6,p7)*zab2(p3,p2,p6,p7)
     &    *zab2(p6,p2,p6,p7)*zab2(p7,p1,p5,p4) )
      a(2,1,2) = a(2,1,2) + prop34**(-1) * ( za(p2,p3)*iza(p1,p5)*iza(
     &    p5,p7)*iza(p6,p7)*izb(p6,p7)*zab2(p6,p1,p5,p7)*zab2(p6,p2,p3,
     &    p4)*zab2(p7,p1,p5,p7)*s234**(-1)*s157**(-1) + za(p2,p3)*iza(
     &    p1,p7)*iza(p5,p7)*iza(p6,p7)*izb(p6,p7)*zab2(p6,p1,p7,p5)*
     &    zab2(p6,p2,p3,p4)*zab2(p7,p1,p5,p7)*s234**(-1)*s157**(-1) -
     &    za(p2,p6)*za(p2,p7)*zb(p1,p4)*iza(p2,p5)*iza(p5,p7)*iza(p6,p7
     &    )*izb(p6,p7)*zab2(p3,p1,p4,p7)*zab2(p6,p2,p5,p7)*s134**(-1)*
     &    s256**(-1) + za(p2,p6)*za(p2,p7)*iza(p1,p7)*iza(p2,p5)*iza(p5
     &    ,p7)*iza(p6,p7)*izb(p6,p7)*zab2(p6,p1,p7,p4)*zab3(p3,p2,p5,p6
     &    ,p7)*s256**(-1) + zb(p1,p4)*iza(p5,p7)*iza(p6,p7)*izb(p2,p6)*
     &    izb(p6,p7)*zab2(p3,p1,p4,p7)*zab2(p6,p2,p5,p5)*zab2(p7,p2,p6,
     &    p7)*s134**(-1)*s256**(-1) - iza(p1,p5)*iza(p5,p7)*iza(p6,p7)*
     &    izb(p2,p6)*izb(p6,p7)*zab2(p3,p2,p6,p7)*zab2(p7,p1,p5,p7)*
     &    zab3(p6,p1,p5,p7,p4)*s157**(-1) - iza(p1,p7)*iza(p5,p7)*iza(
     &    p6,p7)*izb(p2,p6)*izb(p6,p7)*zab2(p3,p2,p6,p5)*zab2(p6,p1,p7,
     &    p4)*zab2(p7,p2,p6,p7)*s256**(-1) )
      a(2,1,2) = a(2,1,2) + prop34**(-1) * (  - iza(p1,p7)*iza(p5,p7)*
     &    iza(p6,p7)*izb(p2,p6)*izb(p6,p7)*zab2(p3,p2,p6,p7)*zab2(p6,p1
     &    ,p7,p5)*zab2(p7,p1,p5,p4)*s157**(-1) )
      a(2,1,2) = a(2,1,2) + prop34**(-1)*Qsum * (  - za(p2,p3)*iza(p1,
     &    p5)*iza(p5,p7)*iza(p6,p7)*izb(p6,p7)*zab2(p6,p1,p5,p7)*zab2(
     &    p6,p2,p3,p4)*zab2(p7,p1,p5,p7)*s234**(-1)*s157**(-1) - za(p2,
     &    p3)*iza(p1,p7)*iza(p5,p7)*iza(p6,p7)*izb(p6,p7)*zab2(p6,p1,p7
     &    ,p5)*zab2(p6,p2,p3,p4)*zab2(p7,p1,p5,p7)*s234**(-1)*
     &    s157**(-1) - za(p2,p6)*za(p2,p7)*zb(p1,p4)*iza(p2,p5)*iza(p5,
     &    p7)*iza(p6,p7)*izb(p6,p7)*zab2(p3,p1,p4,p7)*zab2(p6,p2,p5,p7)
     &    *s134**(-1)*s256**(-1) + za(p2,p6)*za(p2,p7)*iza(p1,p7)*iza(
     &    p2,p5)*iza(p5,p7)*iza(p6,p7)*izb(p6,p7)*zab2(p6,p1,p7,p4)*
     &    zab3(p3,p2,p5,p6,p7)*s256**(-1) + zb(p1,p4)*iza(p5,p7)*iza(p6
     &    ,p7)*izb(p2,p6)*izb(p6,p7)*zab2(p3,p1,p4,p7)*zab2(p6,p2,p5,p5
     &    )*zab2(p7,p2,p6,p7)*s134**(-1)*s256**(-1) + iza(p1,p5)*iza(p5
     &    ,p7)*iza(p6,p7)*izb(p2,p6)*izb(p6,p7)*zab2(p3,p2,p6,p7)*zab2(
     &    p7,p1,p5,p7)*zab3(p6,p1,p5,p7,p4)*s157**(-1) - iza(p1,p7)*
     &    iza(p5,p7)*iza(p6,p7)*izb(p2,p6)*izb(p6,p7)*zab2(p3,p2,p6,p5)
     &    *zab2(p6,p1,p7,p4)*zab2(p7,p2,p6,p7)*s256**(-1) )
      a(2,1,2) = a(2,1,2) + prop34**(-1)*Qsum * ( iza(p1,p7)*iza(p5,p7)
     &    *iza(p6,p7)*izb(p2,p6)*izb(p6,p7)*zab2(p3,p2,p6,p7)*zab2(p6,
     &    p1,p7,p5)*zab2(p7,p1,p5,p4)*s157**(-1) )
      a(2,1,2) = a(2,1,2) + prop345**(-1)*s167**(-1) * (  - 2*za(p1,p6)
     &    *za(p2,p3)*za(p3,p7)*zb(p1,p7)*iza(p1,p7)*iza(p3,p5)*iza(p5,
     &    p7)*iza(p6,p7)*izb(p6,p7)*zab2(p6,p1,p7,p4) )
      a(2,1,2) = a(2,1,2) + prop345**(-1)*s267**(-1) * ( 2*za(p3,p7)*
     &    zb(p1,p4)*iza(p3,p5)*iza(p5,p7)*iza(p6,p7)*izb(p2,p6)*izb(p6,
     &    p7)*zab2(p3,p2,p6,p7)*zab2(p6,p2,p6,p7) )
      a(2,1,2) = a(2,1,2) + prop345**(-1) * ( 2*za(p3,p7)*iza(p1,p7)*
     &    iza(p3,p5)*iza(p5,p7)*iza(p6,p7)*izb(p2,p6)*izb(p6,p7)*zab2(
     &    p3,p2,p6,p7)*zab2(p6,p1,p7,p4) )
      a(2,2,1)= + prop34**(-1)*prop345**(-1)*s167**(-1) * (  - za(p2,p3
     &    )*za(p3,p7)*zb(p1,p6)*zb(p3,p5)*iza(p1,p7)*iza(p5,p7)*iza(p6,
     &    p7)*izb(p1,p7)*izb(p6,p7)*zab2(p7,p1,p6,p4)*zab2(p7,p1,p7,p6)
     &     + za(p2,p3)*za(p3,p7)*zb(p1,p6)*zb(p3,p5)*iza(p5,p7)*iza(p6,
     &    p7)**2*izb(p6,p7)**2*zab2(p7,p1,p6,p4)*zab2(p7,p6,p7,p6) +
     &    za(p2,p3)*za(p3,p7)*zb(p1,p6)*zb(p4,p5)*iza(p1,p7)*iza(p5,p7)
     &    *iza(p6,p7)*izb(p1,p7)*izb(p6,p7)*zab2(p7,p1,p6,p3)*zab2(p7,
     &    p1,p7,p6) - za(p2,p3)*za(p3,p7)*zb(p1,p6)*zb(p4,p5)*iza(p5,p7
     &    )*iza(p6,p7)**2*izb(p6,p7)**2*zab2(p7,p1,p6,p3)*zab2(p7,p6,p7
     &    ,p6) - za(p2,p3)*za(p4,p7)*zb(p1,p6)*zb(p4,p5)*iza(p1,p7)*
     &    iza(p5,p7)*iza(p6,p7)*izb(p1,p7)*izb(p6,p7)*zab2(p7,p1,p6,p4)
     &    *zab2(p7,p1,p7,p6) + za(p2,p3)*za(p4,p7)*zb(p1,p6)*zb(p4,p5)*
     &    iza(p5,p7)*iza(p6,p7)**2*izb(p6,p7)**2*zab2(p7,p1,p6,p4)*
     &    zab2(p7,p6,p7,p6) + za(p2,p3)*zb(p1,p6)*iza(p1,p7)*iza(p5,p7)
     &    *iza(p6,p7)*izb(p1,p7)*izb(p6,p7)*zab2(p7,p1,p6,p4)*zab2(p7,
     &    p1,p7,p6)*zab2(p7,p3,p4,p5) )
      a(2,2,1) = a(2,2,1) + prop34**(-1)*prop345**(-1)*s167**(-1) * (
     &     - za(p2,p3)*zb(p1,p6)*iza(p5,p7)*iza(p6,p7)**2*izb(p6,p7)**2
     &    *zab2(p7,p1,p6,p4)*zab2(p7,p3,p4,p5)*zab2(p7,p6,p7,p6) + za(
     &    p2,p4)*za(p3,p7)*zb(p1,p6)*zb(p4,p5)*iza(p1,p7)*iza(p5,p7)*
     &    iza(p6,p7)*izb(p1,p7)*izb(p6,p7)*zab2(p7,p1,p6,p4)*zab2(p7,p1
     &    ,p7,p6) - za(p2,p4)*za(p3,p7)*zb(p1,p6)*zb(p4,p5)*iza(p5,p7)*
     &    iza(p6,p7)**2*izb(p6,p7)**2*zab2(p7,p1,p6,p4)*zab2(p7,p6,p7,
     &    p6) - za(p2,p5)*za(p3,p7)*zb(p1,p6)*zb(p4,p5)*iza(p1,p7)*iza(
     &    p5,p7)*iza(p6,p7)*izb(p1,p7)*izb(p6,p7)*zab2(p7,p1,p6,p5)*
     &    zab2(p7,p1,p7,p6) + za(p2,p5)*za(p3,p7)*zb(p1,p6)*zb(p4,p5)*
     &    iza(p5,p7)*iza(p6,p7)**2*izb(p6,p7)**2*zab2(p7,p1,p6,p5)*
     &    zab2(p7,p6,p7,p6) + 2*za(p2,p7)*za(p3,p5)*zb(p1,p6)*zb(p4,p5)
     &    *iza(p1,p7)*iza(p5,p7)*iza(p6,p7)*izb(p1,p7)*izb(p6,p7)*zab2(
     &    p7,p1,p6,p5)*zab2(p7,p1,p7,p6) - 2*za(p2,p7)*za(p3,p5)*zb(p1,
     &    p6)*zb(p4,p5)*iza(p5,p7)*iza(p6,p7)**2*izb(p6,p7)**2*zab2(p7,
     &    p1,p6,p5)*zab2(p7,p6,p7,p6) )
      a(2,2,1) = a(2,2,1) + prop34**(-1)*prop345**(-1)*s267**(-1) * (
     &     - 2*za(p2,p7)**2*za(p3,p5)*zb(p1,p5)*zb(p4,p5)*iza(p2,p6)*
     &    iza(p5,p7)*iza(p6,p7)*izb(p6,p7)*zab2(p7,p2,p6,p6) - za(p2,p7
     &    )**2*za(p3,p7)*zb(p1,p3)*zb(p4,p5)*iza(p2,p6)*iza(p5,p7)*iza(
     &    p6,p7)*izb(p6,p7)*zab2(p3,p2,p7,p6) + za(p2,p7)**2*za(p3,p7)*
     &    zb(p1,p4)*zb(p3,p5)*iza(p2,p6)*iza(p5,p7)*iza(p6,p7)*izb(p6,
     &    p7)*zab2(p3,p2,p7,p6) - za(p2,p7)**2*za(p3,p7)*zb(p1,p4)*zb(
     &    p4,p5)*iza(p2,p6)*iza(p5,p7)*iza(p6,p7)*izb(p6,p7)*zab2(p4,p2
     &    ,p7,p6) + za(p2,p7)**2*za(p3,p7)*zb(p1,p5)*zb(p4,p5)*iza(p2,
     &    p6)*iza(p5,p7)*iza(p6,p7)*izb(p6,p7)*zab2(p5,p2,p7,p6) + za(
     &    p2,p7)**2*za(p4,p7)*zb(p1,p4)*zb(p4,p5)*iza(p2,p6)*iza(p5,p7)
     &    *iza(p6,p7)*izb(p6,p7)*zab2(p3,p2,p7,p6) - za(p2,p7)**2*zb(p1
     &    ,p4)*iza(p2,p6)*iza(p5,p7)*iza(p6,p7)*izb(p6,p7)*zab2(p3,p2,
     &    p7,p6)*zab2(p7,p3,p4,p5) + 2*za(p2,p7)*za(p3,p5)*zb(p1,p5)*
     &    zb(p4,p5)*iza(p5,p7)*iza(p6,p7)**2*izb(p6,p7)**2*zab2(p7,p2,
     &    p6,p6)*zab2(p7,p6,p7,p6) )
      a(2,2,1) = a(2,2,1) + prop34**(-1)*prop345**(-1)*s267**(-1) * (
     &    za(p2,p7)*za(p3,p7)*zb(p1,p3)*zb(p4,p5)*iza(p5,p7)*iza(p6,p7)
     &    **2*izb(p6,p7)**2*zab2(p3,p2,p7,p6)*zab2(p7,p6,p7,p6) - za(p2
     &    ,p7)*za(p3,p7)*zb(p1,p4)*zb(p3,p5)*iza(p5,p7)*iza(p6,p7)**2*
     &    izb(p6,p7)**2*zab2(p3,p2,p7,p6)*zab2(p7,p6,p7,p6) + za(p2,p7)
     &    *za(p3,p7)*zb(p1,p4)*zb(p4,p5)*iza(p5,p7)*iza(p6,p7)**2*izb(
     &    p6,p7)**2*zab2(p4,p2,p7,p6)*zab2(p7,p6,p7,p6) - za(p2,p7)*za(
     &    p3,p7)*zb(p1,p5)*zb(p4,p5)*iza(p5,p7)*iza(p6,p7)**2*izb(p6,p7
     &    )**2*zab2(p5,p2,p7,p6)*zab2(p7,p6,p7,p6) - za(p2,p7)*za(p4,p7
     &    )*zb(p1,p4)*zb(p4,p5)*iza(p5,p7)*iza(p6,p7)**2*izb(p6,p7)**2*
     &    zab2(p3,p2,p7,p6)*zab2(p7,p6,p7,p6) + za(p2,p7)*zb(p1,p4)*
     &    iza(p5,p7)*iza(p6,p7)**2*izb(p6,p7)**2*zab2(p3,p2,p7,p6)*
     &    zab2(p7,p3,p4,p5)*zab2(p7,p6,p7,p6) )
      a(2,2,1) = a(2,2,1) + prop34**(-1)*prop345**(-1) * ( za(p2,p3)*
     &    za(p2,p7)*za(p3,p7)*zb(p1,p6)*zb(p3,p5)*iza(p1,p7)*iza(p2,p6)
     &    *iza(p5,p7)*iza(p6,p7)*izb(p1,p7)*izb(p6,p7)*zab2(p7,p1,p7,p4
     &    ) + za(p2,p3)*za(p2,p7)*za(p4,p7)*zb(p1,p6)*zb(p4,p5)*iza(p1,
     &    p7)*iza(p2,p6)*iza(p5,p7)*iza(p6,p7)*izb(p1,p7)*izb(p6,p7)*
     &    zab2(p7,p1,p7,p4) - za(p2,p3)*za(p2,p7)*zb(p1,p6)*iza(p1,p7)*
     &    iza(p2,p6)*iza(p5,p7)*iza(p6,p7)*izb(p1,p7)*izb(p6,p7)*zab2(
     &    p7,p1,p7,p4)*zab2(p7,p3,p4,p5) - 2*za(p2,p7)**2*za(p3,p5)*zb(
     &    p1,p6)*zb(p4,p5)*iza(p1,p7)*iza(p2,p6)*iza(p5,p7)*iza(p6,p7)*
     &    izb(p1,p7)*izb(p6,p7)*zab2(p7,p1,p7,p5) + za(p2,p7)*za(p3,p2)
     &    *za(p3,p7)*zb(p1,p6)*zb(p2,p6)*zb(p4,p5)*iza(p1,p7)*iza(p2,p6
     &    )*iza(p5,p7)*iza(p6,p7)*izb(p1,p7)*izb(p2,p6)*izb(p6,p7)*
     &    zab2(p7,p1,p7,p3) + za(p2,p7)*za(p3,p7)*za(p4,p2)*zb(p1,p6)*
     &    zb(p2,p6)*zb(p4,p5)*iza(p1,p7)*iza(p2,p6)*iza(p5,p7)*iza(p6,
     &    p7)*izb(p1,p7)*izb(p2,p6)*izb(p6,p7)*zab2(p7,p1,p7,p4) - za(
     &    p2,p7)*za(p3,p7)*za(p5,p2)*zb(p1,p6)*zb(p2,p6)*zb(p4,p5)*iza(
     &    p1,p7)*iza(p2,p6)*iza(p5,p7)*iza(p6,p7)*izb(p1,p7)*izb(p2,p6)
     &    *izb(p6,p7)*zab2(p7,p1,p7,p5) )
      a(2,2,1) = a(2,2,1) + prop34**(-1)*s167**(-1) * (  - za(p2,p3)*
     &    za(p2,p7)*zb(p1,p6)*iza(p1,p7)*iza(p2,p5)*iza(p5,p7)*iza(p6,
     &    p7)*izb(p1,p7)*izb(p6,p7)*zab2(p7,p1,p6,p4)*zab2(p7,p1,p7,p6)
     &     + za(p2,p3)*za(p2,p7)*zb(p1,p6)*iza(p2,p5)*iza(p5,p7)*iza(p6
     &    ,p7)**2*izb(p6,p7)**2*zab2(p7,p1,p6,p4)*zab2(p7,p6,p7,p6) +
     &    za(p2,p3)*zb(p1,p6)*iza(p1,p7)*iza(p5,p7)*iza(p6,p7)*izb(p1,
     &    p7)*izb(p6,p7)*zab2(p7,p1,p6,p5)*zab2(p7,p1,p7,p6)*zab2(p7,p2
     &    ,p3,p4)*s234**(-1) - za(p2,p3)*zb(p1,p6)*iza(p5,p7)*iza(p6,p7
     &    )**2*izb(p6,p7)**2*zab2(p7,p1,p6,p5)*zab2(p7,p2,p3,p4)*zab2(
     &    p7,p6,p7,p6)*s234**(-1) )
      a(2,2,1) = a(2,2,1) + prop34**(-1)*s167**(-1)*Qsum * (  - za(p2,
     &    p3)*za(p2,p7)*zb(p1,p6)*iza(p1,p7)*iza(p2,p5)*iza(p5,p7)*iza(
     &    p6,p7)*izb(p1,p7)*izb(p6,p7)*zab2(p7,p1,p6,p4)*zab2(p7,p1,p7,
     &    p6) + za(p2,p3)*za(p2,p7)*zb(p1,p6)*iza(p2,p5)*iza(p5,p7)*
     &    iza(p6,p7)**2*izb(p6,p7)**2*zab2(p7,p1,p6,p4)*zab2(p7,p6,p7,
     &    p6) - za(p2,p3)*zb(p1,p6)*iza(p1,p7)*iza(p5,p7)*iza(p6,p7)*
     &    izb(p1,p7)*izb(p6,p7)*zab2(p7,p1,p6,p5)*zab2(p7,p1,p7,p6)*
     &    zab2(p7,p2,p3,p4)*s234**(-1) + za(p2,p3)*zb(p1,p6)*iza(p5,p7)
     &    *iza(p6,p7)**2*izb(p6,p7)**2*zab2(p7,p1,p6,p5)*zab2(p7,p2,p3,
     &    p4)*zab2(p7,p6,p7,p6)*s234**(-1) )
      a(2,2,1) = a(2,2,1) + prop34**(-1)*s267**(-1) * (  - za(p2,p7)**2
     &    *zb(p1,p4)*iza(p2,p6)*iza(p5,p7)*iza(p6,p7)*izb(p6,p7)*zab2(
     &    p3,p1,p4,p5)*zab2(p7,p2,p6,p6)*s134**(-1) - za(p2,p7)**2*iza(
     &    p1,p5)*iza(p2,p6)*iza(p5,p7)*iza(p6,p7)*izb(p6,p7)*zab2(p3,p2
     &    ,p7,p6)*zab2(p7,p1,p5,p4) + za(p2,p7)*zb(p1,p4)*iza(p5,p7)*
     &    iza(p6,p7)**2*izb(p6,p7)**2*zab2(p3,p1,p4,p5)*zab2(p7,p2,p6,
     &    p6)*zab2(p7,p6,p7,p6)*s134**(-1) + za(p2,p7)*iza(p1,p5)*iza(
     &    p5,p7)*iza(p6,p7)**2*izb(p6,p7)**2*zab2(p3,p2,p7,p6)*zab2(p7,
     &    p1,p5,p4)*zab2(p7,p6,p7,p6) )
      a(2,2,1) = a(2,2,1) + prop34**(-1)*s267**(-1)*Qsum * (  - za(p2,
     &    p7)**2*zb(p1,p4)*iza(p2,p6)*iza(p5,p7)*iza(p6,p7)*izb(p6,p7)*
     &    zab2(p3,p1,p4,p5)*zab2(p7,p2,p6,p6)*s134**(-1) + za(p2,p7)**2
     &    *iza(p1,p5)*iza(p2,p6)*iza(p5,p7)*iza(p6,p7)*izb(p6,p7)*zab2(
     &    p3,p2,p7,p6)*zab2(p7,p1,p5,p4) + za(p2,p7)*zb(p1,p4)*iza(p5,
     &    p7)*iza(p6,p7)**2*izb(p6,p7)**2*zab2(p3,p1,p4,p5)*zab2(p7,p2,
     &    p6,p6)*zab2(p7,p6,p7,p6)*s134**(-1) - za(p2,p7)*iza(p1,p5)*
     &    iza(p5,p7)*iza(p6,p7)**2*izb(p6,p7)**2*zab2(p3,p2,p7,p6)*
     &    zab2(p7,p1,p5,p4)*zab2(p7,p6,p7,p6) )
      a(2,2,1) = a(2,2,1) + prop34**(-1) * ( za(p2,p3)*za(p2,p7)*zb(p1,
     &    p6)*iza(p1,p7)*iza(p2,p6)*iza(p5,p7)*iza(p6,p7)*izb(p1,p7)*
     &    izb(p6,p7)*zab2(p7,p1,p5,p4)*zab2(p7,p1,p7,p5)*s157**(-1) +
     &    za(p2,p3)*za(p2,p7)*iza(p1,p5)*iza(p2,p6)*iza(p5,p7)*iza(p6,
     &    p7)*izb(p6,p7)*zab2(p7,p1,p5,p4)*zab2(p7,p1,p5,p6)*s157**(-1)
     &     + za(p2,p3)*zb(p1,p6)*iza(p1,p7)*iza(p5,p7)*iza(p6,p7)*izb(
     &    p1,p7)*izb(p6,p7)*zab2(p7,p1,p5,p6)*zab2(p7,p1,p7,p5)*zab2(p7
     &    ,p2,p3,p4)*s234**(-1)*s157**(-1) + za(p2,p3)*iza(p1,p5)*iza(
     &    p5,p7)*iza(p6,p7)**2*izb(p6,p7)**2*zab2(p7,p1,p5,p6)*zab2(p7,
     &    p2,p3,p4)*zab2(p7,p6,p7,p6)*s234**(-1) + za(p2,p3)*iza(p1,p5)
     &    *iza(p5,p7)*iza(p6,p7)*izb(p6,p7)*zab2(p7,p1,p5,p6)**2*zab2(
     &    p7,p2,p3,p4)*s234**(-1)*s157**(-1) + za(p2,p7)**2*zb(p1,p4)*
     &    iza(p2,p5)*iza(p5,p7)*iza(p6,p7)**2*izb(p6,p7)**2*zab2(p3,p1,
     &    p4,p6)*zab2(p7,p6,p7,p6)*s134**(-1) - za(p2,p7)**2*zb(p1,p4)*
     &    iza(p2,p5)*iza(p5,p7)*iza(p6,p7)*izb(p6,p7)*zab2(p3,p1,p4,p6)
     &    *zab2(p7,p2,p5,p6)*s134**(-1)*s256**(-1) )
      a(2,2,1) = a(2,2,1) + prop34**(-1) * (  - za(p2,p7)**2*zb(p1,p4)*
     &    iza(p2,p6)*iza(p5,p7)*iza(p6,p7)*izb(p6,p7)*zab2(p3,p1,p4,p6)
     &    *zab2(p7,p2,p6,p5)*s134**(-1)*s256**(-1) + za(p2,p7)**2*zb(p1
     &    ,p6)*iza(p1,p7)*iza(p2,p5)*iza(p5,p7)*iza(p6,p7)*izb(p1,p7)*
     &    izb(p6,p7)*zab2(p3,p2,p5,p6)*zab2(p7,p1,p7,p4)*s256**(-1) +
     &    za(p2,p7)**2*zb(p1,p6)*iza(p1,p7)*iza(p2,p6)*iza(p5,p7)*iza(
     &    p6,p7)*izb(p1,p7)*izb(p6,p7)*zab2(p3,p2,p6,p5)*zab2(p7,p1,p7,
     &    p4)*s256**(-1) )
      a(2,2,1) = a(2,2,1) + prop34**(-1)*Qsum * (  - za(p2,p3)*za(p2,p7
     &    )*zb(p1,p6)*iza(p1,p7)*iza(p2,p6)*iza(p5,p7)*iza(p6,p7)*izb(
     &    p1,p7)*izb(p6,p7)*zab2(p7,p1,p5,p4)*zab2(p7,p1,p7,p5)*
     &    s157**(-1) - za(p2,p3)*za(p2,p7)*iza(p1,p5)*iza(p2,p6)*iza(p5
     &    ,p7)*iza(p6,p7)*izb(p6,p7)*zab2(p7,p1,p5,p4)*zab2(p7,p1,p5,p6
     &    )*s157**(-1) - za(p2,p3)*zb(p1,p6)*iza(p1,p7)*iza(p5,p7)*iza(
     &    p6,p7)*izb(p1,p7)*izb(p6,p7)*zab2(p7,p1,p5,p6)*zab2(p7,p1,p7,
     &    p5)*zab2(p7,p2,p3,p4)*s234**(-1)*s157**(-1) - za(p2,p3)*iza(
     &    p1,p5)*iza(p5,p7)*iza(p6,p7)**2*izb(p6,p7)**2*zab2(p7,p1,p5,
     &    p6)*zab2(p7,p2,p3,p4)*zab2(p7,p6,p7,p6)*s234**(-1) - za(p2,p3
     &    )*iza(p1,p5)*iza(p5,p7)*iza(p6,p7)*izb(p6,p7)*zab2(p7,p1,p5,
     &    p6)**2*zab2(p7,p2,p3,p4)*s234**(-1)*s157**(-1) + za(p2,p7)**2
     &    *zb(p1,p4)*iza(p2,p5)*iza(p5,p7)*iza(p6,p7)**2*izb(p6,p7)**2*
     &    zab2(p3,p1,p4,p6)*zab2(p7,p6,p7,p6)*s134**(-1) - za(p2,p7)**2
     &    *zb(p1,p4)*iza(p2,p5)*iza(p5,p7)*iza(p6,p7)*izb(p6,p7)*zab2(
     &    p3,p1,p4,p6)*zab2(p7,p2,p5,p6)*s134**(-1)*s256**(-1) )
      a(2,2,1) = a(2,2,1) + prop34**(-1)*Qsum * (  - za(p2,p7)**2*zb(p1
     &    ,p4)*iza(p2,p6)*iza(p5,p7)*iza(p6,p7)*izb(p6,p7)*zab2(p3,p1,
     &    p4,p6)*zab2(p7,p2,p6,p5)*s134**(-1)*s256**(-1) + za(p2,p7)**2
     &    *zb(p1,p6)*iza(p1,p7)*iza(p2,p5)*iza(p5,p7)*iza(p6,p7)*izb(p1
     &    ,p7)*izb(p6,p7)*zab2(p3,p2,p5,p6)*zab2(p7,p1,p7,p4)*
     &    s256**(-1) + za(p2,p7)**2*zb(p1,p6)*iza(p1,p7)*iza(p2,p6)*
     &    iza(p5,p7)*iza(p6,p7)*izb(p1,p7)*izb(p6,p7)*zab2(p3,p2,p6,p5)
     &    *zab2(p7,p1,p7,p4)*s256**(-1) )
      a(2,2,1) = a(2,2,1) + prop345**(-1)*s167**(-1) * ( 2*za(p2,p3)*
     &    za(p3,p7)*zb(p1,p6)*iza(p1,p7)*iza(p3,p5)*iza(p5,p7)*iza(p6,
     &    p7)*izb(p1,p7)*izb(p6,p7)*zab2(p7,p1,p6,p4)*zab2(p7,p1,p7,p6)
     &     - 2*za(p2,p3)*za(p3,p7)*zb(p1,p6)*iza(p3,p5)*iza(p5,p7)*iza(
     &    p6,p7)**2*izb(p6,p7)**2*zab2(p7,p1,p6,p4)*zab2(p7,p6,p7,p6) )
      a(2,2,1) = a(2,2,1) + prop345**(-1)*s267**(-1) * (  - 2*za(p2,p7)
     &    **2*za(p3,p7)*zb(p1,p4)*iza(p2,p6)*iza(p3,p5)*iza(p5,p7)*iza(
     &    p6,p7)*izb(p6,p7)*zab2(p3,p2,p7,p6) + 2*za(p2,p7)*za(p3,p7)*
     &    zb(p1,p4)*iza(p3,p5)*iza(p5,p7)*iza(p6,p7)**2*izb(p6,p7)**2*
     &    zab2(p3,p2,p7,p6)*zab2(p7,p6,p7,p6) )
      a(2,2,1) = a(2,2,1) + prop345**(-1) * ( 2*za(p2,p7)*za(p3,p2)*za(
     &    p3,p7)*zb(p1,p6)*zb(p2,p6)*iza(p1,p7)*iza(p2,p6)*iza(p3,p5)*
     &    iza(p5,p7)*iza(p6,p7)*izb(p1,p7)*izb(p2,p6)*izb(p6,p7)*zab2(
     &    p7,p1,p7,p4) )
      a(2,2,2)= + prop34**(-1)*prop345**(-1)*s167**(-1) * ( za(p2,p3)**
     &    2*zb(p1,p6)*zb(p3,p5)*iza(p2,p5)*iza(p2,p7)*iza(p6,p7)*zab3(
     &    p2,p1,p6,p7,p4) - za(p2,p3)**2*zb(p1,p6)*zb(p4,p5)*iza(p2,p5)
     &    *iza(p2,p7)*iza(p6,p7)*zab3(p2,p1,p6,p7,p3) + za(p2,p3)**2*
     &    zb(p1,p7)*zb(p3,p5)*iza(p2,p5)*iza(p2,p6)*iza(p6,p7)*zab3(p2,
     &    p1,p6,p7,p4) - za(p2,p3)**2*zb(p1,p7)*zb(p4,p5)*iza(p2,p5)*
     &    iza(p2,p6)*iza(p6,p7)*zab3(p2,p1,p6,p7,p3) - za(p2,p3)**2*zb(
     &    p3,p5)*iza(p1,p7)*iza(p2,p5)*iza(p2,p6)*iza(p2,p7)*zab2(p2,p1
     &    ,p7,p6)*zab3(p2,p1,p6,p7,p4) + za(p2,p3)**2*zb(p4,p5)*iza(p1,
     &    p7)*iza(p2,p5)*iza(p2,p6)*iza(p2,p7)*zab2(p2,p1,p7,p6)*zab3(
     &    p2,p1,p6,p7,p3) + za(p2,p3)*zb(p1,p6)*zb(p4,p5)*iza(p2,p7)*
     &    iza(p6,p7)*zab3(p2,p1,p6,p7,p5) + za(p2,p3)*zb(p1,p6)*iza(p2,
     &    p5)*iza(p2,p7)*iza(p6,p7)*zab2(p2,p3,p4,p5)*zab3(p2,p1,p6,p7,
     &    p4) + za(p2,p3)*zb(p1,p7)*zb(p4,p5)*iza(p2,p6)*iza(p6,p7)*
     &    zab3(p2,p1,p6,p7,p5) + za(p2,p3)*zb(p1,p7)*iza(p2,p5)*iza(p2,
     &    p6)*iza(p6,p7)*zab2(p2,p3,p4,p5)*zab3(p2,p1,p6,p7,p4) )
      a(2,2,2) = a(2,2,2) + prop34**(-1)*prop345**(-1)*s167**(-1) * (
     &     - za(p2,p3)*zb(p4,p5)*iza(p1,p7)*iza(p2,p6)*iza(p2,p7)*zab2(
     &    p2,p1,p7,p6)*zab3(p2,p1,p6,p7,p5) - za(p2,p3)*iza(p1,p7)*iza(
     &    p2,p5)*iza(p2,p6)*iza(p2,p7)*zab2(p2,p1,p7,p6)*zab2(p2,p3,p4,
     &    p5)*zab3(p2,p1,p6,p7,p4) )
      a(2,2,2) = a(2,2,2) + prop34**(-1)*s167**(-1) * ( za(p2,p3)*zb(p1
     &    ,p6)*iza(p2,p5)*iza(p2,p7)*iza(p6,p7)*zab2(p2,p3,p4,p4)*zab3(
     &    p2,p1,p6,p7,p5)*s234**(-1) + za(p2,p3)*zb(p1,p7)*iza(p2,p5)*
     &    iza(p2,p6)*iza(p6,p7)*zab2(p2,p3,p4,p4)*zab3(p2,p1,p6,p7,p5)*
     &    s234**(-1) - za(p2,p3)*iza(p1,p7)*iza(p2,p5)*iza(p2,p6)*iza(
     &    p2,p7)*zab2(p2,p1,p7,p6)*zab2(p2,p3,p4,p4)*zab3(p2,p1,p6,p7,
     &    p5)*s234**(-1) )
      a(2,2,2) = a(2,2,2) + prop34**(-1)*s167**(-1)*Qsum * (  - za(p2,
     &    p3)*zb(p1,p6)*iza(p2,p5)*iza(p2,p7)*iza(p6,p7)*zab2(p2,p3,p4,
     &    p4)*zab3(p2,p1,p6,p7,p5)*s234**(-1) - za(p2,p3)*zb(p1,p7)*
     &    iza(p2,p5)*iza(p2,p6)*iza(p6,p7)*zab2(p2,p3,p4,p4)*zab3(p2,p1
     &    ,p6,p7,p5)*s234**(-1) + za(p2,p3)*iza(p1,p7)*iza(p2,p5)*iza(
     &    p2,p6)*iza(p2,p7)*zab2(p2,p1,p7,p6)*zab2(p2,p3,p4,p4)*zab3(p2
     &    ,p1,p6,p7,p5)*s234**(-1) )
      a(2,2,2) = a(2,2,2) + prop34**(-1) * (  - za(p2,p3)*iza(p1,p5)*
     &    iza(p2,p5)*iza(p2,p6)*iza(p2,p7)*zab2(p2,p1,p5,p7)*zab2(p2,p3
     &    ,p4,p4)*zab3(p2,p1,p5,p7,p6)*s234**(-1)*s157**(-1) - za(p2,p3
     &    )*iza(p1,p5)*iza(p2,p5)*iza(p2,p6)*iza(p6,p7)*zab2(p2,p1,p5,
     &    p7)*zab2(p2,p3,p4,p4)*s234**(-1) - za(p2,p3)*iza(p1,p5)*iza(
     &    p2,p5)*iza(p2,p7)*iza(p6,p7)*zab2(p2,p1,p5,p6)*zab2(p2,p3,p4,
     &    p4)*s234**(-1) - za(p2,p3)*iza(p1,p7)*iza(p2,p5)*iza(p2,p6)*
     &    iza(p2,p7)*zab2(p2,p1,p7,p5)*zab2(p2,p3,p4,p4)*zab3(p2,p1,p5,
     &    p7,p6)*s234**(-1)*s157**(-1) )
      a(2,2,2) = a(2,2,2) + prop34**(-1)*Qsum * ( za(p2,p3)*iza(p1,p5)*
     &    iza(p2,p5)*iza(p2,p6)*iza(p2,p7)*zab2(p2,p1,p5,p7)*zab2(p2,p3
     &    ,p4,p4)*zab3(p2,p1,p5,p7,p6)*s234**(-1)*s157**(-1) + za(p2,p3
     &    )*iza(p1,p5)*iza(p2,p5)*iza(p2,p6)*iza(p6,p7)*zab2(p2,p1,p5,
     &    p7)*zab2(p2,p3,p4,p4)*s234**(-1) + za(p2,p3)*iza(p1,p5)*iza(
     &    p2,p5)*iza(p2,p7)*iza(p6,p7)*zab2(p2,p1,p5,p6)*zab2(p2,p3,p4,
     &    p4)*s234**(-1) + za(p2,p3)*iza(p1,p7)*iza(p2,p5)*iza(p2,p6)*
     &    iza(p2,p7)*zab2(p2,p1,p7,p5)*zab2(p2,p3,p4,p4)*zab3(p2,p1,p5,
     &    p7,p6)*s234**(-1)*s157**(-1) )
      a(2,2,2) = a(2,2,2) + prop345**(-1)*s167**(-1) * (  - 2*za(p2,p3)
     &    **2*zb(p1,p6)*iza(p2,p5)*iza(p2,p7)*iza(p3,p5)*iza(p6,p7)*
     &    zab3(p2,p1,p6,p7,p4) - 2*za(p2,p3)**2*zb(p1,p7)*iza(p2,p5)*
     &    iza(p2,p6)*iza(p3,p5)*iza(p6,p7)*zab3(p2,p1,p6,p7,p4) + 2*za(
     &    p2,p3)**2*iza(p1,p7)*iza(p2,p5)*iza(p2,p6)*iza(p2,p7)*iza(p3,
     &    p5)*zab2(p2,p1,p7,p6)*zab3(p2,p1,p6,p7,p4) )
      return
      end
