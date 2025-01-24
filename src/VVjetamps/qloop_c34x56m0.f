!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function qloop_c34x56m0(p1,p2,p3,p4,p5,p6,p7,za,zb)
c     This is an implementation of deLaurentis, Campbell and Ellis
c     arXiv:2203.17170, Eq(4.22)
c Wrapper function that implements 1234567 + 1256347 symmetry
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer::p1,p2,p3,p4,p5,p6,p7
      complex(dp)::qloop_c34x56m0,qloop_c34x56m0_symbit

      qloop_c34x56m0=
     &  qloop_c34x56m0_symbit(p1,p2,p3,p4,p5,p6,p7,za,zb)
     & +qloop_c34x56m0_symbit(p1,p2,p5,p6,p3,p4,p7,za,zb)
      return
      end


      function qloop_c34x56m0_symbit(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer::p1,p2,p3,p4,p5,p6,p7
      real(dp)::s3,s34x56,s127,s345,s346,s356,s456,delta
c    & ,delta735
      complex(dp)::qloop_c34x56m0_symbit,zab2,zba2,zaba22,
     & zaba7x1256x7,zaba1x3456x1,zaba1x3456x7,zaba1x5634x7,zaba5x3412x7

c     statement functions
      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p3,p1)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zba2(p1,p2,p3,p4)=zb(p1,p2)*za(p2,p4)+zb(p1,p3)*za(p3,p4)
      zaba22(p1,p2,p3,p4,p5,p6)=
     & za(p1,p2)*zba2(p2,p4,p5,p6)+za(p1,p3)*zba2(p3,p4,p5,p6)
c     end statement functions

      s127=s3(p1,p2,p7)
      s345=s3(p3,p4,p5)
      s356=s3(p3,p5,p6)
      s346=s3(p3,p4,p6)
      s456=s3(p4,p5,p6)
      delta=((s127-s(p3,p4)-s(p5,p6))**2-4._dp*s(p3,p4)*s(p5,p6))
      zaba1x3456x7=zaba22(p1,p3,p4,p5,p6,p7)
      zaba1x5634x7=zaba22(p1,p5,p6,p3,p4,p7)
      zaba7x1256x7=zaba22(p7,p1,p2,p5,p6,p7)
      zaba1x3456x1=zaba22(p1,p3,p4,p5,p6,p1)
      zaba5x3412x7=zaba22(p5,p3,p4,p1,p2,p7)
      s34x56=s(p3,p5)+s(p3,p6)+s(p4,p5)+s(p4,p6)

      qloop_c34x56m0_symbit=
     & za(p1,p7)*za(p3,p7)**2*zab2(p7,p1,p2,p6)**2*delta
     & *(2*za(p1,p7)*s34x56-zaba1x3456x7)
     & /(4*za(p3,p4)*zb(p5,p6)*zaba7x1256x7**4)

     & +za(p1,p7)*za(p3,p7)*zab2(p7,p1,p2,p6)
     & *s34x56**2
     & *(2*za(p1,p7)*za(p3,p5)*zb(p5,p6)-3*za(p1,p7)*za(p3,p4)*zb(p4,p6)
     & -za(p1,p3)*zab2(p7,p1,p2,p6))
     & /(4*za(p3,p4)*zb(p5,p6)*zaba7x1256x7**3)

     & +za(p1,p7)*za(p3,p7)*zab2(p7,p1,p2,p6)*s34x56
     & *(-2*za(p1,p5)*zab2(p7,p1,p2,p4)-za(p1,p7)*zab2(p5,p3,p6,p4))
     & /(2*zaba7x1256x7**3)

     & +za(p1,p7)*za(p3,p7)*zab2(p7,p1,p2,p6)*zaba1x3456x7
     & *(za(p5,p6)*zb(p6,p4)-za(p5,p3)*zb(p3,p4))
     & /(zaba7x1256x7**3)

     & +za(p1,p7)**2*s34x56
     & *(za(p3,p4)**2*zb(p4,p6)**2+za(p3,p5)**2*zb(p5,p6)**2)
     & /(4*za(p3,p4)*zb(p5,p6)*zaba7x1256x7**2)

      qloop_c34x56m0_symbit=qloop_c34x56m0_symbit
     & +s34x56
     & *(za(p1,p3)*zab2(p7,p1,p2,p6)+za(p1,p7)*za(p3,p4)*zb(p4,p6))**2
     & /(4*za(p3,p4)*zb(p5,p6)*zaba7x1256x7**2)

     & -3*zab2(p3,p5,p6,p4)*zab2(p5,p3,p4,p6)*zaba1x3456x1
     & *(s127-s(p3,p4)+s(p5,p6))*(s127+s(p3,p4)-s(p5,p6))
     & /(4*zaba7x1256x7*delta**2)

     & -5*za(p1,p7)**2*zab2(p3,p5,p6,p4)*zab2(p5,p3,p4,p6)
     & *(s127+s(p3,p4)-s(p5,p6))*(s127-s(p3,p4)+s(p5,p6))
     & /(8*zaba7x1256x7**2*delta)

     & -zab2(p3,p5,p6,p4)*zab2(p7,p1,p2,p6)
     & *zaba5x3412x7*zaba1x3456x1
     & /(zaba7x1256x7**2*delta)

     & +za(p1,p3)*za(p5,p7)*s127*(zaba1x3456x7-zaba1x5634x7)
     & *(2*zab2(p5,p3,p4,p6)*zb(p4,p5)-(s345-s346)*zb(p4,p6))
     & /(4*zaba7x1256x7**2*delta)

      qloop_c34x56m0_symbit=qloop_c34x56m0_symbit
     & +zb(p3,p4)*(zaba1x3456x7-zaba1x5634x7)
     & *(2*zab2(p3,p5,p6,p4)*za(p4,p5)-(s456-s356)*za(p3,p5))
     & *(6*za(p1,p3)*zab2(p7,p1,p2,p6)+3*za(p3,p4)*zb(p4,p6)*za(p1,p7))
     & /(4*zaba7x1256x7**2*delta)

     & +za(p1,p3)*(2*zab2(p5,p3,p4,p6)*zb(p4,p5)-(s345-s346)*zb(p4,p6))
     & *(2*za(p5,p6)*zab2(p1,p3,p4,p6)
     & -11*za(p1,p5)*s127-14*za(p1,p5)*za(p3,p4)*zb(p3,p4))
     & /(8*zaba7x1256x7*delta)

     & +(2*zab2(p3,p5,p6,p4)*za(p4,p5)-(s456-s356)*za(p3,p5))
     & *(5*zb(p4,p6)*zaba1x3456x1
     & -16*za(p1,p3)*zb(p3,p4)*zab2(p1,p3,p4,p6))
     & /(8*zaba7x1256x7*delta)

      qloop_c34x56m0_symbit=qloop_c34x56m0_symbit
     & +za(p1,p7)*(-12*za(p3,p5)*zb(p4,p6)*zaba1x3456x7
     & +4*za(p1,p5)*zab2(p7,p1,p2,p6)*zab2(p3,p2,p7,p4)
     & -2*za(p1,p3)*zb(p3,p4)*za(p3,p5)*zab2(p7,p3,p4,p6)
     & +3*za(p1,p7)*za(p3,p4)*zb(p4,p6)**2*za(p5,p6))
     & /(8*zaba7x1256x7**2)

     & +za(p1,p3)*zab2(p7,p1,p2,p6)
     & *(-4*zb(p1,p4)*za(p1,p5)*za(p1,p7)
     & -2*za(p1,p5)*zb(p2,p4)*za(p2,p7)
     & -9*za(p1,p7)*zb(p3,p4)*za(p3,p5))
     & /(4*zaba7x1256x7**2)

     & +za(p1,p3)*za(p1,p5)*zb(p4,p6)/(zaba7x1256x7)

      qloop_c34x56m0_symbit=qloop_c34x56m0_symbit/za(p1,p2)

      return
      end
