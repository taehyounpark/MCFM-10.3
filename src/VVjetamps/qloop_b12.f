!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function qloop_b12(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
c     This is an implementation of deLaurentis, Campbell and Ellis
c     arXiv:2203.17170, Eq(4.42)
c     this function returns
c     qloop_b12 which is symmetric under (3,4)<-->(5,6)
c     i.e. (qloop_b12(1,2,3,4,5,6,7)
c          +qloop_b12(1,2,5,6,3,4,7))
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp):: qloop_b12,qloop_b12_symbit
      qloop_b12=
     & +qloop_b12_symbit(p1,p2,p3,p4,p5,p6,p7,za,zb)
     & +qloop_b12_symbit(p1,p2,p5,p6,p3,p4,p7,za,zb)
      return
      end

      function qloop_b12_symbit(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer::p1,p2,p3,p4,p5,p6,p7,p8
      real(dp)::s3,s567,s12,s34,Delta
      complex(dp)::qloop_b12_symbit,zab2,zba2,zbab22,zaba22,zbaba222,
     & zba21,ab7x12x7,ang7x3456x7,bab7x3456x7

c     statement functions
      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p3,p1)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zba2(p1,p2,p3,p4)=zb(p1,p2)*za(p2,p4)+zb(p1,p3)*za(p3,p4)
      zbab22(p1,p2,p3,p4,p5,p6)=
     & zb(p1,p2)*zab2(p2,p4,p5,p6)+zb(p1,p3)*zab2(p3,p4,p5,p6)
      zaba22(p1,p2,p3,p4,p5,p6)=
     & za(p1,p2)*zba2(p2,p4,p5,p6)+za(p1,p3)*zba2(p3,p4,p5,p6)
      zbaba222(p1,p2,p3,p4,p5,p6,p7,p8)=
     & +zb(p1,p2)*zaba22(p2,p4,p5,p6,p7,p8)
     & +zb(p1,p3)*zaba22(p3,p4,p5,p6,p7,p8)
c     end statement functions

      s567=s3(p5,p6,p7)
      s34=s(p3,p4)
      s12=s(p1,p2)
      Delta=(s567-s12-s34)**2-4*s12*s34
      ang7x3456x7=zaba22(p7,p3,p4,p5,p6,p7)
      bab7x3456x7=zbab22(p7,p3,p4,p5,p6,p7)
      ab7x12x7=zab2(p7,p1,p2,p7)
      zba21=
     & +zbaba222(p2,p1,p2,p3,p4,p5,p6,p1)
     & -zbaba222(p2,p5,p6,p3,p4,p1,p2,p1)

      qloop_b12_symbit=
     & -zb(p1,p2)*za(p1,p7)**2*za(p3,p7)**2*za(p5,p7)**2
     & *bab7x3456x7
     & /(2*za(p3,p4)*za(p5,p6)*ab7x12x7*ang7x3456x7**3)

     & -za(p1,p2)*zb(p2,p7)**2*za(p3,p7)**2*za(p5,p7)**2
     & /(2*za(p3,p4)*za(p5,p6)*ab7x12x7*ang7x3456x7**2)

     & +za(p5,p7)**2*zab2(p7,p1,p2,p4)**2*zba21
     & /(zb(p3,p4)*za(p5,p6)*ab7x12x7*ang7x3456x7**3)

     & +za(p1,p2)*zb(p1,p2)*za(p1,p7)*zb(p2,p7)
     & *zb(p4,p7)*za(p5,p7)**2*zab2(p7,p1,p2,p4)
     & /(2*zb(p3,p4)*za(p5,p6)*ab7x12x7**2*ang7x3456x7**2)

     & +zab2(p7,p1,p2,p4)*zb(p2,p7)*za(p5,p7)
     & *(-3*za(p1,p5)*zab2(p7,p1,p2,p4)
     & +za(p1,p7)*zab2(p5,p1,p2,p4))
     & /(2*zb(p3,p4)*za(p5,p6)*ab7x12x7*ang7x3456x7**2)

     & +zb(p2,p4)*za(p5,p7)
     & *(za(p1,p5)*zab2(p7,p1,p2,p4)
     & +za(p1,p7)*zab2(p5,p1,p2,p4))
     & /(2*zb(p3,p4)*za(p5,p6)*ang7x3456x7**2)

      qloop_b12_symbit=qloop_b12_symbit
     & +za(p1,p2)*zb(p2,p4)**2*zaba22(p5,p3,p4,p1,p2,p5)
     & /(2*zb(p3,p4)*za(p5,p6)*ang7x3456x7*Delta)

     & +za(p1,p2)*zb(p2,p4)*za(p5,p7)*zab2(p5,p1,p2,p4)
     & *zab2(p7,p3,p4,p2)*(s12-s567)
     & /(2*zb(p3,p4)*za(p5,p6)*ang7x3456x7**2*Delta)

     & -za(p1,p2)*zb(p1,p2)*za(p3,p5)**2*zab2(p1,p3,p4,p2)
     & /(2*za(p3,p4)*za(p5,p6)*ang7x3456x7*Delta)

     & -za(p1,p2)*zb(p1,p2)*za(p1,p7)*za(p3,p5)*za(p3,p7)
     & *zab2(p5,p3,p4,p2)*(s12-s567)
     & /(2*za(p3,p4)*za(p5,p6)*ang7x3456x7**2*Delta)

     & +za(p1,p5)*za(p3,p7)*(zb(p1,p2)*za(p1,p5)*za(p3,p7)
     & -2*zb(p1,p2)*za(p1,p7)*za(p3,p5)
     & +2*zb(p2,p4)*za(p3,p4)*za(p5,p7))
     & /(2*za(p3,p4)*za(p5,p6)*ang7x3456x7**2)

     & -3*zab2(p1,p3,p4,p2)*zab2(p3,p1,p2,p4)
     & *zaba22(p5,p3,p4,p1,p2,p5)*(s34-s567-s12)
     & /(za(p5,p6)*ang7x3456x7*Delta**2)

      qloop_b12_symbit=qloop_b12_symbit
     & +za(p1,p7)*zaba22(p5,p3,p4,p1,p2,p5)
     & *(za(p3,p7)*zb(p1,p2)
     & *(za(p1,p3)*zb(p3,p4)+2*za(p1,p2)*zb(p2,p4))
     & +zab2(p3,p1,p2,p4)
     & *(6*zb(p1,p2)*za(p1,p7)-7*zab2(p7,p3,p4,p2)))
     & /(2*za(p5,p6)*ang7x3456x7**2*Delta)

     & +zb(p1,p2)*za(p1,p2)*za(p1,p7)*za(p3,p5)
     & *(2*zb(p1,p2)*za(p1,p3)*zb(p3,p4)*za(p5,p7)
     & -2*zb(p1,p2)*za(p1,p5)*zb(p3,p4)*za(p3,p7)
     & -4*za(p1,p4)*zb(p1,p4)*zb(p2,p4)*za(p5,p7)
     & -4*zb(p1,p3)*za(p1,p5)*zb(p2,p4)*za(p3,p7)
     & -4*zb(p1,p4)*za(p1,p7)*zb(p2,p4)*za(p4,p5)
     & -4*zb(p2,p3)*zb(p2,p4)*za(p2,p5)*za(p3,p7)
     & -4*zb(p2,p4)**2*za(p2,p5)*za(p4,p7)
     & +za(p3,p7)*za(p4,p5)*zb(p2,p4)*zb(p3,p4)
     & -za(p3,p4)*za(p5,p7)*zb(p2,p4)*zb(p3,p4)
     & +za(p3,p5)*za(p3,p7)*zb(p2,p3)*zb(p3,p4))
     & /(2*za(p5,p6)*ang7x3456x7**2*Delta)

      qloop_b12_symbit=qloop_b12_symbit
     & +za(p3,p7)
     & *(zaba22(p1,p1,p2,p3,p4,p7)-zaba22(p1,p3,p4,p1,p2,p7))
     & *(-zb(p1,p2)*za(p1,p5)*zab2(p5,p1,p2,p4)
     & -2*zb(p2,p4)*zaba22(p5,p3,p4,p6,p7,p5))
     & /(2*za(p5,p6)*ang7x3456x7**2*Delta)

     & +(zaba22(p5,p1,p2,p3,p4,p7)-zaba22(p5,p3,p4,p1,p2,p7))
     & *(-4*zb(p1,p2)*za(p1,p3)*zb(p1,p4)*za(p1,p5)*za(p1,p7)
     & +2*za(p1,p3)*zb(p1,p4)*za(p1,p5)*zab2(p7,p3,p4,p2)
     & -zb(p1,p2)*za(p1,p3)*za(p1,p7)*zb(p2,p4)*za(p2,p5)
     & -za(p1,p3)*zb(p2,p3)*zb(p2,p4)*za(p2,p5)*za(p3,p7)
     & -za(p1,p3)*zb(p2,p4)**2*za(p2,p5)*za(p4,p7)
     & -3*za(p1,p5)*za(p2,p3)*zb(p2,p4)*zab2(p7,p5,p6,p2))
     & /(2*za(p5,p6)*ang7x3456x7**2*Delta)

     & +(8*zb(p2,p4)*za(p3,p5)
     & *(za(p1,p2)*zb(p2,p3)*za(p3,p5)+za(p1,p4)*zb(p2,p4)*za(p2,p5))
     & +zb(p1,p2)*za(p1,p5)**2*za(p2,p3)*zb(p2,p4)
     & -za(p1,p5)*za(p2,p4)*zb(p2,p4)**2*za(p3,p5)
     & +zb(p1,p2)*za(p1,p3)*zb(p1,p4)*za(p1,p5)**2
     & -2*zb(p1,p2)*za(p1,p3)*za(p1,p5)*zb(p3,p4)*za(p3,p5)
     & -7*zb(p2,p4)*za(p3,p5)*za(p1,p5)*(s(p2,p3)+s(p1,p3)+s(p1,p4))
     & +5*za(p1,p5)*zb(p2,p4)*za(p3,p4)*zab2(p5,p1,p2,p4))
     & /(2*za(p5,p6)*ang7x3456x7*Delta)

      return
      end
