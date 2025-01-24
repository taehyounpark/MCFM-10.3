!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function qloop_c34x56m2(p1,p2,p3,p4,p5,p6,p7,za,zb)
c     This is an implementation of deLaurentis, Campbell and Ellis
c     arXiv:2203.17170, Eq(4.23)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer::p1,p2,p3,p4,p5,p6,p7
      complex(dp)::qloop_c34x56m2,qloop_c34x56m2_sym
      qloop_c34x56m2=qloop_c34x56m2_sym(p1,p2,p3,p4,p5,p6,p7,za,zb)
     &              +qloop_c34x56m2_sym(p1,p2,p5,p6,p3,p4,p7,za,zb)
      return
      end

      function qloop_c34x56m2_sym(p1,p2,p3,p4,p5,p6,p7,za,zb)
c symmetric part under 12 <-> 56
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer::p1,p2,p3,p4,p5,p6,p7
      complex(dp):: qloop_c34x56m2_sym
      real(dp)::s3,delta,s34x56
      complex(dp)::zab2,zba2,zbab22,zaba22,
     & za7x3456x3,za7x5634x7,zb7x5634x7,zb4x5634x7,za1x3456x1

c     statement functions
      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p3,p1)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zba2(p1,p2,p3,p4)=zb(p1,p2)*za(p2,p4)+zb(p1,p3)*za(p3,p4)
      zbab22(p1,p2,p3,p4,p5,p6)=
     & zb(p1,p2)*zab2(p2,p4,p5,p6)+zb(p1,p3)*zab2(p3,p4,p5,p6)
      zaba22(p1,p2,p3,p4,p5,p6)=
     & za(p1,p2)*zba2(p2,p4,p5,p6)+za(p1,p3)*zba2(p3,p4,p5,p6)
c     end statement functions

      delta =((s3(p1,p2,p7)-s(p3,p4)-s(p5,p6))**2-4*s(p3,p4)*s(p5,p6))
      zb7x5634x7=zbab22(p7,p5,p6,p3,p4,p7)
      zb4x5634x7=zbab22(p4,p5,p6,p3,p4,p7)
      za7x3456x3=zaba22(p7,p3,p4,p5,p6,p3)
      za7x5634x7=zaba22(p7,p5,p6,p3,p4,p7)
      za1x3456x1=zaba22(p1,p3,p4,p5,p6,p1)
      s34x56=s(p3,p5)+s(p3,p6)+s(p4,p5)+s(p4,p6)


      qloop_c34x56m2_sym=delta
     & *(zb(p6,p5)*zab2(p5,p3,p4,p7)*za(p7,p3)-zb(p6,p7)*za7x3456x3)
     & *(-zb(p2,p7)*zb(p3,p4)
     & *(za(p3,p1)*za(p5,p7)-za(p3,p7)*za(p1,p5))
     & -zb(p1,p2)*za(p1,p5)*zab2(p1,p6,p7,p4)
     & +zb(p1,p2)*za(p1,p5)*za(p1,p7)*zb(p4,p7))
     & /(s(p1,p2)*zb7x5634x7*za7x5634x7**2)

     & -za(p1,p6)*zab2(p5,p3,p4,p6)
     & *(zb(p6,p5)*zab2(p5,p3,p4,p7)*za(p7,p3)-zb(p6,p7)*za7x3456x3)
     & *(-2*(s3(p5,p6,p3)-s3(p5,p6,p4))
     & *(za(p1,p5)*zb(p5,p4)-za(p1,p7)*zb(p7,p4))
     & +4*zab2(p3,p5,p6,p4)*(za(p1,p5)*zb(p5,p3)-za(p1,p7)*zb(p7,p3)))
     & /(za(p1,p2)*zb7x5634x7*za7x5634x7**2)

     & -za(p1,p5)*(s3(p3,p4,p5)-s3(p3,p4,p6))
     & *(zb(p6,p5)*zab2(p5,p3,p4,p7)*za(p7,p3)-zb(p6,p7)*za7x3456x3)
     & *(-2*zab2(p1,p6,p7,p3)*zab2(p3,p5,p6,p4)
     & +zab2(p1,p6,p7,p4)*(s3(p5,p6,p3)-s3(p5,p6,p4)))
     & /(za(p1,p2)*zb7x5634x7*za7x5634x7**2)

     & +za(p1,p5)*zb(p2,p4)*za(p3,p7)*zab2(p7,p1,p2,p6)*delta
     & /(s(p1,p2)*za7x5634x7**2)

      qloop_c34x56m2_sym=qloop_c34x56m2_sym
     & +za(p1,p3)*zb(p4,p7)*delta
     & *(-za(p1,p5)*zb(p1,p6)*zb(p2,p7)
     & +2*za(p1,p5)*zb(p1,p7)*zb(p2,p6)+za(p2,p5)*zb(p2,p6)*zb(p2,p7)
     & +0.5_dp*zb(p1,p2)*za(p1,p5)*zb(p6,p7))
     & /(s(p1,p2)*za7x5634x7*zb7x5634x7)

     & -2*za(p1,p3)*zab2(p1,p5,p6,p4)*zab2(p5,p3,p4,p6)
     & /(za(p1,p2)*za7x5634x7)

     & -0.5_dp*zb(p1,p2)*za(p1,p3)*za(p1,p5)*zb(p4,p7)*zb(p6,p7)
     & *(s3(p3,p4,p5)-s3(p3,p4,p6))*(s3(p5,p6,p3)-s3(p5,p6,p4))
     & /(s(p1,p2)*za7x5634x7*zb7x5634x7)

     & +za(p1,p3)*zb(p4,p7)*zab2(p5,p3,p4,p6)
     & *(s3(p5,p6,p3)-s3(p5,p6,p4))
     & *(zb(p2,p7)*s34x56-2*zb(p2,p7)*zab2(p3,p5,p6,p3)
     & +2*zb(p1,p2)*za(p1,p5)*zb(p5,p7))
     & /(s(p1,p2)*za7x5634x7*zb7x5634x7)

     & +zb(p2,p7)*zab2(p3,p5,p6,p4)*zab2(p5,p3,p4,p6)*s34x56
     & *(za(p1,p2)*zb(p2,p7)-4*za(p1,p3)*zb(p3,p7))
     & /(s(p1,p2)*za7x5634x7*zb7x5634x7)

      qloop_c34x56m2_sym=qloop_c34x56m2_sym
     & +za(p1,p3)*zab2(p3,p5,p6,p4)*zab2(p5,p3,p4,p6)
     & *(4*zb(p2,p7)*zb(p3,p7)*(s(p3,p5)+s(p3,p6)-s(p3,p4))
     & -2*zb(p1,p2)*za(p1,p5)*zb(p3,p7)*zb(p5,p7)
     & -4*zb(p3,p4)*zb(p2,p7)*zab2(p4,p1,p2,p7))
     & /(s(p1,p2)*za7x5634x7*zb7x5634x7)

     & +2*za(p1,p7)**2*delta
     & *(zb(p6,p5)*zab2(p5,p3,p4,p7)*za(p7,p3)-zb(p6,p7)*za7x3456x3)
     & *(zb4x5634x7*za(p7,p5)-zb(p4,p7)*zab2(p7,p3,p4,p6)*za(p6,p5))
     & /(za(p1,p2)*zb7x5634x7*za7x5634x7**3)

     & +zb(p2,p7)**2*delta
     & *(zb(p6,p5)*zab2(p5,p3,p4,p7)*za(p7,p3)-zb(p6,p7)*za7x3456x3)
     & *(zb4x5634x7*za(p7,p5)-zb(p4,p7)*zab2(p7,p3,p4,p6)*za(p6,p5))
     & /(zb(p1,p2)*zb7x5634x7**2*za7x5634x7**2)

     & -za1x3456x1*(2*zab2(p3,p5,p6,p4)*zab2(p5,p3,p4,p6)*s34x56
     & +za(p3,p5)*zb(p4,p6)*delta)
     & /(za(p1,p2)*za7x5634x7*delta)

      qloop_c34x56m2_sym=qloop_c34x56m2_sym/(4*s(p3,p4)*s(p5,p6))

      return
      end

