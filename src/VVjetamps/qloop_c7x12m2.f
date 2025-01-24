!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function qloop_c7x12m2(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
c     This is an implementation of deLaurentis, Campbell and Ellis
c     arXiv:2203.17170, Eq(4.32)
c     this function returns part of qloop_c7x12m2 symmetric under (1,2)<-->(3,4)
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp):: qloop_c7x12m2,qloop_c7x12m2_symbit
      qloop_c7x12m2=
     & +qloop_c7x12m2_symbit(p1,p2,p3,p4,p5,p6,p7,za,zb)
     & +qloop_c7x12m2_symbit(p1,p2,p5,p6,p3,p4,p7,za,zb)
      return
      end

      function qloop_c7x12m2_symbit(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer::p1,p2,p3,p4,p5,p6,p7
      complex(dp):: qloop_c7x12m2_symbit
      real(dp):: s34x56
      complex(dp)::zab2,zba2,zbab22,zaba22,
     & za7x5634x7,zb7x5634x7,za7x3456x3,za1x5612x5,za5x1234x7,zb4x1234x7
c     statement functions
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zba2(p1,p2,p3,p4)=zb(p1,p2)*za(p2,p4)+zb(p1,p3)*za(p3,p4)
      zbab22(p1,p2,p3,p4,p5,p6)=
     & zb(p1,p2)*zab2(p2,p4,p5,p6)+zb(p1,p3)*zab2(p3,p4,p5,p6)
      zaba22(p1,p2,p3,p4,p5,p6)=
     & za(p1,p2)*zba2(p2,p4,p5,p6)+za(p1,p3)*zba2(p3,p4,p5,p6)
c     end statement functions

      za7x5634x7=zaba22(p7,p5,p6,p3,p4,p7)
      zb7x5634x7=zbab22(p7,p5,p6,p3,p4,p7)
      za7x3456x3=zaba22(p7,p3,p4,p5,p6,p3)
      za1x5612x5=zaba22(p1,p5,p6,p1,p2,p5)
      za5x1234x7=zaba22(p5,p1,p2,p3,p4,p7)
      zb4x1234x7=zbab22(p4,p1,p2,p3,p4,p7)
      s34x56=s(p3,p5)+s(p3,p6)+s(p4,p5)+s(p4,p6)

      qloop_c7x12m2_symbit=
     & +4*za(p1,p7)*za(p3,p4)*za(p5,p7)
     & *zb(p5,p6)*zb(p2,p7)*zb(p4,p7)**2*zab2(p5,p1,p2,p7)
     & /(zab2(p7,p1,p2,p7)*za7x5634x7*zb7x5634x7)

     & +4*za(p1,p7)**2*za(p3,p7)*za(p5,p7)
     & *(zb(p6,p5)*zab2(p5,p3,p4,p7)*za(p7,p3)
     & -zb(p6,p7)*za7x3456x3)*s34x56*zb(p4,p3)
     & /(za(p1,p2)*za7x5634x7**3)

     & +za(p1,p7)*za(p5,p7)*zb(p2,p7)*zb(p4,p7)
     & *(zb(p6,p5)*zab2(p5,p3,p4,p7)*za(p7,p3)
     & -zb(p6,p7)*za7x3456x3)*s34x56
     & /(za7x5634x7**2*zb7x5634x7)

     & +2*za(p1,p7)**2*za(p3,p7)*zab2(p5,p3,p4,p6)
     & *(zb(p4,p7)*s34x56+zb4x1234x7)
     & /(za(p1,p2)*za7x5634x7**2)

     & +zb(p2,p6)*zb(p2,p7)*za(p3,p7)*za(p5,p7)
     & *zab2(p7,p1,p2,p4)*s34x56
     & /(zb(p1,p2)*za7x5634x7**2)

     & -za(p1,p7)*za(p3,p7)
     & *(2*zb(p2,p4)*zb(p6,p7)*za5x1234x7
     & +zb(p1,p2)*zb(p4,p6)*za1x5612x5)
     & /(za7x5634x7**2)

     & +za(p1,p7)*za(p3,p7)*zb(p1,p2)
     & *(zab2(p1,p3,p7,p4)*zab2(p5,p3,p4,p6)
     & -2*za(p1,p3)*zb(p3,p4)*zab2(p5,p3,p4,p6)
     & +za(p1,p2)*zb(p2,p4)*zab2(p5,p1,p2,p6)
     & -za(p1,p5)*zb(p5,p6)*zab2(p5,p3,p6,p4))
     & /(za7x5634x7**2)

      qloop_c7x12m2_symbit=qloop_c7x12m2_symbit
     & +2*zb(p2,p7)**2*za(p5,p6)*zb(p6,p7)*zab2(p3,p1,p2,p7)*zb(p4,p3)
     & *(zb(p6,p5)*zab2(p5,p3,p4,p7)*za(p7,p3)-zb(p6,p7)*za7x3456x3)
     & /(zb(p1,p2)*za7x5634x7*zb7x5634x7**2)

     & +zb(p2,p7)*zb(p4,p7)*za(p3,p7)*zb(p6,p5)
     & *(-2*zab2(p5,p3,p4,p2)*zab2(p5,p1,p2,p7)
     & +3*za(p1,p5)*zb(p1,p2)*zab2(p5,p1,p2,p7)
     & -za(p1,p5)**2*zb(p1,p2)*zb(p1,p7))
     & /(zb(p1,p2)*za7x5634x7*zb7x5634x7)

     & +zb(p2,p7)*zb(p4,p7)*za(p3,p7)
     & *(3*za(p1,p5)*zb(p6,p7)*s34x56
     & +2*zb(p6,p7)*za1x5612x5
     & +za(p1,p5)*za(p2,p5)*zb(p2,p7)*zb(p5,p6)
     & +6*za(p1,p5)*s(p5,p6)*zb(p6,p7)
     & -2*za(p5,p7)*zb(p6,p7)*zab2(p1,p5,p6,p7))
     & /(za7x5634x7*zb7x5634x7)

     & -zb(p2,p7)*za(p3,p7)*zb(p4,p6)
     & *(za(p5,p3)*zb(p3,p2)+za(p5,p4)*zb(p4,p2)-3*za(p5,p6)*zb(p6,p2))
     & /(zb(p1,p2)*za7x5634x7)

     & +za(p1,p3)*za(p5,p7)
     & *(zb(p2,p4)*zb(p6,p7)+3*zb(p2,p6)*zb(p4,p7))
     & /(za7x5634x7)

      qloop_c7x12m2_symbit=qloop_c7x12m2_symbit/(4*s(p3,p4)*s(p5,p6))

      return
      end
