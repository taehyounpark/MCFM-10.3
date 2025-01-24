!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function qloop_d12x34x56m2_sym(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
c     This is an implementation of deLaurentis, Campbell and Ellis
c     arXiv:2203.17170, Eq(4.9)
c     this function returns
c     part of qloop_d12x34x56m2 symmetric under (1,2)<-->(5,6)
c     i.e. +1/2*qloop_d12x34x56m2(1,2,3,4,5,6,7)
c          +1/2*qloop_d12x34x56m2(5,6,3,4,1,2,7))
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp):: qloop_d12x34x56m2_sym,qloop_d12x34x56m2_symbit
      qloop_d12x34x56m2_sym=
     & +qloop_d12x34x56m2_symbit(p1,p2,p3,p4,p5,p6,p7,za,zb)
     & +qloop_d12x34x56m2_symbit(p5,p6,p3,p4,p1,p2,p7,za,zb)
      return
      end

      function qloop_d12x34x56m2_symbit(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer p1,p2,p3,p4,p5,p6,p7,p8
      complex(dp):: qloop_d12x34x56m2_symbit,zab2,zaba22,zbab22,
     & zabab222,zbaba222,zba41,zba23,ang7x34x56x7,sqr7x34x56x7
c--- statement functions
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zaba22(p1,p2,p3,p4,p5,p6)=
     & zab2(p1,p2,p3,p4)*za(p4,p6)+zab2(p1,p2,p3,p5)*za(p5,p6)
      zbab22(p1,p2,p3,p4,p5,p6)=
     & zb(p1,p2)*zab2(p2,p4,p5,p6)+zb(p1,p3)*zab2(p3,p4,p5,p6)
      zabab222(p1,p2,p3,p4,p5,p6,p7,p8)=
     & +zaba22(p1,p2,p3,p4,p5,p6)*zb(p6,p8)
     & +zaba22(p1,p2,p3,p4,p5,p7)*zb(p7,p8)
      zbaba222(p1,p2,p3,p4,p5,p6,p7,p8)=
     & +zbab22(p1,p2,p3,p4,p5,p6)*za(p6,p8)
     & +zbab22(p1,p2,p3,p4,p5,p7)*za(p7,p8)
c--- end statement functions

      zba41=zbaba222(p4,p1,p2,p3,p4,p5,p6,p1)
     &     -zbaba222(p4,p5,p6,p3,p4,p1,p2,p1)
      zba23=zbaba222(p2,p1,p2,p3,p4,p5,p6,p3)
     &     -zbaba222(p2,p5,p6,p3,p4,p1,p2,p3)

      ang7x34x56x7=zaba22(p7,p3,p4,p5,p6,p7)
      sqr7x34x56x7=zbab22(p7,p3,p4,p5,p6,p7)

      qloop_d12x34x56m2_symbit=
     & -za(p5,p7)**2*zab2(p7,p1,p2,p4)*zaba22(p1,p3,p4,p5,p6,p7)*zba41
     & /(za(p1,p2)*zb(p3,p4)*za(p5,p6)*ang7x34x56x7**3)

     & -za(p1,p7)*zb(p6,p7)**2*zab2(p7,p5,p6,p4)*zba23
     & /(4*zb(p5,p6)*ang7x34x56x7**2*sqr7x34x56x7)

     & +zb(p6,p7)**2*zab2(p3,p1,p2,p7)*zbab22(p2,p3,p4,p5,p6,p7)
     & *zba23
     & /(4*zb(p1,p2)*za(p3,p4)*zb(p5,p6)*ang7x34x56x7*sqr7x34x56x7**2)

     & +za(p1,p3)*zb(p2,p7)*zb(p6,p7)**2*zab2(p3,p1,p2,p7)
     & /(8*s(p1,p2)*za(p3,p4)*zb(p5,p6)*sqr7x34x56x7)

     & +za(p1,p3)*zb(p4,p6)*za(p5,p7)*zaba22(p1,p3,p4,p5,p6,p7)
     & *zabab222(p7,p1,p2,p3,p4,p5,p6,p7)
     & /(4*za(p1,p2)*s(p3,p4)*s(p5,p6)*ang7x34x56x7**2)

     & +za(p1,p5)*zb(p4,p6)*zab2(p7,p5,p6,p4)*zab2(p7,p1,p2,p7)
     & *zaba22(p1,p3,p4,p5,p6,p7)
     & /(4*za(p1,p2)*zb(p3,p4)*s(p5,p6)*ang7x34x56x7**2)

     & -za(p1,p3)*za(p3,p5)*zb(p6,p7)*zaba22(p1,p3,p4,p5,p6,p7)
     & /(8*za(p1,p2)*za(p3,p4)*za(p5,p6)*zb(p5,p6)*ang7x34x56x7)

     & -za(p1,p5)*zb(p4,p6)*zb(p4,p7)*zaba22(p1,p3,p4,p5,p6,p7)
     & /(8*za(p1,p2)*zb(p3,p4)*za(p5,p6)*zb(p5,p6)*ang7x34x56x7)

     & +za(p1,p5)*zb(p2,p7)*zb(p4,p6)*zab2(p7,p5,p6,p4)
     & /(8*zb(p3,p4)*za(p5,p6)*zb(p5,p6)*ang7x34x56x7)

     & -zb(p2,p6)*zb(p2,p7)*za(p3,p7)**2*zb(p6,p7)
     & /(8*zb(p1,p2)*za(p3,p4)*zb(p5,p6)*ang7x34x56x7)

      qloop_d12x34x56m2_symbit=qloop_d12x34x56m2_symbit
     & -za(p3,p5)**2*zaba22(p1,p3,p4,p5,p6,p7)**2
     & /(4*za(p1,p2)*za(p3,p4)*za(p5,p6)*ang7x34x56x7**2)

     & -zb(p4,p6)**2*zaba22(p1,p3,p4,p5,p6,p7)**2
     & /(2*za(p1,p2)*zb(p3,p4)*zb(p5,p6)*ang7x34x56x7**2)

     & +zb(p2,p4)*za(p5,p7)*zab2(p5,p3,p6,p4)
     & *zaba22(p1,p3,p4,p5,p6,p7)
     & /(2*zb(p3,p4)*za(p5,p6)*ang7x34x56x7**2)

     & -za(p3,p5)*za(p1,p3)*zab2(p7,p1,p2,p6)
     & *zaba22(p1,p3,p4,p5,p6,p7)
     & /(2*za(p1,p2)*za(p3,p4)*ang7x34x56x7**2)

     & -zb(p2,p6)*za(p3,p7)**2*zb(p6,p7)*zaba22(p1,p3,p4,p5,p6,p7)
     & /(4*za(p3,p4)*zb(p5,p6)*ang7x34x56x7**2)

     & +za(p3,p5)*zaba22(p1,p3,p4,p5,p6,p7)
     & *(za(p1,p7)*zb(p3,p4)*za(p3,p5)-za(p1,p5)*zab2(p7,p1,p2,p4))
     & /(4*za(p1,p2)*za(p5,p6)*ang7x34x56x7**2)

     & +za(p1,p7)*za(p3,p4)*zb(p4,p6)**2*zaba22(p1,p3,p4,p5,p6,p7)
     & /(2*za(p1,p2)*zb(p5,p6)*ang7x34x56x7**2)

     & +za(p1,p5)*zb(p4,p6)*zab2(p7,p1,p2,p4)*zaba22(p1,p3,p4,p5,p6,p7)
     & /(2*za(p1,p2)*zb(p3,p4)*ang7x34x56x7**2)

     & -zb(p1,p2)*za(p1,p7)*zb(p4,p6)**2*zaba22(p1,p3,p4,p5,p6,p7)
     & /(2*zb(p3,p4)*zb(p5,p6)*ang7x34x56x7**2)

     & -za(p1,p3)*za(p1,p5)*zb(p4,p6)
     & /(za(p1,p2)*ang7x34x56x7)

      qloop_d12x34x56m2_symbit=qloop_d12x34x56m2_symbit
     & -3*za(p1,p3)*zb(p2,p6)*za(p3,p5)
     & /(4*za(p3,p4)*ang7x34x56x7)

     & -za(p1,p3)*(za(p1,p3)*zab2(p5,p1,p2,p6)
     & -za(p1,p3)*zab2(p5,p3,p4,p6)
     & +za(p1,p5)*za(p3,p5)*zb(p5,p6))
     & /(2*za(p1,p2)*za(p3,p4)*ang7x34x56x7)

     & -3*za(p1,p5)*za(p1,p3)
     & *(2*zab2(p5,p1,p2,p4)-zb(p3,p4)*za(p3,p5))
     & /(16*za(p1,p2)*za(p5,p6)*ang7x34x56x7)

     & -za(p1,p3)*zab2(p1,p5,p7,p6)*zab2(p3,p5,p7,p6)
     & /(2*za(p1,p2)*za(p3,p4)*zb(p5,p6)*ang7x34x56x7)

     & -za(p1,p5)*(8*za(p1,p2)*zb(p2,p4)*zb(p4,p6)*za(p5,p6)
     & +6*zb(p1,p4)*za(p1,p5)*za(p1,p6)*zb(p4,p6)
     & +3*zb(p1,p4)*za(p1,p5)**2*zb(p4,p5)
     & +3*za(p1,p5)*zb(p2,p4)*za(p2,p6)*zb(p4,p6))
     & /(16*za(p1,p2)*zb(p3,p4)*za(p5,p6)*ang7x34x56x7)

     & +zb(p4,p6)**2*zab2(p1,p3,p4,p2)
     & /(2*zb(p3,p4)*zb(p5,p6)*ang7x34x56x7)
      return
      end
