!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function qloop_d12x34x56m2_asy(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
c     This is an implementation of deLaurentis, Campbell and Ellis
c     arXiv:2203.17170, Eq(4.10)
c     this function returns
c     part of qloop_d12x34x56m2 antisymmetric under (1,2)<-->(5,6)
c     i.e. +1/2*qloop_d12x34x56m2(1,2,3,4,5,6,7)
c          -1/2*qloop_d12x34x56m2(5,6,3,4,1,2,7))
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp):: qloop_d12x34x56m2_asy,qloop_d12x34x56m2_asybit
      qloop_d12x34x56m2_asy=
     & +qloop_d12x34x56m2_asybit(p1,p2,p3,p4,p5,p6,p7,za,zb)
     & -qloop_d12x34x56m2_asybit(p5,p6,p3,p4,p1,p2,p7,za,zb)
      return
      end

      function qloop_d12x34x56m2_asybit(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer p1,p2,p3,p4,p5,p6,p7,p8
      complex(dp):: qloop_d12x34x56m2_asybit,zab2,zaba22,zbab22,
     & zabab222,ang7x12x56x7,sqr7x12x56x7
c--- statement functions
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zaba22(p1,p2,p3,p4,p5,p6)=
     & zab2(p1,p2,p3,p4)*za(p4,p6)+zab2(p1,p2,p3,p5)*za(p5,p6)
      zbab22(p1,p2,p3,p4,p5,p6)=
     & zb(p1,p2)*zab2(p2,p4,p5,p6)+zb(p1,p3)*zab2(p3,p4,p5,p6)
      zabab222(p1,p2,p3,p4,p5,p6,p7,p8)=
     & +zaba22(p1,p2,p3,p4,p5,p6)*zb(p6,p8)
     & +zaba22(p1,p2,p3,p4,p5,p7)*zb(p7,p8)
c--- end statement functions

      ang7x12x56x7=zaba22(p7,p1,p2,p5,p6,p7)
      sqr7x12x56x7=zbab22(p7,p1,p2,p5,p6,p7)
      qloop_d12x34x56m2_asybit=
     & -0.25_dp*za(p1,p7)*za(p3,p5)*zb(p4,p6)*zaba22(p1,p3,p4,p5,p6,p7)
     & *zabab222(p7,p1,p2,p3,p4,p5,p6,p7)
     & /(za(p1,p2)*s(p3,p4)*s(p5,p6)*ang7x12x56x7**2)

     & +za(p1,p3)*zb(p2,p7)*zaba22(p5,p3,p4,p1,p2,p7)
     & *(0.25_dp*zb(p1,p4)*za(p1,p5)+0.25_dp*zb(p2,p4)*za(p2,p5)
     & +0.125_dp*zb(p3,p4)*za(p3,p5))
     & /(s(p1,p2)*s(p3,p4)*za(p5,p6)*ang7x12x56x7)

     & +0.125_dp*za(p1,p5)*zb(p4,p7)*zab2(p3,p1,p4,p2)
     & *zaba22(p5,p3,p4,p1,p2,p7)
     & /(s(p1,p2)*s(p3,p4)*za(p5,p6)*ang7x12x56x7)

     & +0.125_dp*za(p1,p5)*zb(p2,p7)*zb(p4,p6)*zab2(p7,p5,p6,p4)
     & /(zb(p3,p4)*s(p5,p6)*ang7x12x56x7)

     & +0.125_dp*za(p1,p3)*zb(p2,p7)*zb(p6,p7)**2*zab2(p3,p1,p2,p7)
     & /(s(p1,p2)*za(p3,p4)*zb(p5,p6)*sqr7x12x56x7)

      return
      end
