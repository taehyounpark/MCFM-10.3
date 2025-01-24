!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function qloop_diff(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
c     This is an implementation of deLaurentis, Campbell and Ellis
c     arXiv:2203.17170, Eq(4.16)
c     this function returns
c     qloop_d12x34x56m2(1,2,3,4,5,6,7)
c     - symmetric part of qloop_d34x56x12m2_sym(5,6,1,2,3,4,7)
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp):: qloop_diff,
     & diff_symbit

      qloop_diff=
     & +diff_symbit(p1,p2,p3,p4,p5,p6,p7,za,zb)
     & +diff_symbit(p1,p2,p5,p6,p3,p4,p7,za,zb)
      return
      end

      function diff_symbit(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer p1,p2,p3,p4,p5,p6,p7,p8
      complex(dp):: diff_symbit,zab2,zaba22,zbab22,zabab222,
     & ang7x34x56x7,sqr7x34x56x7,zabab

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
      ang7x34x56x7=zaba22(p7,p3,p4,p5,p6,p7)
      sqr7x34x56x7=zbab22(p7,p3,p4,p5,p6,p7)
      zabab=zabab222(p7,p3,p4,p1,p2,p5,p6,p7)
      diff_symbit=
     & +0.25_dp*za(p1,p3)*zb(p2,p4)*za(p5,p7)*zaba22(p5,p1,p2,p3,p4,p7)
     & *zabab
     & /(s(p1,p2)*s(p3,p4)*za(p5,p6)*ang7x34x56x7**2)

     & +0.5_dp*za(p1,p2)*zb(p2,p4)*za(p3,p5)*za(p5,p7)*zab2(p7,p3,p4,p2)
     & *zabab
     & /(s(p1,p2)*s(p3,p4)*za(p5,p6)*ang7x34x56x7**2)

     & -0.25_dp*za(p3,p5)*zb(p4,p7)*zab2(p1,p3,p4,p2)
     & *zaba22(p5,p1,p2,p3,p4,p7)
     & /(s(p1,p2)*s(p3,p4)*za(p5,p6)*ang7x34x56x7)

     & -3/8._dp*za(p1,p3)*zb(p2,p7)*za(p3,p5)*zaba22(p5,p1,p2,p3,p4,p7)
     & /(s(p1,p2)*za(p3,p4)*za(p5,p6)*ang7x34x56x7)

     & -1/8._dp*za(p1,p3)*za(p1,p5)*zb(p4,p7)*zaba22(p5,p1,p2,p3,p4,p7)
     & /(za(p1,p2)*s(p3,p4)*za(p5,p6)*ang7x34x56x7)

     & -3/8._dp*zb(p2,p4)*zb(p2,p7)*za(p3,p5)*zaba22(p5,p1,p2,p3,p4,p7)
     & /(zb(p1,p2)*s(p3,p4)*za(p5,p6)*ang7x34x56x7)

     & -3/8._dp*zb(p2,p4)*za(p3,p5)*zb(p6,p7)*zab2(p7,p3,p4,p2)
     & /(zb(p1,p2)*s(p3,p4)*ang7x34x56x7)

     & +1/16._dp*za(p1,p3)*za(p1,p5)*zb(p4,p7)*zb(p6,p7)
     & /(za(p1,p2)*s(p3,p4)*s(p5,p6))

     & +1/16._dp*za(p3,p5)*zb(p4,p7)*zb(p6,p7)
     & *zab2(p1,p3,p4,p7)*zab2(p1,p5,p6,p7)
     & /(za(p1,p2)*s(p3,p4)*s(p5,p6)*sqr7x34x56x7)
      return
      end
