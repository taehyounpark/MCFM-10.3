!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function qloop_dt12x34x56m2(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
c     This is an implementation of deLaurentis, Campbell and Ellis
c     arXiv:2203.17170, Eq(4.12)
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer p1,p2,p3,p4,p5,p6,p7,p8
      complex(dp):: qloop_dt12x34x56m2,zab2,zaba22,zbab22,zabab222,
     & ang7x3456x7,sqr7x3456x7
      real(dp):: s12,s34,s56
c---  statement functions
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zaba22(p1,p2,p3,p4,p5,p6)=
     & zab2(p1,p2,p3,p4)*za(p4,p6)+zab2(p1,p2,p3,p5)*za(p5,p6)
      zbab22(p1,p2,p3,p4,p5,p6)=
     & zb(p1,p2)*zab2(p2,p4,p5,p6)+zb(p1,p3)*zab2(p3,p4,p5,p6)
      zabab222(p1,p2,p3,p4,p5,p6,p7,p8)=
     & +zaba22(p1,p2,p3,p4,p5,p6)*zb(p6,p8)
     & +zaba22(p1,p2,p3,p4,p5,p7)*zb(p7,p8)
c---  end statement functions
      s12=s(p1,p2)
      s34=s(p3,p4)
      s56=s(p5,p6)
      ang7x3456x7=zaba22(p7,p3,p4,p5,p6,p7)
      sqr7x3456x7=zbab22(p7,p3,p4,p5,p6,p7)
      qloop_dt12x34x56m2=
     & -2*za(p1,p7)*za(p3,p5)*zb(p1,p2)*zb(p4,p6)
     & *zaba22(p1,p3,p4,p5,p6,p7)*zabab222(p7,p1,p2,p3,p4,p5,p6,p7)
     & /ang7x3456x7**2

     & +za(p1,p2)*za(p3,p5)*zb(p2,p6)*zb(p2,p7)
     & *zabab222(p7,p1,p2,p3,p4,p5,p6,p4)/ang7x3456x7

     & -za(p1,p3)*za(p1,p5)*zb(p1,p2)*zb(p4,p6)
     & *zabab222(p7,p1,p2,p3,p4,p5,p6,p7)/ang7x3456x7

     & +za(p3,p5)*zb(p1,p2)*zb(p4,p7)*zab2(p1,p2,p7,p6)
     & *zaba22(p1,p3,p4,p5,p6,p7)/ang7x3456x7

     & -za(p1,p2)*za(p3,p5)*zb(p3,p4)*zb(p2,p7)**2*zb(p6,p7)
     & *zab2(p3,p5,p6,p7)/sqr7x3456x7

     & -za(p1,p2)*za(p3,p5)*zb(p2,p7)**2*zb(p4,p6)
      qloop_dt12x34x56m2=qloop_dt12x34x56m2/(2*s12*s34*s56)
      return
      end

