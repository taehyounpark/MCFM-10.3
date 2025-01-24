!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function qloop_d12x34x56m4(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
c     This is an implementation of deLaurentis, Campbell and Ellis
c     arXiv:2203.17170, Eq(4.11)
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6,p7,p8
      complex(dp):: qloop_d12x34x56m4,zab2,zaba22,zbab22,zbaba222,
     & zba21,zba65,zba43,ang7x3456x7,sqr7x3456x7
c--- statement functions
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zaba22(p1,p2,p3,p4,p5,p6)=
     & zab2(p1,p2,p3,p4)*za(p4,p6)+zab2(p1,p2,p3,p5)*za(p5,p6)
      zbab22(p1,p2,p3,p4,p5,p6)=
     & zb(p1,p2)*zab2(p2,p4,p5,p6)+zb(p1,p3)*zab2(p3,p4,p5,p6)
      zbaba222(p1,p2,p3,p4,p5,p6,p7,p8)=
     & +zb(p1,p2)*zaba22(p2,p4,p5,p6,p7,p8)
     & +zb(p1,p3)*zaba22(p3,p4,p5,p6,p7,p8)
c--- end statement functions
      zba21=zbaba222(p2,p3,p4,p1,p2,p5,p6,p1)
     &     -zbaba222(p2,p5,p6,p1,p2,p3,p4,p1)
      zba43=zbaba222(p4,p3,p4,p1,p2,p5,p6,p3)
     &     -zbaba222(p4,p5,p6,p1,p2,p3,p4,p3)
      zba65=zbaba222(p6,p3,p4,p1,p2,p5,p6,p5)
     &     -zbaba222(p6,p5,p6,p1,p2,p3,p4,p5)
      ang7x3456x7=zaba22(p7,p3,p4,p5,p6,p7)
      sqr7x3456x7=zbab22(p7,p3,p4,p5,p6,p7)

      qloop_d12x34x56m4=zba21
     & /(s(p1,p2)*s(p3,p4)*s(p5,p6)*ang7x3456x7**2*sqr7x3456x7)
     & *(zba43*zba65-za(p3,p5)*zb(p4,p6)*ang7x3456x7*sqr7x3456x7)

      return
      end

