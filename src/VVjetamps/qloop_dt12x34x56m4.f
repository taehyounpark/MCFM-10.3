!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function qloop_dt12x34x56m4(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
c     This is an implementation of deLaurentis, Campbell and Ellis
c     arXiv:2203.17170, Eq(4.13)
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: p1,p2,p3,p4,p5,p6,p7,p8
      complex(dp):: qloop_dt12x34x56m4,zab2,zaba22,zbab22,zbaba222

c--- statement functions
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zaba22(p1,p2,p3,p4,p5,p6)=
     & zab2(p1,p2,p3,p4)*za(p4,p6)+zab2(p1,p2,p3,p5)*za(p5,p6)
      zbab22(p1,p2,p3,p4,p5,p6)=
     & zb(p1,p2)*zab2(p2,p4,p5,p6)+zb(p1,p3)*zab2(p3,p4,p5,p6)
      zbaba222(p1,p2,p3,p4,p5,p6,p7,p8)=
     & +zbab22(p1,p2,p3,p4,p5,p6)*za(p6,p8)
     & +zbab22(p1,p2,p3,p4,p5,p7)*za(p7,p8)
c--- end statement functions

      qloop_dt12x34x56m4=2*za(p3,p5)*zb(p4,p6)
     & *(zbaba222(p2,p3,p4,p1,p2,p5,p6,p1)
     &  -zbaba222(p2,p5,p6,p1,p2,p3,p4,p1))
     & /(s(p1,p2)*s(p3,p4)*s(p5,p6)*zaba22(p7,p3,p4,p5,p6,p7))
      return
      end
