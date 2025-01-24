!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function qloop_dt56x12x34m2(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: p1,p2,p3,p4,p5,p6,p7
      complex(dp):: qloop_dt56x12x34m2,dt56x12x34m2_sym
      qloop_dt56x12x34m2=
     &  dt56x12x34m2_sym(p1,p2,p3,p4,p5,p6,p7,za,zb)
     & +dt56x12x34m2_sym(p1,p2,p5,p6,p3,p4,p7,za,zb)
      return
      end

      function dt56x12x34m2_sym(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
c     This is an implementation of deLaurentis, Campbell and Ellis
c     arXiv:2203.17170y, Eq(4.18)
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: p1,p2,p3,p4,p5,p6,p7,p8
      complex(dp):: zab2,zaba22,zbab22,zabab222,
     & dt56x12x34m2_sym,za7x3456x7,za7x3412x5,za7x561234x7
      real(dp):: s34,s56

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
      s34=s(p3,p4)
      s56=s(p5,p6)
      za7x3456x7=zaba22(p7,p3,p4,p5,p6,p7)
      za7x3412x5=zaba22(p7,p3,p4,p1,p2,p5)
      za7x561234x7=zabab222(p7,p5,p6,p1,p2,p3,p4,p7)
      dt56x12x34m2_sym=
     & za(p3,p5)*zb(p4,p6)*za7x561234x7
     & *(2*zab2(p7,p3,p4,p2)**2-zb(p1,p2)**2*za(p1,p7)**2)
     & /(zb(p1,p2)*za7x3456x7**2)

     & +2*za(p1,p3)**2*zb(p3,p4)*zb(p6,p7)*za7x3412x5
     & /(za(p1,p2)*za7x3456x7)

     & -2*za(p1,p3)*za(p5,p6)*zb(p4,p7)*zb(p6,p7)**2*zab2(p1,p3,p4,p7)
     & /(za(p1,p2)*zbab22(p7,p1,p2,p5,p6,p7))

     & +za(p3,p7)*zb(p4,p7)*zb(p6,p7)
     & *(za(p1,p3)*zb(p2,p3)*za(p5,p7)-za(p1,p5)*zab2(p7,p1,p5,p2))
     & /(za7x3456x7)

     & -2*za(p1,p3)*zb(p6,p7)
     & *(+zb(p2,p6)*zb(p3,p4)*za(p3,p7)*za(p5,p6)
     &   -za(p5,p7)*zb(p2,p4)*(s(p3,p4)+s(p3,p7)))
     & /(za7x3456x7)

     & +za(p3,p5)
     & *(-zb(p2,p4)*zb(p2,p6)*za7x561234x7
     &   -zb(p4,p7)*zb(p6,p7)*zab2(p7,p3,p4,p2)**2)
     & /(zb(p1,p2)*za7x3456x7)

      dt56x12x34m2_sym=dt56x12x34m2_sym/(4*s34*s56)
      return
      end
