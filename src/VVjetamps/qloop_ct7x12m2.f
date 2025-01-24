!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function qloop_ct7x12m2(p1,p2,p3,p4,p5,p6,p7,za,zb)
c symmetric part under 12 <-> 56
c     This is an implementation of deLaurentis, Campbell and Ellis
c     arXiv:2203.17170, Eq(4.33)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer::p1,p2,p3,p4,p5,p6,p7
      complex(dp):: qloop_ct7x12m2
      real(dp)::s34,s56
      complex(dp)::zab2,zba2,zaba22,zbab22,za7x3456x7,zb7x3456x7

c     statement functions
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zba2(p1,p2,p3,p4)=zb(p1,p2)*za(p2,p4)+zb(p1,p3)*za(p3,p4)
      zaba22(p1,p2,p3,p4,p5,p6)=
     & za(p1,p2)*zba2(p2,p4,p5,p6)+za(p1,p3)*zba2(p3,p4,p5,p6)
      zbab22(p1,p2,p3,p4,p5,p6)=
     & zb(p1,p2)*zab2(p2,p4,p5,p6)+zb(p1,p3)*zab2(p3,p4,p5,p6)
c     end statement functions


      za7x3456x7=zaba22(p7,p3,p4,p5,p6,p7)
      zb7x3456x7=zbab22(p7,p3,p4,p5,p6,p7)
      s34=s(p3,p4)
      s56=s(p5,p6)

      qloop_ct7x12m2=
     & +za(p1,p7)**2*za(p3,p5)*zb(p4,p6)*zab2(p7,p1,p2,p7)
     & *(s(p3,p5)+s(p3,p6)+s(p4,p5)+s(p4,p6))
     & /(za(p1,p2)*za7x3456x7**2)
     & -zb(p2,p4)*zb(p2,p6)*za(p3,p5)*zab2(p7,p1,p2,p7)
     & /(zb(p1,p2)*za7x3456x7)
     & +za(p1,p3)*za(p1,p5)*zb(p4,p6)*zab2(p7,p1,p2,p7)
     & /(za(p1,p2)*za7x3456x7)
     & +za(p1,p7)**2*za(p3,p5)*zb(p4,p7)*zb(p6,p7)
     & /(za(p1,p2)*za7x3456x7)
     & -zb(p2,p7)**2*za(p3,p5)*zb(p4,p7)*zb(p6,p7)
     & /(zb(p1,p2)*zb7x3456x7)
      qloop_ct7x12m2=qloop_ct7x12m2/(s34*s56)

      return
      end
